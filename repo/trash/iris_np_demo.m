function iris_np_demo()
% IRIS_NP_DEMO  IRIS-NP: convex collision-free regions in N-D configuration space.
%
% Provenance: Petersen & Tedrake, "Growing Convex Collision-Free Regions in
% Configuration Space using Nonlinear Programming", arXiv:2303.14737.
%   Sec. II      : problem (1), counterexample search (2), hyperplanes (3),(4),
%                  Algorithms 1 and 2, inscribed ellipse (5), termination
%   Sec. III-A   : collision-pair ordering by task-space distance at the seed
%   Sec. III-B   : probabilistic certification, hit-and-run restarts
%   Sec. III-C   : additional nonlinear constraints g(q) <= 0, Eq. (6)
%   Sec. III-D   : multi-region coverage (previous regions as C-space obstacles)
%   Sec. III-E   : region refinement for runtime obstacles
%
% Solver substitutions (paper -> here):
%   SNOPT (counterexample NLP)  -> fmincon, SQP algorithm
%   Mosek/Gurobi (ellipse SDP)  -> fmincon on a lower-triangular factor L with
%       analytic gradients. Problem: max sum(log(diag(L))) s.t.
%       ||L' a_i||_2 + a_i' d <= b_i. This is convex in (d,L), so the local
%       solver returns the global optimum. Ellipse: E = {d + L u : ||u|| <= 1},
%       metric C'C = inv(L*L') (paper's C equals inv(L)). Volume ~ det(L).
%
% Requirements: Optimization Toolbox (fmincon, linprog, quadprog).
% Geometry model (this file): robot links are line segments in the workspace,
% obstacles are convex polytopes {x : A x <= b} in the workspace. Link
% thickness is handled by inflating obstacles. Any forward-kinematics map
% can be used by supplying link.fk(q) -> [P0 P1] (m x 2 segment endpoints).
% Forward-kinematics gradients are obtained by fmincon finite differences.
%
% Demo: planar 4-link arm, 3 box obstacles, self-collision pairs, and a
% nonlinear end-effector orientation constraint |sum(q)| <= 1 rad.

warning('on','all');
rng(1);

%% ---------------- Scenario ----------------
linkLen = [1 1 1 1];
n       = numel(linkLen);
links   = makePlanarArm(linkLen, [0;0]);
obs     = {boxObstacle([1.0 1.6],[3.0 2.4]), ...
           boxObstacle([2.0 -1.5],[2.4 -0.2]), ...
           boxObstacle([-0.5 1.2],[0.3 2.0])};
geoms   = [links, obs];
pairs   = buildPairs(geoms);

% Extra constraints g(q) <= 0 (Sec. III-C): end-effector orientation bound
pairs{end+1} = struct('kind','extra','g',@(q) sum(q)-1.0);
pairs{end+1} = struct('kind','extra','g',@(q) -sum(q)-1.0);

q0 = [0.2; 0.1; -0.1; 0.0];
ql = -0.9*pi*ones(n,1);
qu =  0.9*pi*ones(n,1);

opts = struct('maxInfeasible',3,'maxIter',5,'tolGrowth',0.02,'delta',0.01);

%% ---------------- Region growing ----------------
t0 = tic;
region = irisNP(q0, ql, qu, geoms, pairs, opts);
fprintf('Region 1: %d faces, log-volume(det L)=%.4f, iterations=%d, time=%.1f s\n', ...
    size(region.A,1), region.logVol(end), region.nIter, toc(t0));

%% ---------------- Verification ----------------
verifyRegion(region, q0, geoms, pairs, 500);

%% ---------------- Second region (coverage, Sec. III-D) ----------------
q1 = [1.2; -0.8; 0.9; -0.5];
if all(isCollisionAny(pairs, geoms, q1, region.opts)==false) && all(q1 > ql) && all(q1 < qu)
    opts2 = opts; opts2.prevRegions = {region};
    region2 = irisNP(q1, ql, qu, geoms, pairs, opts2);
    fprintf('Region 2 (non-overlap heuristic): %d faces, log-volume=%.4f\n', ...
        size(region2.A,1), region2.logVol(end));
    verifyRegion(region2, q1, geoms, pairs, 300);
else
    warning('Second seed q1 is in collision or out of limits; coverage demo skipped.');
end

%% ---------------- Refinement for a new obstacle (Sec. III-E) ----------------
newObs   = boxObstacle([0.6 -0.9],[1.0 -0.4]);
geoms2   = [geoms, {newObs}];
idxNew   = numel(geoms2);
newPairs = {};
for k = 1:numel(links)
    newPairs{end+1} = struct('kind','linkpoly','i',k,'j',idxNew); %#ok<AGROW>
end
allPairs2 = [pairs, newPairs];
tR = tic;
refined = irisNPRefine(region, q0, geoms2, newPairs, opts);
fprintf('Refined region: %d faces (was %d), time=%.1f s\n', ...
    size(refined.A,1), size(region.A,1), toc(tR));
verifyRegion(refined, q0, geoms2, allPairs2, 500);

%% ---------------- Plot ----------------
plotDemo(region, links, obs, q0, linkLen);
end


%% =====================================================================
%  CORE: Algorithm 1 (IRIS-NP)
%% =====================================================================

function region = irisNP(q0, ql, qu, geoms, pairs, opts)
% Algorithm 1. Returns polytope {q : A q <= b}, inscribed ellipse (L,d),
% and diagnostics. Seed q0 must be collision free.
opts = setDefaults(opts);
q0 = q0(:);
n  = numel(q0);

if isempty(opts.A0)
    ql = ql(:); qu = qu(:);
    if any(q0 < ql) || any(q0 > qu)
        error('irisNP:seedOutsideLimits','Seed lies outside the joint limits.');
    end
    A0 = [eye(n); -eye(n)];
    b0 = [qu; -ql];
else
    A0 = opts.A0; b0 = opts.b0(:);
    if any(A0*q0 > b0 + 1e-9)
        error('irisNP:seedOutsidePolytope','Seed lies outside the supplied initial polytope.');
    end
end

for k = 1:numel(pairs)
    if isCollision(pairs{k}, geoms, q0, opts)
        error('irisNP:seedInCollision','Seed violates collision/constraint pair %d.', k);
    end
end

% Sec. III-A: order pairs by task-space distance at the seed
if opts.sortPairs
    keys = zeros(numel(pairs),1);
    for k = 1:numel(pairs)
        keys(k) = pairDistance(pairs{k}, geoms, q0, opts);
    end
    [~,ord] = sort(keys,'ascend');
    pairs = pairs(ord);
end

% Initial ellipse: hypersphere of radius epsilon about the seed
L = opts.epsilon*eye(n);
d = q0;
A = A0; b = b0;
logVol = n*log(opts.epsilon);
facesHist = size(A0,1);
nSolvesTotal = 0;
converged = false;
nIter = 0;

for it = 1:opts.maxIter
    M = inv(L*L'); M = (M+M')/2;

    [Anew,bnew,info] = addSeparatingHyperplanes(A0,b0,d,M,pairs,geoms,opts);
    nSolvesTotal = nSolvesTotal + info.nSolves;
    if info.centreCut
        warning('irisNP:centreCut',['Iteration %d: a hyperplane (margin delta=%g) excludes the ' ...
            'ellipse centre. Reduce delta or epsilon. Returning previous iterate.'], it, opts.delta);
        break;
    end
    if opts.requireSeed && any(Anew*q0 > bnew + 1e-9)
        warning('irisNP:seedExcluded','Iteration %d: seed excluded by new hyperplanes. Stopping (seed-containment criterion).', it);
        break;
    end

    [Ln,dn,ok] = inscribedEllipse(Anew,bnew,L,d,opts);
    if ~ok
        warning('irisNP:ellipseFailed','Iteration %d: inscribed-ellipse step failed. Returning previous iterate.', it);
        break;
    end

    newLogVol = sum(log(diag(Ln)));
    if newLogVol < logVol - 1e-6
        warning('irisNP:nonMonotone','Iteration %d: ellipse volume decreased (%.6g -> %.6g).', it, logVol, newLogVol);
    end
    growth = exp(newLogVol - logVol) - 1;

    A = Anew; b = bnew; L = Ln; d = dn; logVol(end+1) = newLogVol; %#ok<AGROW>
    facesHist(end+1) = size(A,1); %#ok<AGROW>
    nIter = it;
    if opts.verbose
        fprintf('  iter %d: faces=%d, solves=%d, volume growth=%.2f%%\n', it, size(A,1), info.nSolves, 100*growth);
    end
    if growth < opts.tolGrowth
        converged = true;
        break;
    end
end

if nIter == 0
    error('irisNP:noRegion','No valid region was produced (first iteration failed).');
end
if ~converged
    warning('irisNP:notConverged', ...
        'Terminated by iteration limit/failure before volume growth fell below %.3g.', opts.tolGrowth);
end

if opts.prune
    [A,b] = pruneRedundant(A,b,1e-9);
end

region = struct('A',A,'b',b,'L',L,'d',d,'M',inv(L*L'),'seed',q0, ...
    'logVol',logVol,'faces',facesHist,'nIter',nIter,'converged',converged, ...
    'nSolves',nSolvesTotal,'seedContained',all(A*q0 <= b+1e-9),'opts',opts);
end


%% =====================================================================
%  Algorithm 2: AddSeparatingHyperplanes
%% =====================================================================

function [A,b,info] = addSeparatingHyperplanes(A,b,d,M,pairs,geoms,opts)
info = struct('centreCut',false,'nSolves',0,'nCuts',0);

% Sec. III-D: previously generated regions act as C-space obstacles
for r = 1:numel(opts.prevRegions)
    [a,bj,ok] = regionSeparator(opts.prevRegions{r}, d, M, opts);
    if ok
        A = [A; a]; b = [b; bj]; %#ok<AGROW>
        if a*d > bj - 1e-12
            info.centreCut = true; return;
        end
    end
end

for k = 1:numel(pairs)
    pair   = pairs{k};
    fails  = 0;
    cuts   = 0;
    walker = d;
    while fails < opts.maxInfeasible
        qInit = hitAndRun(A,b,walker,opts.mixSteps);
        [qs,ok] = solveCounterexample(pair, geoms, qInit, A, b, M, d, opts);
        info.nSolves = info.nSolves + 1;
        if ok
            v  = M*(qs-d);
            nv = norm(v);
            if nv < 1e-12
                error('addSeparatingHyperplanes:degenerate', ...
                    'Counterexample coincides with ellipse centre (pair %d); seed is in collision.', k);
            end
            a  = (v/nv)';
            bj = a*qs - opts.delta;               % Eq. (4)
            A  = [A; a]; b = [b; bj]; %#ok<AGROW>
            cuts = cuts + 1; info.nCuts = info.nCuts + 1;
            fails = 0;
            if a*d > bj - 1e-12
                info.centreCut = true; return;
            end
            walker = d;
            if cuts >= opts.maxCutsPerPair
                warning('addSeparatingHyperplanes:maxCuts', ...
                    'Pair %d reached the %d-hyperplane cap; its separation may be incomplete.', k, opts.maxCutsPerPair);
                break;
            end
        else
            fails = fails + 1;
            walker = qInit;                        % continue the Markov chain
        end
    end
end
end


%% =====================================================================
%  Counterexample search, Eq. (2) and Eq. (6)
%% =====================================================================

function [qs,ok] = solveCounterexample(pair, geoms, qInit, A, b, M, d, opts)
n  = numel(d);
qs = []; ok = false;
N  = size(A,1);
obj = @(z) quadCost(z, d, M, n);

switch pair.kind
    case 'extra'                                   % g(q) >= violTol, Aq <= b
        nl  = @(z) deal(opts.violTol - pair.g(z(1:n)), []);
        z0  = qInit;
        [z,~,ef] = fmincon(obj, z0, A, b, [], [], [], [], nl, opts.nlpOpts);
        if ef > 0 && all(A*z <= b + 1e-7) && pair.g(z(1:n)) >= opts.violTol - 1e-8
            qs = z(1:n); ok = true;
        end

    case 'linkpoly'                                % FK(q,t) = x, x in obstacle
        gl = geoms{pair.i}; go = geoms{pair.j};
        m  = size(go.A,2);
        Ai = [A, zeros(N,1+m); zeros(size(go.A,1),n+1), go.A];
        bi = [b; go.b];
        lb = [-inf(n,1); 0; -inf(m,1)];
        ub = [ inf(n,1); 1;  inf(m,1)];
        z0 = [qInit; 0.5; go.xc];
        nl = @(z) linkPolyCon(z, n, gl.fk);
        [z,~,ef] = fmincon(obj, z0, Ai, bi, [], [], lb, ub, nl, opts.nlpOpts);
        if ef > 0 && all(Ai*z <= bi + 1e-7) && norm(linkPolyEq(z,n,gl.fk),inf) < opts.eqTol
            qs = z(1:n); ok = true;
        end

    case 'linklink'                                % FK_i(q,t) = FK_j(q,s)
        gi = geoms{pair.i}; gj = geoms{pair.j};
        Ai = [A, zeros(N,2)];
        lb = [-inf(n,1); 0; 0];
        ub = [ inf(n,1); 1; 1];
        z0 = [qInit; 0.5; 0.5];
        nl = @(z) linkLinkCon(z, n, gi.fk, gj.fk);
        [z,~,ef] = fmincon(obj, z0, Ai, b, [], [], lb, ub, nl, opts.nlpOpts);
        if ef > 0 && all(Ai*z <= b + 1e-7) && norm(linkLinkEq(z,n,gi.fk,gj.fk),inf) < opts.eqTol
            qs = z(1:n); ok = true;
        end

    otherwise
        error('solveCounterexample:kind','Unknown pair kind "%s".', pair.kind);
end
end

function [f,g] = quadCost(z, d, M, n)
e = z(1:n) - d;
f = e'*M*e;
g = zeros(size(z));
g(1:n) = 2*M*e;
end

function [c,ceq] = linkPolyCon(z,n,fk)
c = []; ceq = linkPolyEq(z,n,fk);
end
function r = linkPolyEq(z,n,fk)
P = fk(z(1:n)); t = z(n+1); x = z(n+2:end);
r = P(:,1) + t*(P(:,2)-P(:,1)) - x(:);
end
function [c,ceq] = linkLinkCon(z,n,fka,fkb)
c = []; ceq = linkLinkEq(z,n,fka,fkb);
end
function r = linkLinkEq(z,n,fka,fkb)
Pa = fka(z(1:n)); Pb = fkb(z(1:n)); t = z(n+1); s = z(n+2);
r = (Pa(:,1)+t*(Pa(:,2)-Pa(:,1))) - (Pb(:,1)+s*(Pb(:,2)-Pb(:,1)));
end


%% =====================================================================
%  Inscribed ellipse, Eq. (5)
%% =====================================================================

function [L,d,ok] = inscribedEllipse(A,b,L0,d0,opts)
[N,n] = size(A);
L = L0; d = d0; ok = false;

% Chebyshev centre (strictly feasible starting point)
nrm = sqrt(sum(A.^2,2));
[xr,~,ef] = linprog([zeros(n,1);-1], [A nrm], b, [], [], [-inf(n,1);0], [], opts.lpOpts);
if ef <= 0 || xr(end) < 1e-9
    return;                                         % empty interior
end
dC = xr(1:n); rC = xr(end);
cands = {rC*eye(n), dC};
if max(sqrt(sum((A*L0).^2,2)) + A*d0 - b) <= 1e-9
    if sum(log(diag(L0))) > n*log(rC)
        cands = {L0, d0};                           % warm start keeps monotonicity
    end
end
Ls = cands{1}; ds = cands{2};

mask = tril(true(n)); nl = nnz(mask);
[rr,cc] = find(mask);
diagIdx = find(rr == cc);
x0 = [ds; Ls(mask)];
lb = -inf(n+nl,1); lb(n+diagIdx) = 1e-12;
obj = @(x) ellObj(x, n, diagIdx);
nlc = @(x) ellCon(x, A, b, n, mask, rr, cc, N);
[x,~,ef2] = fmincon(obj, x0, [], [], [], [], lb, [], nlc, opts.ellOpts);

Lk = zeros(n); Lk(mask) = x(n+1:end); dk = x(1:n);
viol = max(sqrt(sum((A*Lk).^2,2)) + A*dk - b);
if ef2 <= 0
    warning('inscribedEllipse:exitflag','Ellipse solver exit flag %d (max violation %.2e).', ef2, viol);
end
if viol <= 1e-8 && all(diag(Lk) > 0)
    L = Lk; d = dk; ok = true;
elseif ef2 <= 0
    L = Ls; d = ds; ok = true;                      % feasible Chebyshev/warm fallback
    warning('inscribedEllipse:fallback','Using feasible fallback ellipse (not maximal).');
end
end

function [f,g] = ellObj(x,n,diagIdx)
dv = x(n+diagIdx);
f = -sum(log(dv));
g = zeros(size(x));
g(n+diagIdx) = -1./dv;
end

function [c,ceq,gc,gceq] = ellCon(x,A,b,n,mask,rr,cc,N)
d = x(1:n);
L = zeros(n); L(mask) = x(n+1:end);
W  = A*L;                                           % row i = (L' a_i)'
nW = max(sqrt(sum(W.^2,2)),1e-14);
c  = nW + A*d - b;
ceq = []; gceq = [];
if nargout > 2
    nl = numel(rr);
    gc = zeros(n+nl,N);
    gc(1:n,:) = A';
    for k = 1:nl
        gc(n+k,:) = (A(:,rr(k)).*W(:,cc(k))./nW)';
    end
end
end


%% =====================================================================
%  Sampling, redundancy pruning, multi-region separator
%% =====================================================================

function x = hitAndRun(A,b,x0,nSteps)
% Hit-and-run MCMC inside {A x <= b} (Belisle, Boneh, Caron 1998).
x = x0(:); n = numel(x);
for s = 1:nSteps
    u = randn(n,1); u = u/norm(u);
    Au = A*u;
    slack = max(b - A*x, 0);
    pos = Au >  1e-12;
    neg = Au < -1e-12;
    if ~any(pos) || ~any(neg)
        error('hitAndRun:unbounded','Polytope is unbounded along the sampled direction.');
    end
    tmax = min(slack(pos)./Au(pos));
    tmin = max(slack(neg)./Au(neg));
    x = x + (tmin + (tmax-tmin)*rand)*u;
end
end

function [A,b] = pruneRedundant(A,b,tol)
lpOpts = optimoptions('linprog','Display','none');
keep = true(size(A,1),1);
for i = 1:size(A,1)
    idx = keep; idx(i) = false;
    [~,fval,ef] = linprog(-A(i,:)', A(idx,:), b(idx), [], [], [], [], lpOpts);
    if ef == 1 && -fval <= b(i) + tol
        keep(i) = false;
    end
end
A = A(keep,:); b = b(keep);
end

function [a,bj,ok] = regionSeparator(reg,d,M,opts)
% Tangent hyperplane separating the ellipse metric from a previous region
% (Sec. III-D; Deits & Tedrake closest-point construction in C-space).
a = []; bj = []; ok = false;
if all(reg.A*d <= reg.b + 1e-9)
    warning('regionSeparator:inside','Ellipse centre lies inside a previous region; no separator added.');
    return;
end
[xs,~,ef] = quadprog(2*M, -2*M*d, reg.A, reg.b, [], [], [], [], [], opts.qpOpts);
if ef <= 0
    warning('regionSeparator:qp','Closest-point QP failed (exit flag %d); separator skipped.', ef);
    return;
end
v = M*(xs-d); nv = norm(v);
if nv < 1e-12, return; end
a = (v/nv)'; bj = a*xs; ok = true;
end


%% =====================================================================
%  Refinement for runtime obstacles (Sec. III-E)
%% =====================================================================

function reg2 = irisNPRefine(region, q0, geoms, changedPairs, opts)
% One pass of hyperplane addition with the final ellipse as metric,
% keeping the original hyperplanes; only changed pairs are examined.
opts = setDefaults(opts);
for k = 1:numel(changedPairs)
    if isCollision(changedPairs{k}, geoms, q0, opts)
        error('irisNPRefine:seedInCollision','Seed collides with changed geometry; refinement undefined.');
    end
end
[A,b,info] = addSeparatingHyperplanes(region.A, region.b, region.d, region.M, changedPairs, geoms, opts);
if info.centreCut
    warning('irisNPRefine:centreCut','New geometry intersects the ellipse centre neighbourhood; regrowing from seed inside the original region.');
    o2 = opts; o2.A0 = region.A; o2.b0 = region.b;
    reg2 = irisNP(q0, [], [], geoms, changedPairs, o2);
    return;
end
if any(A*q0 > b + 1e-9)
    warning('irisNPRefine:seedExcluded', ...
        'Seed excluded by refinement; growing a new region from the seed inside the original region (changed pairs only).');
    o2 = opts; o2.A0 = region.A; o2.b0 = region.b;
    reg2 = irisNP(q0, [], [], geoms, changedPairs, o2);
    return;
end
[L,d,ok] = inscribedEllipse(A,b,region.L,region.d,opts);
if ~ok
    warning('irisNPRefine:ellipse','Ellipse recomputation failed; keeping previous ellipse metadata.');
    L = region.L; d = region.d;
end
if opts.prune, [A,b] = pruneRedundant(A,b,1e-9); end
reg2 = region;
reg2.A = A; reg2.b = b; reg2.L = L; reg2.d = d; reg2.M = inv(L*L');
reg2.logVol = [region.logVol, sum(log(diag(L)))];
reg2.faces = [region.faces, size(A,1)];
reg2.nIter = region.nIter + 1;
reg2.nSolves = region.nSolves + info.nSolves;
reg2.seedContained = all(A*q0 <= b+1e-9);
reg2.opts = opts;
end


%% =====================================================================
%  Collision queries and task-space distances
%% =====================================================================

function tf = isCollision(pair, geoms, q, opts)
if strcmp(pair.kind,'extra')
    tf = pair.g(q(:)) > opts.violTol;
else
    tf = pairDistance(pair, geoms, q, opts) <= opts.collTol;
end
end

function tf = isCollisionAny(pairs, geoms, q, opts)
opts = setDefaults(opts);
tf = false(numel(pairs),1);
for k = 1:numel(pairs)
    tf(k) = isCollision(pairs{k}, geoms, q, opts);
end
end

function dist = pairDistance(pair, geoms, q, opts)
% Exact minimum distance between two convex task-space geometries (convex QP).
switch pair.kind
    case 'extra'
        dist = max(-pair.g(q(:)),0);                % surrogate (feasibility margin)
    case 'linkpoly'
        P = geoms{pair.i}.fk(q(:)); go = geoms{pair.j};
        m = size(go.A,2);
        D = P(:,2)-P(:,1);
        R = [D, -eye(m)];
        [z,~,ef] = quadprog(2*(R'*R), 2*R'*P(:,1), [zeros(size(go.A,1),1), go.A], go.b, ...
            [], [], [0;-inf(m,1)], [1;inf(m,1)], [], opts.qpOpts);
        if ef <= 0, error('pairDistance:qp','Distance QP failed (exit flag %d).', ef); end
        dist = norm(R*z + P(:,1));
    case 'linklink'
        Pa = geoms{pair.i}.fk(q(:)); Pb = geoms{pair.j}.fk(q(:));
        R = [Pa(:,2)-Pa(:,1), -(Pb(:,2)-Pb(:,1))];
        r0 = Pa(:,1)-Pb(:,1);
        [z,~,ef] = quadprog(2*(R'*R), 2*R'*r0, [], [], [], [], [0;0], [1;1], [], opts.qpOpts);
        if ef <= 0, error('pairDistance:qp','Distance QP failed (exit flag %d).', ef); end
        dist = norm(R*z + r0);
    otherwise
        error('pairDistance:kind','Unknown pair kind "%s".', pair.kind);
end
end


%% =====================================================================
%  Verification: Sec. III-B probabilistic certification + structural checks
%% =====================================================================

function verifyRegion(region, seed, geoms, pairs, nSamples)
opts = setDefaults(region.opts);
A = region.A; b = region.b; L = region.L; d = region.d;

% 1. Ellipse containment: max_{||u||<=1} a_i'(d+Lu) = a_i'd + ||L'a_i|| <= b_i
slackE = b - (A*d + sqrt(sum((A*L).^2,2)));
if min(slackE) < -1e-7
    warning('verifyRegion:ellipse','Ellipse not contained in polytope (worst slack %.3e).', min(slackE));
end
% 2. Seed containment
if ~all(A*seed <= b + 1e-9)
    warning('verifyRegion:seed','Seed is not contained in the region.');
end
% 3. Volume monotonicity
if any(diff(region.logVol) < -1e-6)
    warning('verifyRegion:monotone','log-volume history is not monotone: %s', mat2str(region.logVol,5));
end
% 4. Probabilistic certification: collision/violation fraction among
%    hit-and-run samples from the region (approximately uniform).
walker = d; nBad = 0; nPair = zeros(numel(pairs),1);
for s = 1:nSamples
    walker = hitAndRun(A,b,walker,opts.mixSteps);
    bad = false;
    for k = 1:numel(pairs)
        if isCollision(pairs{k}, geoms, walker, opts)
            nPair(k) = nPair(k) + 1; bad = true;
        end
    end
    nBad = nBad + bad;
end
fprintf('  Certification: faces=%d, samples=%d, colliding/violating=%d (%.3f%%), converged=%d\n', ...
    size(A,1), nSamples, nBad, 100*nBad/nSamples, region.converged);
if nBad > 0
    [~,kmax] = max(nPair);
    warning('verifyRegion:collisions', ...
        ['%d/%d samples violate a pair (worst pair index %d: %d). IRIS-NP certification is ' ...
         'probabilistic; raise opts.maxInfeasible or lower the margin delta.'], nBad, nSamples, kmax, nPair(kmax));
end
end


%% =====================================================================
%  Scenario helpers
%% =====================================================================

function links = makePlanarArm(len, base)
% Serial planar arm; link k is the segment between joint k and joint k+1.
n = numel(len);
links = cell(1,n);
for k = 1:n
    links{k} = struct('kind','link','fk',@(q) planarLinkEnds(q,len,base,k));
end
end

function P = planarLinkEnds(q,len,base,k)
ang = cumsum(q(:)');
pts = base(:) + cumsum([len(:)'.*cos(ang); len(:)'.*sin(ang)], 2);
J   = [base(:), pts];
P   = J(:,[k k+1]);
end

function o = boxObstacle(lo,hi)
lo = lo(:); hi = hi(:); m = numel(lo);
o = struct('kind','poly','A',[eye(m);-eye(m)],'b',[hi;-lo],'xc',(lo+hi)/2);
end

function pairs = buildPairs(geoms)
% link-polytope pairs for every link/obstacle; link-link pairs for
% non-adjacent links (adjacent links share a joint by construction).
isLink = cellfun(@(g) strcmp(g.kind,'link'), geoms);
iL = find(isLink); iO = find(~isLink);
pairs = {};
for a = iL
    for o = iO
        pairs{end+1} = struct('kind','linkpoly','i',a,'j',o); %#ok<AGROW>
    end
end
for p = 1:numel(iL)
    for r = p+2:numel(iL)
        pairs{end+1} = struct('kind','linklink','i',iL(p),'j',iL(r)); %#ok<AGROW>
    end
end
end

function opts = setDefaults(opts)
def = struct('epsilon',1e-2,'delta',1e-2,'maxIter',5,'tolGrowth',0.02, ...
    'maxInfeasible',3,'mixSteps',25,'maxCutsPerPair',200,'requireSeed',true, ...
    'sortPairs',true,'prune',true,'violTol',1e-6,'eqTol',1e-6,'collTol',1e-6, ...
    'prevRegions',{{}},'A0',[],'b0',[],'verbose',true);
f = fieldnames(def);
for k = 1:numel(f)
    if ~isfield(opts,f{k}), opts.(f{k}) = def.(f{k}); end
end
if ~isfield(opts,'nlpOpts') || isempty(opts.nlpOpts)
    opts.nlpOpts = optimoptions('fmincon','Algorithm','sqp','Display','none', ...
        'SpecifyObjectiveGradient',true,'ConstraintTolerance',1e-8, ...
        'OptimalityTolerance',1e-8,'StepTolerance',1e-10,'MaxIterations',300);
end
if ~isfield(opts,'ellOpts') || isempty(opts.ellOpts)
    opts.ellOpts = optimoptions('fmincon','Algorithm','sqp','Display','none', ...
        'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true, ...
        'ConstraintTolerance',1e-9,'OptimalityTolerance',1e-9,'StepTolerance',1e-12, ...
        'MaxIterations',500);
end
if ~isfield(opts,'lpOpts') || isempty(opts.lpOpts)
    opts.lpOpts = optimoptions('linprog','Display','none');
end
if ~isfield(opts,'qpOpts') || isempty(opts.qpOpts)
    opts.qpOpts = optimoptions('quadprog','Display','none');
end
end


%% =====================================================================
%  Plot (planar demo)
%% =====================================================================

function plotDemo(region, links, obs, q0, linkLen)
figure('Name','IRIS-NP planar 4-link arm');
subplot(1,2,1); hold on; axis equal; grid on;
for k = 1:numel(obs)
    lo = -obs{k}.b(numel(obs{k}.b)/2+1:end); hi = obs{k}.b(1:numel(obs{k}.b)/2);
    fill([lo(1) hi(1) hi(1) lo(1)],[lo(2) lo(2) hi(2) hi(2)],[0.3 0.3 0.3]);
end
walker = region.d;
cm = lines(25);
for s = 1:25
    walker = hitAndRun(region.A, region.b, walker, 40);
    pts = zeros(2,numel(links)+1);
    for k = 1:numel(links)
        P = links{k}.fk(walker);
        if k == 1, pts(:,1) = P(:,1); end
        pts(:,k+1) = P(:,2);
    end
    plot(pts(1,:),pts(2,:),'-o','Color',[cm(s,:) 0.5],'MarkerSize',3);
end
pts = zeros(2,numel(links)+1);
for k = 1:numel(links)
    P = links{k}.fk(q0);
    if k == 1, pts(:,1) = P(:,1); end
    pts(:,k+1) = P(:,2);
end
plot(pts(1,:),pts(2,:),'k-o','LineWidth',2.5,'MarkerFaceColor','k');
xlim([-sum(linkLen) sum(linkLen)]); ylim([-sum(linkLen) sum(linkLen)]);
title('Configurations sampled from region (black: seed)');
xlabel('x'); ylabel('y');

subplot(1,2,2);
plot(0:numel(region.logVol)-1, region.logVol,'-o','LineWidth',1.5); grid on;
xlabel('IRIS-NP iteration'); ylabel('log det L (ellipse log-volume + const)');
title(sprintf('Convergence, final faces = %d', size(region.A,1)));
end
