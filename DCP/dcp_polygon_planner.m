function plan = dcp_polygon_planner(x0, xf, obs, par)
%DCP_POLYGON_PLANNER  Collision-free waypoint trajectory among multiple convex polygons.
%   plan = dcp_polygon_planner(x0, xf, obs)
%   plan = dcp_polygon_planner(x0, xf, obs, par)
%
%   x0, xf : 2x1 start / goal
%   obs    : cell array, each cell an Nx2 vertex matrix of a CONVEX polygon (any vertex order)
%   par    : optional struct, missing fields take the defaults below
%       N 60 segments | Ts 0.1 sample time | vmax 2 max speed (Inf disables) | margin 0.15 clearance
%       wAcc 0 second-difference weight | tau0 10, mu 2, tauMax 1e6 penalty schedule
%       maxIt 40 | slackTol 1e-6 | objTol 1e-5
%
%   plan fields: x (2xN+1) position, t (1xN+1), v (2xN+1) feedforward velocity, L path length,
%                feasible, collisions, iters, time, refFun(tq) -> 2xM linear interpolation, par.
%
%   Requires CVX (Grant & Boyd).
%
%   Method (convex-concave procedure, Boyd's DCP view):
%   Avoidance of a convex polygon O = {x : A x <= b} is the disjunction exists k: a_k'x >= b_k+margin.
%   Each SEGMENT is assigned ONE face k per obstacle (from the previous iterate). Both endpoints are
%   constrained to that same half-space, which certifies the whole segment. Violations carry an exact
%   L1 slack penalty with growing weight (penalty CCP). Objective: sum of squared steps = shortest
%   path with uniform spacing (Cauchy-Schwarz), i.e. a constant-speed reference for a tracker.
%   Local solution only: several homotopy starts are solved and the shortest feasible one is kept.

assert(exist('cvx_begin', 'file') == 2, 'CVX not on path. Run cvx_setup first.');
if nargin < 4, par = struct(); end
par = fillDefaults(par, struct('N', 60, 'Ts', 0.1, 'vmax', 2.0, 'margin', 0.15, 'wAcc', 0, ...
    'tau0', 10, 'mu', 2, 'tauMax', 1e6, 'maxIt', 40, 'slackTol', 1e-6, 'objTol', 1e-5));
x0 = x0(:); xf = xf(:);
assert(norm(xf - x0) > 0, 'Start and goal coincide.');

H = cellfun(@polyToHalfspaces, obs(:)', 'UniformOutput', false);
inits = buildInits(x0, xf, obs, par);

t0 = tic;
best = [];
for c = 1:numel(inits)
    r = ccpSolve(inits{c}, x0, xf, H, par);
    if isempty(best) || betterThan(r, best), best = r; end
end
if ~best.feasible
    warning('dcp_polygon_planner:infeasible', ...
        'No start produced a collision-free path (%d segment/obstacle violations).', best.collisions);
end

N = par.N;
plan.x = best.x;
plan.t = (0:N) * par.Ts;
plan.v = [diff(best.x, 1, 2) / par.Ts, zeros(2,1)];
plan.L = best.L;
plan.feasible = best.feasible;
plan.collisions = best.collisions;
plan.iters = best.iters;
plan.time = toc(t0);
plan.par = par;
plan.refFun = @(tq) interp1(plan.t, plan.x', tq(:), 'linear', 'extrap')';  % shortcut: linear, upgrade to spline for C1 refs
end

%% ------------------------------------------------------------------------------------------
function par = fillDefaults(par, def)
f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(par, f{i}), par.(f{i}) = def.(f{i}); end
end
end

function tf = betterThan(a, b)
if a.feasible ~= b.feasible, tf = a.feasible; return; end
if a.feasible, tf = a.L < b.L; else, tf = a.collisions < b.collisions; end
end

function r = ccpSolve(z, x0, xf, H, par)
N = par.N; m = numel(H);
Lref = norm(xf - x0);
tau = par.tau0 / Lref;                  % keeps penalty above the objective gradient (~2/Lref)
Iprev = [];
Jprev = inf;
iters = 0;
for k = 1:par.maxIt
    iters = k;
    I = selectFaces(z, H, par.margin, Iprev);

    cvx_begin quiet
        variable x(2, N+1)
        variable s(N, m) nonnegative
        dx = x(:, 2:N+1) - x(:, 1:N);
        obj = (N / Lref^2) * sum_square(dx(:)) + tau * sum(s(:));
        if par.wAcc > 0
            ddx = dx(:, 2:N) - dx(:, 1:N-1);
            obj = obj + par.wAcc * (N^3 / Lref^2) * sum_square(ddx(:));
        end
        minimize(obj)
        subject to
            x(:,1) == x0;
            x(:,N+1) == xf;
            if isfinite(par.vmax)
                norms(dx, 2, 1) <= par.vmax * par.Ts;
            end
            for j = 1:m
                Nr = H{j}.A(I(:,j), :);                 % N x 2, face normal per segment
                dj = H{j}.b(I(:,j)) + par.margin;       % N x 1
                Nr(:,1).*x(1,1:N)'   + Nr(:,2).*x(2,1:N)'   - dj >= -s(:,j);
                Nr(:,1).*x(1,2:N+1)' + Nr(:,2).*x(2,2:N+1)' - dj >= -s(:,j);
            end
    cvx_end

    if ~any(strcmp(cvx_status, {'Solved', 'Inaccurate/Solved'}))
        warning('dcp_polygon_planner:solver', 'CVX status "%s" at iteration %d; stopping this start.', cvx_status, k);
        break;
    end

    z = x;
    maxS = max(s(:));
    J = (N / Lref^2) * sum(sum((z(:,2:end) - z(:,1:end-1)).^2));
    facesFixed = isequal(I, Iprev);
    relDec = abs(Jprev - J) <= par.objTol * abs(Jprev);
    if maxS <= par.slackTol && (facesFixed || relDec), break; end
    if maxS > par.slackTol, tau = min(par.mu * tau, par.tauMax); end
    Iprev = I; Jprev = J;
end

r.x = z;
r.L = sum(sqrt(sum(diff(z, 1, 2).^2, 1)));
r.collisions = countCollisions(z, H, par.margin);
r.feasible = (r.collisions == 0);
r.iters = iters;
end

function I = selectFaces(z, H, margin, Iprev)
% Face per (segment, obstacle): maximise the worse endpoint's signed clearance; keep the
% previous face while it is still satisfied (hysteresis against face flipping).
N = size(z, 2) - 1; m = numel(H);
I = zeros(N, m);
for j = 1:m
    G = H{j}.A * z - (H{j}.b + margin);          % K x (N+1)
    S = min(G(:,1:N), G(:,2:N+1));               % K x N
    [~, kmax] = max(S, [], 1);
    I(:,j) = kmax(:);
    if ~isempty(Iprev)
        keep = S(sub2ind(size(S), Iprev(:,j)', 1:N)) >= 0;
        I(keep, j) = Iprev(keep, j);
    end
end
end

function n = countCollisions(z, H, margin)
% Exact test (Cyrus-Beck clipping): does any segment enter the interior of an inflated polygon?
tol = 1e-6;
N = size(z, 2) - 1; n = 0;
for j = 1:numel(H)
    A = H{j}.A; b = H{j}.b + margin - tol;
    for i = 1:N
        a = z(:,i); d = z(:,i+1) - a;
        g0 = A*a - b; gs = A*d;                  % inside iff g0 + t*gs < 0 for all faces
        t0 = 0; t1 = 1; hit = true;
        for k = 1:numel(g0)
            if abs(gs(k)) < eps
                if g0(k) >= 0, hit = false; break; end
            else
                tb = -g0(k) / gs(k);
                if gs(k) > 0, t1 = min(t1, tb); else, t0 = max(t0, tb); end
            end
        end
        if hit && t0 < t1, n = n + 1; end
    end
end
end

function h = polyToHalfspaces(P)
% Convex polygon (vertices Nx2, any order) -> {x : A x <= b}, unit outward normals.
k = convhull(P(:,1), P(:,2));                    % counter-clockwise, closed
V = P(k(1:end-1), :);
if size(V,1) < size(P,1)
    warning('dcp_polygon_planner:nonconvex', 'Obstacle is not convex or has interior vertices; using its convex hull.');
end
e = circshift(V, -1, 1) - V;
A = [e(:,2), -e(:,1)];
A = A ./ sqrt(sum(A.^2, 2));
h.A = A;
h.b = sum(A .* V, 2);
end

function inits = buildInits(x0, xf, obs, par)
% Homotopy starts: straight line, each obstacle passed left/right, all-left, all-right.
N = par.N; m = numel(obs);
dir = (xf - x0) / norm(xf - x0);
perp = [-dir(2); dir(1)];
c = zeros(2, m); R = zeros(1, m);
for j = 1:m
    P = obs{j}';
    c(:,j) = mean(P, 2);
    R(j) = max(sqrt(sum((P - c(:,j)).^2, 1))) + par.margin;
end

sets = {struct('idx', [], 'sg', [])};
sets{1} = struct('idx', 1, 'sg', -1); %#ok<AGROW>
sets{2} = struct('idx', 1, 'sg', 1); %#ok<AGROW>
for j = 2:m
    for sg = [-1 1]
        sets{end+1} = struct('idx', j, 'sg', sg); %#ok<AGROW>
    end
end
if m > 1
    sets{end+1} = struct('idx', 1:m, 'sg',  ones(1,m));
    sets{end+1} = struct('idx', 1:m, 'sg', -ones(1,m));
end

inits = cell(1, numel(sets));
for q = 1:numel(sets)
    st = sets{q};
    via = c(:, st.idx) + perp * (st.sg .* R(st.idx));
    if ~isempty(via)
        [~, ord] = sort(dir' * via);
        via = via(:, ord);
    end
    pts = [x0, via, xf];
    arc = [0, cumsum(sqrt(sum(diff(pts, 1, 2).^2, 1)))];
    inits{q} = interp1(arc, pts', linspace(0, arc(end), N+1))';
end
end
