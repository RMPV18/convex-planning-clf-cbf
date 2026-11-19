function out = double_integrator_online_firi_clf_cbf(targets, obstacles, workspace, x0, z0, params, opts)
%DOUBLE_INTEGRATOR_ONLINE_FIRI_CLF_CBF  Secondary version of the safe-chain
%   tracker: the safe convex set of every mode is computed ONLINE with FIRI.
%
%   out = double_integrator_online_firi_clf_cbf(targets, obstacles, workspace, x0, z0)
%   out = double_integrator_online_firi_clf_cbf(..., params, opts)
%
%   targets    2 x N targets (offline: crossing points translated into the
%              seed quadrilaterals of build_safe_chain(..., 'seed'), so that
%              consecutive targets are mutually visible with clearance >= mu);
%              targets(:,N) = goal
%   obstacles  cell of convex obstacle pieces (convex_obstacle_pieces of the
%              inflated obstacles)
%   workspace  [xmin xmax; ymin ymax]
%   x0, z0     initial position and velocity
%   params     overrides of safe_chain_params
%   opts       delta (0.12): radius of the disc around x and the target
%                     included in the seed;
%              boxPad (2.0): padding of the heuristic bounding box around
%                     the seed (Gonzalez-Garcia et al.), clipped to the
%                     workspace;
%              dtCheck (0.05): period of the line-of-sight test;
%              firi: options of firi_polytope
%
%   MODE i = (S_i, targets(:,i)), S_i = FIRI(seed = conv(D(x) ∪ D(g_i)) ∩
%   workspace), D(p) = octagon circumscribing the disc of radius delta at p:
%   the largest FIRI polytope containing the current state and the target,
%   so the CLF-CBF backstepping controller (safe_chain_control) converges to
%   g_i inside S_i. Every dtCheck the after-next target g_{i+1} is tested for
%   line of sight from x (segment clearance >= circumradius of D); when in
%   sight, S_{i+1} = FIRI(conv(D(x) ∪ D(g_{i+1}))) is computed and the mode
%   jumps if h_{S_{i+1}}(x) >= sigma (beta_h reset Eq. 84). Since |x - g_i| -> 0
%   and the segment g_i -> g_{i+1} has clearance >= mu > circumradius, the
%   test succeeds in finite time; at most N-1 jumps.
%
%   out: as double_integrator_safe_chain_clf_cbf plus sets (cell of the
%   online polytopes), setsV (their vertices), firiTime (total online FIRI
%   time), firiCalls, wallTime (simulation wall-clock time, FIRI included).

if nargin < 6, params = struct(); end
if nargin < 7, opts = struct(); end
P = safe_chain_params(params);
def = struct('delta', 0.12, 'boxPad', 2.0, 'dtCheck', 0.05, 'firi', struct());
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(opts, fn{k}), opts.(fn{k}) = def.(fn{k}); end
end
t0 = tic;
x0 = x0(:);  z0 = z0(:);  n = numel(x0);
assert(n == 2, 'online_firi:dim', 'The online version is 2-D (FIRI seed and clearance tests).');
N = size(targets, 2);
rC = opts.delta / cos(pi / 8);                  % circumradius of the octagon
G = targets;  goal = G(:, N);

firiTime = 0;  firiCalls = 0;
[S, V, ft] = online_set(x0, G(:, 1), obstacles, workspace, opts);
firiTime = firiTime + ft;  firiCalls = firiCalls + 1;
sets = {S};  setsV = {V};
if ~segment_clear(x0, G(:, 1), obstacles, rC)
    error('online_firi:sight', 'The first target is not in line of sight of x0.');
end
if lse_cbf(x0, S.A, S.b, P.kappa) <= 0
    error('online_firi:seed', 'x0 is not inside the LSE safe set of the first online polytope.');
end
i = 1;  t = 0;  s = [x0; z0];
betaH = update_beta_h(s, G(:, i), S, P);
jumps = zeros(0, 4);
Tc = {};  Sc = {};  Ic = {};  Bc = {};  Xc = {};
converged = false;  nRejected = 0;
while t < P.t_max
    x = s(1:n);
    if i < N && segment_clear(x, G(:, i+1), obstacles, rC)
        [Sn, Vn, ft] = online_set(x, G(:, i+1), obstacles, workspace, opts);
        firiTime = firiTime + ft;  firiCalls = firiCalls + 1;
        if lse_cbf(x, Sn.A, Sn.b, P.kappa) >= P.sigma
            i = i + 1;  S = Sn;  sets{end+1} = Sn;  setsV{end+1} = Vn; %#ok<AGROW>
            betaH = update_beta_h(s, G(:, i), S, P);
            jumps(end+1, :) = [t, i-1, i, betaH]; %#ok<AGROW>
        else
            nRejected = nRejected + 1;
        end
    end
    if i == N && norm([x - goal; s(n+1:end)]) <= P.tol_final
        converged = true;  break
    end
    ii = numel(sets);  gg = G(:, i);  bb = betaH;
    ev = @(tt, ss) conv_event(ss, gg, goal, n, P);
    o = odeset('RelTol', P.rel_tol, 'AbsTol', P.abs_tol, 'MaxStep', P.max_step, 'Events', ev);
    [ts, ss, ~, ~, ie] = ode45(@(tt, ss) [ss(n+1:2*n); safe_chain_control(ss, gg, sets{ii}, bb, P)], ...
        [t, min(t + opts.dtCheck, P.t_max)], s', o);
    Tc{end+1} = ts';  Sc{end+1} = ss';  Ic{end+1} = ii * ones(1, numel(ts)); %#ok<AGROW>
    Bc{end+1} = betaH * ones(1, numel(ts));  Xc{end+1} = repmat(gg, 1, numel(ts)); %#ok<AGROW>
    t = ts(end);  s = ss(end, :)';
    if ~isempty(ie), converged = true; break; end
end
if ~converged
    warning('online_firi:timeout', 'Reached t_max = %g in mode %d without convergence.', P.t_max, i);
end

% ---- assemble and verify ------------------------------------------------
T = [Tc{:}];  Ssm = [Sc{:}];  I = [Ic{:}];  Bh = [Bc{:}];  X = [Xc{:}];
K = numel(T);
u = zeros(n, K);  h = zeros(1, K);  h1 = zeros(1, K);  mg = zeros(1, K);
maxBs = -inf;  nFb = 0;
for k = 1:K
    [u(:, k), aux] = safe_chain_control(Ssm(:, k), X(:, k), sets{I(k)}, Bh(k), P);
    h(k) = aux.h;  h1(k) = aux.h1;
    mg(k) = min(sets{I(k)}.b - sets{I(k)}.A * Ssm(1:n, k));
    nFb = nFb + aux.fallback;  maxBs = max(maxBs, aux.bsViol);
end
ver = struct('minH1', min(h1), 'minMargin', min(mg), 'maxBackstepViolation', maxBs, ...
    'nQpFallback', nFb, 'finalError', norm(Ssm(1:n, end) - goal), ...
    'nJumps', size(jumps, 1), 'converged', converged, 'nRejectedSets', nRejected);
if ver.minH1 < -1e-9
    warning('online_firi:safety', 'h1 became negative (min %.3e).', ver.minH1);
end
if ver.minMargin < -1e-9
    warning('online_firi:safety', 'The state left the active safe set (min margin %.3e).', ver.minMargin);
end
if nFb > 0
    warning('online_firi:qpFallback', '%d samples with incompatible CLF/CBF constraints.', nFb);
end
out = struct('t', T, 's', Ssm, 'x', Ssm(1:n, :), 'z', Ssm(n+1:end, :), 'u', u, 'mode', I, ...
    'target', X, 'h', h, 'h1', h1, 'betaH', Bh, 'margin', mg, 'jumps', jumps, ...
    'verification', ver, 'params', P, 'sets', {sets}, 'setsV', {setsV}, ...
    'firiTime', firiTime, 'firiCalls', firiCalls, 'wallTime', toc(t0));
end

% =========================================================================
function [S, V, ft] = online_set(x, g, obstacles, ws, opts)
% FIRI polytope seeded with conv(D(x) ∪ D(g)) ∩ workspace.
t0 = tic;
th = (0:7)' * pi / 4 + pi / 8;
oct = opts.delta / cos(pi / 8) * [cos(th) sin(th)];
H = [x' + oct; g' + oct];
H = H(convhull(H(:, 1), H(:, 2)), :);
box = polyshape([ws(1,1) ws(1,2) ws(1,2) ws(1,1)], [ws(2,1) ws(2,1) ws(2,2) ws(2,2)]);
seed = intersect(polyshape(H(:, 1), H(:, 2)), box);
seedV = seed.Vertices;
lo = max(min(seedV, [], 1) - opts.boxPad, ws(:, 1)');
hi = min(max(seedV, [], 1) + opts.boxPad, ws(:, 2)');
[A, b] = firi_polytope(seedV, obstacles, [eye(2); -eye(2)], [hi'; -lo'], opts.firi);
S = struct('A', A, 'b', b);
V = polytope_vertices_2d(A, b);
ft = toc(t0);
end

function ok = segment_clear(p, q, obstacles, r)
% Segment p-q at distance >= r from every convex obstacle piece (2-D).
ok = true;
for k = 1:numel(obstacles)
    if seg_poly_dist(p(:)', q(:)', obstacles{k}) < r, ok = false; return; end
end
end

function d = seg_poly_dist(p, q, U)
% Distance between segment p-q and the convex polygon conv(U) (2-D).
K = convhull(U(:, 1), U(:, 2));  U = U(K(1:end-1), :);
m = size(U, 1);
if inpolygon(p(1), p(2), U(:, 1), U(:, 2)), d = 0; return; end
d = inf;
for k = 1:m
    a = U(k, :);  b = U(mod(k, m) + 1, :);
    if seg_cross(p, q, a, b), d = 0; return; end
    d = min([d, pt_seg(p, a, b), pt_seg(q, a, b), pt_seg(a, p, q)]);
end
end

function tf = seg_cross(p1, p2, q1, q2)
c = @(o, a, b) (a(1) - o(1)) * (b(2) - o(2)) - (a(2) - o(2)) * (b(1) - o(1));
d1 = c(q1, q2, p1);  d2 = c(q1, q2, p2);  d3 = c(p1, p2, q1);  d4 = c(p1, p2, q2);
tf = d1 * d2 <= 0 && d3 * d4 <= 0;
end

function d = pt_seg(x, a, b)
ab = b - a;  t = max(0, min(1, ((x - a) * ab') / max(ab * ab', eps)));
d = norm(x - (a + t * ab));
end

function betaH = update_beta_h(s, xh, S, P)
n = numel(xh);
[~, aux] = safe_chain_control(s, xh, S, P.beta_h_nom, P);
if aux.h <= 0
    error('online_firi:safety', 'Cannot set beta_h: h = %.3e <= 0 at the mode update.', aux.h);
end
e = s(n+1:end) - aux.k;
betaH = max(P.beta_h_nom, 1.05 * (e' * e) / (2 * aux.h));
end

function [val, term, dir] = conv_event(s, g, goal, n, P)
val = 1;  term = 1;  dir = -1;
if isequal(g, goal), val = norm([s(1:n) - goal; s(n+1:2*n)]) - P.tol_final; end
end
