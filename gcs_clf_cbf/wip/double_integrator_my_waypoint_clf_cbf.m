function out = double_integrator_my_waypoint_clf_cbf(waypoints, regions, x0, z0, params)
%DOUBLE_INTEGRATOR_WAYPOINT_CLF_CBF  Hybrid CLF-CBF QP waypoint follower with
%   joint CLF-CBF backstepping, double integrator.
%
%   out = double_integrator_waypoint_clf_cbf(waypoints, regions, x0, z0)
%   out = double_integrator_waypoint_clf_cbf(waypoints, regions, x0, z0, params)
%
% PLANT  xdot = z, zdot = u,  x, z, u in R^n  (Matias & Silvestre Problem 2, r = 1)
%
% INPUTS
%   waypoints  n x M   waypoints w_1..w_M (last = desired point)
%   regions    1 x M   struct array A, b: safe set while targeting w_k
%                      convex: A (m x n), b (m x 1); union of 2: A = {A1,A2}, b = {b1,b2}
%   x0, z0     n x 1   initial position and velocity
%   params     struct  overrides of the fields of WP_PARAMS (end of file)
%
% TOP LEVEL (first-order system xdot = v), closed-form CLF-CBF QP:
%   V = 1/2 |x - w_k|^2,   h = h_region,k(x)                (WP_REGION_CBF)
%   k(x) = argmin 1/2|v|^2  s.t.  gV'v <= -2 gamma_bar V
%                                 gh'v >= -alpha_bar h
%   (CLF_CBF_QP_CLOSED_FORM). Compatibility: Eq. 59 (alpha_bar >= gamma_bar) holds
%   for convex regions; for a union it is monitored, the CBF is prioritised and
%   warning('wp:qpFallback') is raised.
%
% BACKSTEPPING (Eq. 72-77)   e = z - k(x),
%   V1 = V + |e|^2/(2 beta_v),   h1 = h - |e|^2/(2 beta_h)
%   kdot = (dk/dx) z   (central differences; k is Lipschitz, piecewise smooth,
%                       not C^1 at regime boundaries)
%   Both constraints reduce to  e'u <= B_V  and  e'u <= B_h:
%     B_V = -beta_v [ gV'z + 2 gamma_bar gamma1_ratio V1 ] + e'kdot
%     B_h =  beta_h [ gh'z + alpha_bar alpha1_ratio h1 ] + e'kdot
%   Closed-form QP min 1/2|u|^2 s.t. e'u <= min(B_V,B_h) = B:
%     u = 0 if B >= 0, else u = (B/|e|^2) e.
%
% SWITCHING k -> k+1 when |x - w_k| <= r_k (WP_PREPARE_REGIONS). Target and safe
%   set change together; beta_h+ = max(beta_h_nom, 1.05 |e+|^2/(2 h_k+1(x)))
%   (Eq. 84) keeps h1 >= 0. Sight (h_k+1(x) >= sigma and segment x -> w_{k+1}
%   inside the new region) is verified at every jump.
%
% OUTPUT struct: t, s, y, p, u, mode, h, h1, betaH, bodyMargin, jumps
%   [t from to beta_h], switchRadius, verification, regions, waypoints, params.
%   verification: minH1, minBodyMargin, (top-level CBF is satisfied by the QP by construction),
%   maxBackstepViolation, nQpFallback, finalOutputError, nJumps, converged.
%
%
% SELF-CONTAINED: all helpers are local functions at the end of this file
%   (parameters, region preparation, region CBF, geometry, closed-form QP).
%
% LIMITATIONS  Guarantees apply to the interior of each safe set. Union regions
%   are non-convex: deadlock-freedom is not guaranteed; a stall is reported by
%   warning('wp:timeout').

if nargin < 5
    params = struct();
end
P = wp_params(params);
x0 = x0(:);
z0 = z0(:);
n = size(waypoints, 1);
if numel(x0) ~= n || numel(z0) ~= n
    error('wp:input', 'x0 and z0 must have %d elements.', n);
end
[regions, r] = wp_prepare_regions(waypoints, regions, 0, x0, P);
M = size(waypoints, 2);

% ======================= hybrid simulation ================================
k = 1;
t = 0;
s = [x0; z0];
betaH = update_beta_h(s, waypoints(:, k), regions(k), n, P);
jumps = zeros(0, 4);
Tc = {};
Sc = {};
Kc = {};
Bc = {};
converged = false;

while true
    while k < M && norm(s(1:n) - waypoints(:, k)) <= r(k)
        [k, betaH, jumps] = perform_jump(s, t, k, betaH, waypoints, regions, jumps, n, P);
    end
    if k == M && norm([s(1:n) - waypoints(:, M); s(n + 1:2 * n)]) <= P.tol_final
        Tc{end + 1} = t; Sc{end + 1} = s; Kc{end + 1} = k; Bc{end + 1} = betaH; %#ok<AGROW>
        converged = true;
        break
    end

    w = waypoints(:, k);
    reg = regions(k);
    kk = k;
    bh = betaH;
    opts = odeset('RelTol', P.rel_tol, 'AbsTol', P.abs_tol, 'MaxStep', P.max_step, ...
        'Events', @(tt, ss) event_fcn(ss, w, kk, M, r, n, P));
    [ts, ss, te, se, ie] = ode45(@(tt, ss) closed_loop(ss, w, reg, bh, n, P), [t, P.t_max], s.', opts);

    Tc{end + 1} = ts(:).'; %#ok<AGROW>
    Sc{end + 1} = ss.'; %#ok<AGROW>
    Kc{end + 1} = k * ones(1, numel(ts)); %#ok<AGROW>
    Bc{end + 1} = betaH * ones(1, numel(ts)); %#ok<AGROW>

    if isempty(ie)
        warning('wp:timeout', 'Reached t_max = %g in mode %d without a switching/convergence event.', P.t_max, k);
        break
    end
    if k == M
        converged = true;
        break
    end
    t = te(end);
    s = se(end, :).';
    [k, betaH, jumps] = perform_jump(s, t, k, betaH, waypoints, regions, jumps, n, P);
end

out = postprocess_and_verify(Tc, Sc, Kc, Bc, jumps, converged, waypoints, regions, r, n, P);
end

% =========================================================================
function ds = closed_loop(s, w, reg, betaH, n, P)
u = control_law(s, w, reg, betaH, n, P);
ds = [s(n + 1:2 * n); u];
end

% =========================================================================
function [u, aux] = control_law(s, w, reg, betaH, n, P)
%CONTROL_LAW  Backstepping CLF-CBF QP in closed form.
x = s(1:n);
z = s(n + 1:2 * n);
[kx, tq] = top_level_qp(x, w, reg, P);
J = zeros(n, n);                                   % dk/dx, central differences
for j = 1:n
    d = P.fd_step * (1 + abs(x(j)));
    xp = x;
    xm = x;
    xp(j) = xp(j) + d;
    xm(j) = xm(j) - d;
    J(:, j) = (top_level_qp(xp, w, reg, P) - top_level_qp(xm, w, reg, P)) / (2 * d);
end
kdot = J * z;
e = z - kx;
ee = e.' * e;
V1 = tq.V + ee / (2 * P.beta_v);
h1 = tq.h - ee / (2 * betaH);
BV = -P.beta_v * (tq.gV.' * z + 2 * P.gamma_bar * P.gamma1_ratio * V1) + e.' * kdot;
Bh = betaH * (tq.gh.' * z + P.alpha_bar * P.alpha1_ratio * h1) + e.' * kdot;
B = min(BV, Bh);
if B >= 0 || ee < P.e_tol^2
    u = zeros(n, 1);
else
    u = (B / ee) * e;
end
aux = struct('h', tq.h, 'h1', h1, 'fallback', ~tq.feasible, 'bsViol', e.' * u - B);
end

% =========================================================================
function [v, aux] = top_level_qp(x, w, reg, P)
%TOP_LEVEL_QP  Closed-form CLF-CBF QP for xdot = v.
[h, gh] = wp_region_cbf(x, reg, P);
gV = x - w;
V = 0.5 * (gV.' * gV);
[v, ~, feasible] = clf_cbf_qp_closed_form(0, gV.', 2 * P.gamma_bar * V, 0, gh.', P.alpha_bar * h, P.tol_qp);
aux = struct('h', h, 'V', V, 'gV', gV, 'gh', gh, 'feasible', feasible);
end

% =========================================================================
function [val, term, dir] = event_fcn(s, w, k, M, r, n, P)
%EVENT_FCN  Jump event (k < M) or convergence event (k = M).
if k < M
    val = norm(s(1:n) - w) - r(k);
else
    val = norm([s(1:n) - w; s(n + 1:2 * n)]) - P.tol_final;
end
term = 1;
dir = -1;
end

% =========================================================================
function betaH = update_beta_h(s, w, reg, n, P)
%UPDATE_BETA_H  Eq. (84): guarantees h1 >= 0 at the mode update.
[kx, tq] = top_level_qp(s(1:n), w, reg, P);
if tq.h <= 0
    error('wp:safety', 'Cannot set beta_h: h = %.3e <= 0 at the mode update.', tq.h);
end
e = s(n + 1:2 * n) - kx;
betaH = max(P.beta_h_nom, 1.05 * (e.' * e) / (2 * tq.h));
end

% =========================================================================
function [k, betaH, jumps] = perform_jump(s, t, k, betaH, waypoints, regions, jumps, n, P)
%PERFORM_JUMP  Switch CLF target and CBF safe set from mode k to k+1.
x = s(1:n);
hn = wp_region_cbf(x, regions(k + 1), P);
if hn < P.sigma - 1e-9
    error('wp:sight', 'Sight condition failed at jump %d -> %d (h_next = %.3e < sigma = %.3e).', k, k + 1, hn, P.sigma);
end
if ~wp_segment_in_region(x, waypoints(:, k + 1), regions(k + 1), -1e-9)
    error('wp:sight', 'Segment from x to w_%d leaves the safe set of region %d at the jump.', k + 1, k + 1);
end
from = k;
k = k + 1;
betaH = update_beta_h(s, waypoints(:, k), regions(k), n, P);
jumps(end + 1, :) = [t, from, k, betaH]; %#ok<AGROW>
end

% =========================================================================
function out = postprocess_and_verify(Tc, Sc, Kc, Bc, jumps, converged, waypoints, regions, r, n, P)
%POSTPROCESS_AND_VERIFY  Assemble the output and run the a-posteriori checks.
M = size(waypoints, 2);
T = [Tc{:}];
S = [Sc{:}];
modeIdx = [Kc{:}];
betaHs = [Bc{:}];
N = numel(T);
u = zeros(n, N);
h = zeros(1, N);
h1 = zeros(1, N);
bodyMargin = zeros(1, N);
maxBs = -inf;
nFallback = 0;
for i = 1:N
    kk = modeIdx(i);
    [u(:, i), aux] = control_law(S(:, i), waypoints(:, kk), regions(kk), betaHs(i), n, P);
    h(i) = aux.h;
    h1(i) = aux.h1;
    bodyMargin(i) = wp_region_margin(S(1:n, i), regions(kk));
    nFallback = nFallback + aux.fallback;
    maxBs = max(maxBs, aux.bsViol);
end
ver = struct('minH1', min(h1), 'minBodyMargin', min(bodyMargin), ...
    'maxBackstepViolation', maxBs, 'nQpFallback', nFallback, ...
    'finalOutputError', norm(S(1:n, end) - waypoints(:, M)), ...
    'nJumps', size(jumps, 1), 'converged', converged);

tolSafe = 1e-9;
if ver.minH1 < -tolSafe
    warning('wp:safety', 'h1 became negative (min = %.3e).', ver.minH1);
end
if ver.minBodyMargin < -tolSafe
    warning('wp:safety', 'Position left the active region (min margin = %.3e).', ver.minBodyMargin);
end
if ver.maxBackstepViolation > 1e-7
    warning('wp:constraint', 'Backstepping constraint violated by %.3e.', ver.maxBackstepViolation);
end
if nFallback > 0
    warning('wp:qpFallback', '%d samples had incompatible CLF/CBF constraints; the CBF was prioritised.', nFallback);
end
if converged && ver.nJumps ~= M - 1
    warning('wp:jumps', 'Number of jumps %d differs from M-1 = %d.', ver.nJumps, M - 1);
end
out = struct('t', T, 's', S, 'y', S(1:n, :), 'p', S(1:n, :), 'u', u, 'mode', modeIdx, ...
    'h', h, 'h1', h1, 'betaH', betaHs, 'bodyMargin', bodyMargin, 'jumps', jumps, ...
    'switchRadius', r, 'verification', ver, 'regions', regions, ...
    'waypoints', waypoints, 'params', P);
end

% =========================================================================
% SHARED HELPERS (local copies)
% =========================================================================

% =========================================================================
% PARAMETERS
% =========================================================================
function P = wp_params(user)
%WP_PARAMS  Defaults and validation of the controller parameters (shared).
%   P = wp_params(struct('kappa', 30))  overrides single fields.
%   beta_v, beta_h_nom, gamma1_ratio, alpha1_ratio, fd_step, e_tol are used by
%   the second-order controller only.
if nargin < 1
    user = struct();
end
P = struct( ...
    'kappa', 20, ...         % LSE smoothing of the polytope CBF
    'kappa_union', 20, ...   % smooth-max sharpness for a union of 2 polytopes
    'gamma_bar', 1.0, ...    % CLF rate: gamma(s) = 2 gamma_bar s
    'alpha_bar', 2.0, ...    % CBF rate: alpha(s) = alpha_bar s, alpha_bar >= gamma_bar
    'beta_v', 1.0, ...       % backstepping CLF weight
    'beta_h_nom', 1.0, ...   % nominal backstepping CBF weight
    'gamma1_ratio', 0.9, ... % gamma1 = ratio * gamma, ratio in (0,1]
    'alpha1_ratio', 1.1, ... % alpha1 = ratio * alpha, ratio >= 1
    'sigma', 0.05, ...       % sight margin h_eff,k+1(y) >= sigma at jumps
    'r_sw', 1.0, ...         % maximum switching radius
    'fd_step', 1e-6, ...     % central-difference step for dk/dx
    'tol_qp', 1e-12, ...     % gradient-norm tolerance in the closed-form QP
    'e_tol', 1e-12, ...      % |e| below which the backstepping input is 0
    'tol_final', 1e-4, ...   % termination residual in the last mode
    't_max', 200, ...        % absolute time limit
    'max_step', 0.05, ...    % ode45 maximum step
    'rel_tol', 1e-9, ...
    'abs_tol', 1e-11);
fn = fieldnames(user);
for i = 1:numel(fn)
    if ~isfield(P, fn{i})
        error('wp:params', 'Unknown parameter "%s".', fn{i});
    end
    P.(fn{i}) = user.(fn{i});
end
if ~(P.gamma_bar > 0 && P.alpha_bar >= P.gamma_bar)
    error('wp:params', 'Require 0 < gamma_bar <= alpha_bar (Matias & Silvestre Eq. 59).');
end
if ~(P.gamma1_ratio > 0 && P.gamma1_ratio <= 1 && P.alpha1_ratio >= 1)
    error('wp:params', 'Require gamma1_ratio in (0,1] and alpha1_ratio >= 1.');
end
if ~(P.sigma > 0 && P.r_sw > 0 && P.kappa > 0 && P.kappa_union > 0 && P.beta_v > 0 ...
        && P.beta_h_nom > 0 && P.fd_step > 0)
    error('wp:params', 'sigma, r_sw, kappa, kappa_union, beta_v, beta_h_nom, fd_step must be positive.');
end
end

% =========================================================================
% REGION PREPARATION
% =========================================================================
function [regions, r] = wp_prepare_regions(waypoints, regions, rho, y0, P)
%WP_PREPARE_REGIONS  Validate regions, normalise them and compute switching radii.
%
%   [regions, r] = wp_prepare_regions(waypoints, regions, rho, y0, P)
%
% INPUT regions(k).A, regions(k).b : safe set while targeting waypoints(:,k)
%     convex polytope  : A (m x n), b (m x 1)           {A y <= b}
%     union of 2       : A = {A1, A2}, b = {b1, b2}
%   rho  safety radius (CBF h_eff = h_region - rho);  y0  initial controlled output.
% OUTPUT
%   regions(k).sets  cell of 1 or 2 polytopes with unit-norm rows
%   r(k)             switching radius of the jump k -> k+1, k = 1..M-1
%
% CHECKS (error on violation), h_eff = wp_region_cbf - rho:
%   h_eff,k(w_k) > 0,   h_eff,k+1(w_k) > sigma,   h_eff,1(y0) > 0, and for a union
%   region k+1 the segment w_k -> w_{k+1} lies in the union.
%
% SWITCHING RADIUS  r_k <= min(r_sw, h_eff,k+1(w_k) - sigma)  (1-Lipschitz h gives
%   h_eff,k+1(y) >= sigma for |y - w_k| <= r_k). Convex region k+1: convexity
%   gives sight y -> w_{k+1}. Union region k+1: r_k is the largest radius (bisection,
%   factor 0.99) such that segment w_k -> w_{k+1} lies in the union of the
%   polytopes shrunk by rho + r_k; every segment y -> w_{k+1} with |y - w_k| <= r_k
%   stays within r_k of it. Conservative.
n = size(waypoints, 1);
M = size(waypoints, 2);
if numel(regions) ~= M
    error('wp:input', 'regions must have one entry per waypoint (%d expected, %d given).', M, numel(regions));
end
nreg = struct('sets', cell(1, M));
for k = 1:M
    nreg(k).sets = normalise_region(regions(k), n, k);
end
regions = nreg;

r = zeros(1, max(M - 1, 0));
for k = 1:M
    hk = wp_region_cbf(waypoints(:, k), regions(k), P) - rho;
    if hk <= 0
        error('wp:input', ['Waypoint %d is not in the interior of its safe set (h_eff = %.3g <= 0). ' ...
            'Increase kappa, reduce rho or move it inward.'], k, hk);
    end
    if k < M
        % hn = wp_region_cbf(waypoints(:, k), regions(k + 1), P) - rho;
        % if hn <= P.sigma
        %     error('wp:input', ['Waypoint %d is not deep enough in region %d (h_eff = %.3g <= sigma = %.3g). ' ...
        %         'Regions must overlap around w_%d.'], k, k + 1, hn, P.sigma, k);
        % end
        r(k) = switch_radius(waypoints(:, k), waypoints(:, k + 1), regions(k + 1), hn, rho, P, k);
    end
end
h0 = wp_region_cbf(y0, regions(1), P) - rho;
if h0 <= 0
    error('wp:input', 'Initial output is not in the interior of the first safe set (h_eff = %.3g).', h0);
end
end

function sets = normalise_region(reg, n, k)
A = reg.A;
b = reg.b;
if ~iscell(A)
    A = {A};
end
if ~iscell(b)
    b = {b};
end
if numel(A) ~= numel(b) || numel(A) < 1 || numel(A) > 2
    error('wp:input', 'regions(%d): A and b must be a matrix/vector or cells of 1 or 2 polytopes.', k);
end
sets = cell(1, numel(A));
for i = 1:numel(A)
    Ai = A{i};
    bi = b{i}(:);
    if size(Ai, 2) ~= n || size(Ai, 1) ~= numel(bi)
        error('wp:input', 'regions(%d), set %d: A must be m x %d and b m x 1.', k, i, n);
    end
    nr = sqrt(sum(Ai.^2, 2));
    if any(nr <= 0)
        error('wp:input', 'regions(%d), set %d has a zero row in A.', k, i);
    end
    sets{i} = struct('A', Ai ./ nr, 'b', bi ./ nr);   % unit normals -> 1-Lipschitz h
end
end

function rk = switch_radius(wk, wn, regNext, hn, rho, P, k)
rmax = min(P.r_sw, hn - P.sigma);
cond = @(q) wp_segment_in_region(wk, wn, regNext, rho + q);
if ~cond(0)
    error('wp:input', 'Segment w_%d -> w_%d is not inside the union region %d: no sight between consecutive waypoints.', k, k + 1, k + 1);
end
if cond(rmax)
    rk = 0.99 * rmax;
    return
end
lo = 0;
hi = rmax;
for it = 1:60
    mid = 0.5 * (lo + hi);
    if cond(mid)
        lo = mid;
    else
        hi = mid;
    end
end
rk = 0.99 * lo;
if rk <= 0
    error('wp:input', 'Switching radius for jump %d -> %d is zero: segment w_%d -> w_%d touches the boundary of union region %d.', k, k + 1, k, k + 1, k + 1);
end
end

% =========================================================================
% REGION CBF AND GEOMETRY
% =========================================================================
function [h, gh] = wp_region_cbf(y, reg, P)
%WP_REGION_CBF  CBF of a normalised region (see WP_PREPARE_REGIONS).
%   One polytope : h = -(1/kappa) ln sum_i exp(kappa (a_i'y - b_i))
%                  concave, 1-Lipschitz (unit rows), h >= 0 => y in the polytope.
%   Two polytopes: h = smoothmax(h_1,h_2) - ln(2)/kappa_union <= max(h_1,h_2)
%                  (inner approximation of the union, 1-Lipschitz, not concave).
%   gh = gradient of h with respect to y.
ns = numel(reg.sets);
hs = zeros(1, ns);
gs = zeros(numel(y), ns);
for i = 1:ns
    [hs(i), gs(:, i)] = lse_cbf(y, reg.sets{i}.A, reg.sets{i}.b, P.kappa);
end
if ns == 1
    h = hs;
    gh = gs;
    return
end
m = max(hs);
wt = exp(P.kappa_union * (hs - m));
S = sum(wt);
h = m + log(S) / P.kappa_union - log(2) / P.kappa_union;
gh = gs * (wt.' / S);
end

function [h, gh] = lse_cbf(x, A, b, kappa)
g = A * x - b;
gmax = max(g);
ex = exp(kappa * (g - gmax));
S = sum(ex);
h = -(gmax + log(S) / kappa);
gh = -A.' * (ex / S);
end

function ok = wp_segment_in_region(p0, p1, reg, erode)
%WP_SEGMENT_IN_REGION  Exact test: segment p0->p1 inside the union of the region's
%   polytopes, each shrunk by "erode" (unit-norm rows). Each polytope clips the
%   segment to an interval of t in [0,1]; inside iff the intervals cover [0,1].
d = p1 - p0;
iv = zeros(0, 2);
for i = 1:numel(reg.sets)
    A = reg.sets{i}.A;
    c = reg.sets{i}.b - erode - A * p0;
    ad = A * d;
    lo = 0;
    hi = 1;
    feas = true;
    for j = 1:numel(c)
        if abs(ad(j)) < 1e-14
            if c(j) < -1e-12
                feas = false;
            end
        elseif ad(j) > 0
            hi = min(hi, c(j) / ad(j));
        else
            lo = max(lo, c(j) / ad(j));
        end
    end
    if feas && lo <= hi + 1e-12
        iv(end + 1, :) = [lo, hi]; %#ok<AGROW>
    end
end
ok = false;
if isempty(iv)
    return
end
iv = sortrows(iv);
if iv(1, 1) > 1e-12
    return
end
reach = iv(1, 2);
for i = 2:size(iv, 1)
    if iv(i, 1) <= reach + 1e-12
        reach = max(reach, iv(i, 2));
    end
end
ok = reach >= 1 - 1e-12;
end

function mg = wp_region_margin(p, reg)
%WP_REGION_MARGIN  max over the region's polytopes of min_i (b_i - a_i'p).
%   mg >= 0 iff p lies in the region.
mg = -inf;
for i = 1:numel(reg.sets)
    mg = max(mg, min(reg.sets{i}.b - reg.sets{i}.A * p));
end
end

% =========================================================================
% CLOSED-FORM CLF-CBF QP
% =========================================================================
function [u, regime, feasible] = clf_cbf_qp_closed_form(LfV, LGV, gammaV, Lfh, LGh, alphah, tol)
%CLF_CBF_QP_CLOSED_FORM  Exact solution of the hard-constrained CLF-CBF QP
%
%       u = argmin 1/2 |u|^2
%           s.t.  LfV + LGV u <= -gammaV        (CLF)
%                 Lfh + LGh u >= -alphah        (CBF)
%
%   for a single CLF/CBF pair (Matias & Silvestre, Eq. 60 = Eq. 9 with p = inf;
%   closed form Eq. 10-15 with p^{-1} = 0).
%
%   [u, regime, feasible] = clf_cbf_qp_closed_form(LfV, LGV, gammaV, Lfh, LGh, alphah)
%
% INPUTS
%   LfV, Lfh   scalars  Lie derivatives along the drift f
%   LGV, LGh   1 x q    Lie derivatives along G (row vectors)
%   gammaV     scalar   gamma(V(x))
%   alphah     scalar   alpha(h(x))
%   tol        scalar   optional numerical tolerance (default 1e-12)
%
% OUTPUTS
%   u         q x 1   minimiser
%   regime    1..4    active set (Eq. 14): 1 CLF only (S1), 2 CBF only (S2),
%                     3 both (S3), 4 none (S4, u = 0)
%   feasible  logical false if the two constraints are incompatible
%                     (aligned gradients, Eq. 48-53) or the CLF gradient
%                     vanishes while the CLF constraint is violated. In that
%                     case the CBF-only solution (safety priority) is returned.
%
% With FV = LfV + gammaV, Fh = Lfh + alphah, L = LGV LGh',
% Delta = L^2 - |LGV|^2 |LGh|^2:
%   S4: FV <= 0 and Fh >= 0                         u = 0
%   S1: u = -FV LGV' / |LGV|^2        (CBF satisfied)
%   S2: u = -Fh LGh' / |LGh|^2        (CLF satisfied)
%   S3: u = -l1 LGV' + l2 LGh',
%       l1 = (L Fh - |LGh|^2 FV)/Delta, l2 = (|LGV|^2 Fh - L FV)/Delta, l1,l2 >= 0.
% The candidates are tested in the order S4, S1, S2, S3, which yields the
% unique KKT point of the strictly convex QP.

if nargin < 7 || isempty(tol)
    tol = 1e-12;
end
a = LGV(:);
g = LGh(:);
q = numel(a);
if numel(g) ~= q
    error('clf_cbf_qp:dim', 'LGV and LGh must have the same length.');
end
FV = LfV + gammaV;
Fh = Lfh + alphah;
na = a.' * a;
ng = g.' * g;
feasible = true;

% S4: both constraints inactive at u = 0
if FV <= 0 && Fh >= 0
    u = zeros(q, 1);
    regime = 4;
    return
end

u1 = [];
u2 = [];
% S1: CLF only
if na > tol^2 && FV > 0
    u1 = -(FV / na) * a;
    if g.' * u1 >= -Fh - 1e-12
        u = u1;
        regime = 1;
        return
    end
end
% S2: CBF only
if ng > tol^2 && Fh < 0
    u2 = -(Fh / ng) * g;
    if a.' * u2 <= -FV + 1e-12
        u = u2;
        regime = 2;
        return
    end
end
% S3: both active
L = a.' * g;
Delta = L^2 - na * ng;
if abs(Delta) > 1e-12 * max(1, na * ng)
    l1 = (L * Fh - ng * FV) / Delta;
    l2 = (na * Fh - L * FV) / Delta;
    if l1 >= -1e-9 && l2 >= -1e-9
        u = -l1 * a + l2 * g;
        regime = 3;
        return
    end
end

% incompatible constraints: keep safety, report infeasibility
feasible = false;
regime = 2;
if isempty(u2)
    if ng > tol^2 && Fh < 0
        u2 = -(Fh / ng) * g;
    else
        u2 = zeros(q, 1);
    end
end
u = u2;
end
