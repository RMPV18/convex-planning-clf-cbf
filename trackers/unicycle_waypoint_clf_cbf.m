function out = unicycle_waypoint_clf_cbf(waypoints, regions, s0, l, params)
%UNICYCLE_WAYPOINT_CLF_CBF  Hybrid CLF-CBF QP waypoint follower, kinematic unicycle.
%
%   out = unicycle_waypoint_clf_cbf(waypoints, regions, s0, l)
%   out = unicycle_waypoint_clf_cbf(waypoints, regions, s0, l, params)
%
% PLANT  pdot = v (cos th, sin th), thdot = w, u = (v, w), s = (px, py, th).
%   Mixed relative degree -> first-order fully actuated system by controlling the
%   offset point y = p + l (cos th, sin th):   ydot = R(th) diag(1, l) u = G(s) u
%   (Matias & Silvestre Sec. V-D). The body is protected by rho = l.
%   The body p converges to a ball of radius l around the last waypoint.
%
% INPUTS
%   waypoints  2 x M   waypoints w_1..w_M (last = desired point of y)
%   regions    1 x M   struct array A, b: safe set while targeting w_k
%                      convex: A (m x 2), b (m x 1); union of 2: A = {A1,A2}, b = {b1,b2}
%   s0         3 x 1   initial state (px, py, th)
%   l          offset distance (> 0)
%   params     struct  overrides of the fields of WP_PARAMS (end of file)
%
% MODE k = (CLF target w_k, safe set Q_k)
%   V = 1/2 |y - w_k|^2,   h = h_region,k(y) - l        (WP_REGION_CBF)
%   u = argmin 1/2|u|^2  s.t.  LGV u <= -2 gamma_bar V
%                              LGh u >= -alpha_bar h
%   (closed form, CLF_CBF_QP_CLOSED_FORM; LGV = (y-w)'G, LGh = grad h' G)
%   Incompatible constraints are detected, the CBF is prioritised and
%   warning('wp:qpFallback') is raised.
%
% SWITCHING k -> k+1 when |y - w_k| <= r_k (WP_PREPARE_REGIONS). Target and safe
%   set change together. Sight is verified at every jump: h_eff,k+1(y) >= sigma and
%   segment y -> w_{k+1} inside the new region shrunk by l.
%
% OUTPUT struct: t, s, y, p (body), u, mode, h, bodyMargin, jumps [t from to],
%   switchRadius, verification, regions, waypoints, params.
%   verification: minH1 (= min h_eff), minBodyMargin, maxCbfViolation,
%   maxClfViolation, nQpFallback, finalOutputError, nJumps, converged.
%
%
% SELF-CONTAINED: all helpers are local functions at the end of this file
%   (parameters, region preparation, region CBF, geometry, closed-form QP).
%
% LIMITATION  Union regions are non-convex: deadlock-freedom is not guaranteed;
%   a stall is reported by warning('wp:timeout').

if nargin < 5
    params = struct();
end
P = wp_params(params);
s0 = s0(:);
if numel(s0) ~= 3 || l <= 0
    error('wp:input', 's0 must be 3 x 1 and l > 0.');
end
[regions, r] = wp_prepare_regions(waypoints, regions, l, offset_point(s0, l), P);
M = size(waypoints, 2);

% ======================= hybrid simulation ================================
k = 1;
t = 0;
s = s0;
jumps = zeros(0, 3);
Tc = {};
Sc = {};
Kc = {};
converged = false;

while true
    while k < M && norm(offset_point(s, l) - waypoints(:, k)) <= r(k)
        [k, jumps] = perform_jump(s, t, k, waypoints, regions, jumps, l, P);
    end
    if k == M && norm(offset_point(s, l) - waypoints(:, M)) <= P.tol_final
        Tc{end + 1} = t; Sc{end + 1} = s; Kc{end + 1} = k; %#ok<AGROW>
        converged = true;
        break
    end

    w = waypoints(:, k);
    reg = regions(k);
    kk = k;
    opts = odeset('RelTol', P.rel_tol, 'AbsTol', P.abs_tol, 'MaxStep', P.max_step, ...
        'Events', @(tt, ss) event_fcn(ss, w, kk, M, r, l, P));
    [ts, ss, te, se, ie] = ode45(@(tt, ss) closed_loop(ss, w, reg, l, P), [t, P.t_max], s.', opts);

    Tc{end + 1} = ts(:).'; %#ok<AGROW>
    Sc{end + 1} = ss.'; %#ok<AGROW>
    Kc{end + 1} = k * ones(1, numel(ts)); %#ok<AGROW>

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
    [k, jumps] = perform_jump(s, t, k, waypoints, regions, jumps, l, P);
end

out = postprocess_and_verify(Tc, Sc, Kc, jumps, converged, waypoints, regions, r, l, P);
end

% =========================================================================
function ds = closed_loop(s, w, reg, l, P)
u = control_law(s, w, reg, l, P);
ds = [u(1) * cos(s(3)); u(1) * sin(s(3)); u(2)];
end

% =========================================================================
function [u, aux] = control_law(s, w, reg, l, P)
%CONTROL_LAW  Closed-form CLF-CBF QP for ydot = G(s) u.
th = s(3);
G = [cos(th) -sin(th); sin(th) cos(th)] * diag([1, l]);
y = offset_point(s, l);
[hp, gh] = wp_region_cbf(y, reg, P);
h = hp - l;
gV = y - w;
V = 0.5 * (gV.' * gV);
LGV = gV.' * G;
LGh = gh.' * G;
gammaV = 2 * P.gamma_bar * V;
alphah = P.alpha_bar * h;
[u, ~, feasible] = clf_cbf_qp_closed_form(0, LGV, gammaV, 0, LGh, alphah, P.tol_qp);
aux = struct('h', h, 'fallback', ~feasible, ...
    'clfViol', LGV * u + gammaV, 'cbfViol', -(LGh * u + alphah));
end

% =========================================================================
function y = offset_point(s, l)
y = s(1:2) + l * [cos(s(3)); sin(s(3))];
end

% =========================================================================
function [val, term, dir] = event_fcn(s, w, k, M, r, l, P)
%EVENT_FCN  Jump event (k < M) or convergence event (k = M).
if k < M
    val = norm(offset_point(s, l) - w) - r(k);
else
    val = norm(offset_point(s, l) - w) - P.tol_final;
end
term = 1;
dir = -1;
end

% =========================================================================
function [k, jumps] = perform_jump(s, t, k, waypoints, regions, jumps, l, P)
%PERFORM_JUMP  Switch CLF target and CBF safe set from mode k to k+1.
y = offset_point(s, l);
hn = wp_region_cbf(y, regions(k + 1), P) - l;
if hn < P.sigma - 1e-9
    error('wp:sight', 'Sight condition failed at jump %d -> %d (h_next = %.3e < sigma = %.3e).', k, k + 1, hn, P.sigma);
end
if ~wp_segment_in_region(y, waypoints(:, k + 1), regions(k + 1), l - 1e-9)
    error('wp:sight', 'Segment from the output to w_%d leaves the safe set of region %d at the jump.', k + 1, k + 1);
end
jumps(end + 1, :) = [t, k, k + 1]; %#ok<AGROW>
k = k + 1;
end

% =========================================================================
function out = postprocess_and_verify(Tc, Sc, Kc, jumps, converged, waypoints, regions, r, l, P)
%POSTPROCESS_AND_VERIFY  Assemble the output and run the a-posteriori checks.
M = size(waypoints, 2);
T = [Tc{:}];
S = [Sc{:}];
modeIdx = [Kc{:}];
N = numel(T);
u = zeros(2, N);
Y = zeros(2, N);
h = zeros(1, N);
bodyMargin = zeros(1, N);
maxCbf = -inf;
maxClf = -inf;
nFallback = 0;
for i = 1:N
    kk = modeIdx(i);
    [u(:, i), aux] = control_law(S(:, i), waypoints(:, kk), regions(kk), l, P);
    Y(:, i) = offset_point(S(:, i), l);
    h(i) = aux.h;
    bodyMargin(i) = wp_region_margin(S(1:2, i), regions(kk));
    nFallback = nFallback + aux.fallback;
    maxCbf = max(maxCbf, aux.cbfViol);
    if ~aux.fallback
        maxClf = max(maxClf, aux.clfViol);
    end
end
ver = struct('minH1', min(h), 'minBodyMargin', min(bodyMargin), ...
    'maxCbfViolation', maxCbf, 'maxClfViolation', maxClf, 'nQpFallback', nFallback, ...
    'finalOutputError', norm(Y(:, end) - waypoints(:, M)), ...
    'nJumps', size(jumps, 1), 'converged', converged);

tolSafe = 1e-9;
if ver.minH1 < -tolSafe
    warning('wp:safety', 'h_eff became negative (min = %.3e).', ver.minH1);
end
if ver.minBodyMargin < -tolSafe
    warning('wp:safety', 'Body position left the active region (min margin = %.3e).', ver.minBodyMargin);
end
if ver.maxCbfViolation > 1e-7 || ver.maxClfViolation > 1e-7
    warning('wp:constraint', 'Constraint violation: CBF %.3e, CLF %.3e.', ver.maxCbfViolation, ver.maxClfViolation);
end
if nFallback > 0
    warning('wp:qpFallback', '%d samples had incompatible CLF/CBF constraints; the CBF was prioritised.', nFallback);
end
if converged && ver.nJumps ~= M - 1
    warning('wp:jumps', 'Number of jumps %d differs from M-1 = %d.', ver.nJumps, M - 1);
end
out = struct('t', T, 's', S, 'y', Y, 'p', S(1:2, :), 'u', u, 'mode', modeIdx, 'h', h, ...
    'bodyMargin', bodyMargin, 'jumps', jumps, 'switchRadius', r, 'verification', ver, ...
    'regions', regions, 'waypoints', waypoints, 'params', P);
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
if nargin < 1
    user = struct();
end
P = struct( ...
    'kappa', 20, ...         % LSE smoothing of the polytope CBF
    'kappa_union', 20, ...   % smooth-max sharpness for a union of 2 polytopes
    'gamma_bar', 1.0, ...    % CLF rate: gamma(s) = 2 gamma_bar s
    'alpha_bar', 2.0, ...    % CBF rate: alpha(s) = alpha_bar s, alpha_bar >= gamma_bar
    'sigma', 0.05, ...       % sight margin h_eff,k+1(y) >= sigma at jumps
    'r_sw', 1.0, ...         % maximum switching radius
    'tol_qp', 1e-12, ...     % gradient-norm tolerance in the closed-form QP
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
if ~(P.sigma > 0 && P.r_sw > 0 && P.kappa > 0 && P.kappa_union > 0)
    error('wp:params', 'sigma, r_sw, kappa, kappa_union must be positive.');
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
        hn = wp_region_cbf(waypoints(:, k), regions(k + 1), P) - rho;
        if hn <= P.sigma
            error('wp:input', ['Waypoint %d is not deep enough in region %d (h_eff = %.3g <= sigma = %.3g). ' ...
                'Regions must overlap around w_%d.'], k, k + 1, hn, P.sigma, k);
        end
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
if numel(regNext.sets) == 1
    rk = rmax;
    return
end
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
