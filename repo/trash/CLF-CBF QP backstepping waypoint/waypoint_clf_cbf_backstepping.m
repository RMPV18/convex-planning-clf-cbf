function out = waypoint_clf_cbf_backstepping(waypoints, regions, model, s0, params)
%WAYPOINT_CLF_CBF_BACKSTEPPING  Hybrid CLF-CBF QP waypoint follower.
%   The top-level (first-order) input is the closed-form CLF-CBF QP; for
%   second-order models it is the virtual control of a joint CLF-CBF
%   backstepping. The CLF target and the CBF safe set switch together.
%
%   out = waypoint_clf_cbf_backstepping(waypoints, regions, model, s0)
%   out = waypoint_clf_cbf_backstepping(waypoints, regions, model, s0, params)
%
% INPUTS
%   waypoints  n x M   waypoints w_1..w_M (w_M = final desired point).
%   regions    1 x M   struct array, fields A (m_k x n), b (m_k x 1): convex
%                      polytope Q_k = {y : A y <= b} active while the target is
%                      w_k (IRIS / IRIS-NP / FIRI / GCS region sequence).
%   model      struct  see WAYPOINT_MODELS (single/double integrator, damped
%                      second order, unicycle, or custom handles).
%   s0         full initial state of the model (column).
%   params     struct  optional overrides, see FILL_PARAMS.
%
% Requirements (checked, error on violation), with h_eff = h_poly - model.rho:
%   h_eff,k(w_k) > 0,  h_eff,k+1(w_k) > sigma,  h_eff,1(y(s0)) > 0.
%
% TOP LEVEL: CLF-CBF QP in closed form (Matias & Silvestre Eq. 10-15, p = inf)
%   y = x:  xdot = f0 + G0 v     (v = u for first-order models, v = z otherwise)
%   V = 1/2 |y - w_k|^2,    h = h_poly,k(y) - rho,
%   h_poly(y) = -(1/kappa) ln sum_i exp(kappa (a_i'y - b_i))  (concave, 1-Lipschitz)
%   k(s) = argmin 1/2|v|^2  s.t.  LfV + LGV v <= -2 gamma_bar V
%                                 Lfh + LGh v >= -alpha_bar h
%   evaluated by CLF_CBF_QP_CLOSED_FORM. Compatibility for G0 = I is Eq. 59
%   (alpha_bar >= gamma_bar). For general G0(s) incompatibility (aligned LGV, LGh)
%   is detected, reported with warning('waypoint_clf_cbf:qpFallback') and resolved
%   in favour of the CBF.
%
% SECOND-ORDER MODELS (joint CLF-CBF backstepping, Eq. 72-77, strict feedback)
%       xdot = f0 + G0 z,   zdot = f1 + G1 u
%   e = z - k(s),  V1 = V + |e|^2/(2 beta_v),  h1 = h - |e|^2/(2 beta_h),
%   kdot = (dk/dx)(f0 + G0 z)   (dk/dx by central differences: the QP controller
%   is locally Lipschitz and piecewise smooth, not C^1 at regime boundaries).
%   With a = G1'e, both constraints become  a'u <= B_V  and  a'u <= B_h:
%       B_V = -beta_v [ gV'(f0+G0 z) + gamma1(V1) ] - e'(f1 - kdot)
%       B_h =  beta_h [ gh'(f0+G0 z) + alpha1(h1) ] - e'(f1 - kdot)
%   gamma1 = gamma1_ratio * gamma (< gamma), alpha1 = alpha1_ratio * alpha (>= alpha)
%   so that both are feasible at e = 0. Closed-form QP min 1/2|u|^2 s.t. a'u <= B:
%       u = 0 if B = min(B_V, B_h) >= 0,  else  u = (B/|a|^2) a.
%
% SWITCHING (target and safe set change together, k -> k+1)
%   Jump when |y - w_k| <= r_k,  r_k = min(r_sw, h_eff,k+1(w_k) - sigma).
%   Sight: h_eff,k+1(y) >= sigma. {h_eff,k+1 >= 0} is convex and contains
%   w_{k+1}, hence the segment y -> w_{k+1} lies in the new safe set. The
%   1-Lipschitz property of h_poly guarantees the condition at the jump; it is
%   also verified at runtime. For second-order models beta_h is reset by
%       beta_h+ = max(beta_h_nom, 1.05 |e+|^2 / (2 h_eff,k+1(y)))   (Eq. 84)
%   which preserves h1 >= 0.
%
% OUTPUT struct
%   t, s, y, u, mode        time, full state, controlled output, input, index k
%   h, h1, betaH            CBF h_eff, backstepping CBF (= h for first order)
%   bodyMargin              min_i (b_i - a_i'p) of the active polytope, with
%                           p = model.position(s) (physical reference point)
%   jumps                   J x 4 [time, from, to, beta_h after jump]
%   switchRadius, verification, regions, waypoints, params, model
%   verification: minH1, minBodyMargin, maxCbfViolation, maxClfViolation,
%                 maxBackstepViolation, nQpFallback, finalOutputError,
%                 finalResidual, nJumps, converged.
%   Violations raise warning('waypoint_clf_cbf:...') (never suppressed).
%
% REFERENCES
%   Matias & Silvestre, arXiv:2504.09760v2: Sec. II-B, III, IV, V-D.
%   Molnar & Ames (smooth-max CBF; [36] therein).

if nargin < 5
    params = struct();
end
P = fill_params(params);
validate_model(model);

n = model.n;
M = size(waypoints, 2);
if size(waypoints, 1) ~= n
    error('waypoint_clf_cbf:input', 'waypoints must be %d x M.', n);
end
if numel(regions) ~= M
    error('waypoint_clf_cbf:input', 'regions must have one entry per waypoint (%d expected, %d given).', M, numel(regions));
end
s0 = s0(:);

% ---- normalise polytope rows (unit normals -> h_poly is 1-Lipschitz) -----
for k = 1:M
    A = regions(k).A;
    b = regions(k).b(:);
    if size(A, 2) ~= n || size(A, 1) ~= numel(b)
        error('waypoint_clf_cbf:input', 'regions(%d): A must be m x %d and b m x 1.', k, n);
    end
    nr = sqrt(sum(A.^2, 2));
    if any(nr <= 0)
        error('waypoint_clf_cbf:input', 'regions(%d) has a zero row in A.', k);
    end
    regions(k).A = A ./ nr;
    regions(k).b = b ./ nr;
end

% ---- consistency of the region / waypoint chain --------------------------
r = zeros(1, max(M - 1, 0));
for k = 1:M
    hk = lse_cbf(waypoints(:, k), regions(k).A, regions(k).b, P.kappa) - model.rho;
    if hk <= 0
        error('waypoint_clf_cbf:input', ...
            'Waypoint %d is not in the interior of its safe set (h_eff = %.3g <= 0). Increase kappa, reduce rho or move it inward.', k, hk);
    end
    if k < M
        hn = lse_cbf(waypoints(:, k), regions(k + 1).A, regions(k + 1).b, P.kappa) - model.rho;
        if hn <= P.sigma
            error('waypoint_clf_cbf:input', ...
                'Waypoint %d is not deep enough in region %d (h_eff = %.3g <= sigma = %.3g). Regions must overlap around w_%d.', k, k + 1, hn, P.sigma, k);
        end
        r(k) = min(P.r_sw, hn - P.sigma);
    end
end
y0 = model.output(s0);
h0 = lse_cbf(y0, regions(1).A, regions(1).b, P.kappa) - model.rho;
if h0 <= 0
    error('waypoint_clf_cbf:input', 'Initial output is not in the interior of the first safe set (h_eff = %.3g).', h0);
end

% ---- hybrid simulation ----------------------------------------------------
k = 1;
t = 0;
s = s0;
betaH = update_beta_h(s, waypoints(:, k), regions(k), model, P);
jumps = zeros(0, 4);
Tc = {};
Sc = {};
Kc = {};
Bc = {};
converged = false;

while true
    while k < M && norm(model.output(s) - waypoints(:, k)) <= r(k)
        [k, betaH, jumps] = perform_jump(s, t, k, betaH, waypoints, regions, jumps, model, P);
    end
    if k == M && final_residual(s, waypoints(:, M), model) <= P.tol_final
        Tc{end + 1} = t; %#ok<AGROW>
        Sc{end + 1} = s; %#ok<AGROW>
        Kc{end + 1} = k; %#ok<AGROW>
        Bc{end + 1} = betaH; %#ok<AGROW>
        converged = true;
        break
    end

    w = waypoints(:, k);
    reg = regions(k);
    kk = k;
    bh = betaH;
    opts = odeset('RelTol', P.rel_tol, 'AbsTol', P.abs_tol, 'MaxStep', P.max_step, ...
        'Events', @(tt, ss) event_fcn(ss, w, kk, M, r, model, P));
    [ts, ss, te, se, ie] = ode45(@(tt, ss) closed_loop(ss, w, reg, bh, model, P), [t, P.t_max], s.', opts);

    Tc{end + 1} = ts(:).'; %#ok<AGROW>
    Sc{end + 1} = ss.'; %#ok<AGROW>
    Kc{end + 1} = k * ones(1, numel(ts)); %#ok<AGROW>
    Bc{end + 1} = betaH * ones(1, numel(ts)); %#ok<AGROW>

    if isempty(ie)
        warning('waypoint_clf_cbf:timeout', ...
            'Simulation reached t_max = %g in mode %d without triggering the switching/convergence event.', P.t_max, k);
        break
    end
    if k == M
        converged = true;
        break
    end
    t = te(end);
    s = se(end, :).';
    [k, betaH, jumps] = perform_jump(s, t, k, betaH, waypoints, regions, jumps, model, P);
end

T = [Tc{:}];
S = [Sc{:}];
modeIdx = [Kc{:}];
betaHs = [Bc{:}];
N = numel(T);

% ---- post-processing and verification ------------------------------------
u = zeros(model.m, N);
Y = zeros(n, N);
h = zeros(1, N);
h1 = zeros(1, N);
bodyMargin = zeros(1, N);
maxCbf = -inf;
maxClf = -inf;
maxBs = -inf;
nFallback = 0;
for i = 1:N
    kk = modeIdx(i);
    [u(:, i), aux] = control_law(S(:, i), waypoints(:, kk), regions(kk), betaHs(i), model, P);
    Y(:, i) = model.output(S(:, i));
    h(i) = aux.h;
    h1(i) = aux.h1;
    p = model.position(S(:, i));
    bodyMargin(i) = min(regions(kk).b - regions(kk).A * p);
    nFallback = nFallback + aux.fallback;
    if strcmp(model.type, 'first_order')
        maxCbf = max(maxCbf, aux.cbfViol);
        if ~aux.fallback
            maxClf = max(maxClf, aux.clfViol);
        end
    else
        maxBs = max(maxBs, aux.bsViol);
    end
end

ver = struct();
ver.minH1 = min(h1);
ver.minBodyMargin = min(bodyMargin);
ver.maxCbfViolation = maxCbf;
ver.maxClfViolation = maxClf;
ver.maxBackstepViolation = maxBs;
ver.nQpFallback = nFallback;
ver.finalOutputError = norm(Y(:, end) - waypoints(:, M));
ver.finalResidual = final_residual(S(:, end), waypoints(:, M), model);
ver.nJumps = size(jumps, 1);
ver.converged = converged;

tolSafe = 1e-9;
if ver.minH1 < -tolSafe
    warning('waypoint_clf_cbf:safety', 'h1/h_eff became negative (min = %.3e).', ver.minH1);
end
if ver.minBodyMargin < -tolSafe
    warning('waypoint_clf_cbf:safety', 'Body reference point left the active polytope (min margin = %.3e).', ver.minBodyMargin);
end
if ver.maxCbfViolation > 1e-7 || ver.maxBackstepViolation > 1e-7 || ver.maxClfViolation > 1e-7
    warning('waypoint_clf_cbf:constraint', 'Constraint violation: CBF %.3e, CLF %.3e, backstepping %.3e.', ...
        ver.maxCbfViolation, ver.maxClfViolation, ver.maxBackstepViolation);
end
if nFallback > 0
    warning('waypoint_clf_cbf:qpFallback', ...
        '%d samples had incompatible CLF/CBF constraints; the CBF was prioritised.', nFallback);
end
if converged && ver.nJumps ~= M - 1
    warning('waypoint_clf_cbf:jumps', 'Number of jumps %d differs from M-1 = %d.', ver.nJumps, M - 1);
end

out = struct('t', T, 's', S, 'y', Y, 'u', u, 'mode', modeIdx, 'h', h, 'h1', h1, ...
    'betaH', betaHs, 'bodyMargin', bodyMargin, 'jumps', jumps, 'switchRadius', r, ...
    'verification', ver, 'regions', regions, 'waypoints', waypoints, ...
    'params', P, 'model', model);
end

% =========================================================================
function P = fill_params(user)
%FILL_PARAMS  Defaults and validation of controller parameters.
P = struct( ...
    'kappa', 20, ...         % LSE smoothing of the polytope CBF
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
    'e_tol', 1e-12, ...      % |a| below which the backstepping input is 0
    'tol_final', 1e-4, ...   % termination residual in the last mode
    't_max', 200, ...        % absolute time limit
    'max_step', 0.05, ...    % ode45 maximum step
    'rel_tol', 1e-9, ...
    'abs_tol', 1e-11);
fn = fieldnames(user);
for i = 1:numel(fn)
    if ~isfield(P, fn{i})
        error('waypoint_clf_cbf:params', 'Unknown parameter "%s".', fn{i});
    end
    P.(fn{i}) = user.(fn{i});
end
if ~(P.gamma_bar > 0 && P.alpha_bar >= P.gamma_bar)
    error('waypoint_clf_cbf:params', 'Require 0 < gamma_bar <= alpha_bar (Eq. 59).');
end
if ~(P.gamma1_ratio > 0 && P.gamma1_ratio <= 1 && P.alpha1_ratio >= 1)
    error('waypoint_clf_cbf:params', 'Require gamma1_ratio in (0,1] and alpha1_ratio >= 1.');
end
if ~(P.sigma > 0 && P.r_sw > 0 && P.kappa > 0 && P.beta_v > 0 && P.beta_h_nom > 0 && P.fd_step > 0)
    error('waypoint_clf_cbf:params', 'sigma, r_sw, kappa, beta_v, beta_h_nom, fd_step must be positive.');
end
end

% =========================================================================
function validate_model(model)
req = {'type', 'n', 'm', 'output', 'position', 'f0', 'G0', 'dyn', 'rho'};
for i = 1:numel(req)
    if ~isfield(model, req{i})
        error('waypoint_clf_cbf:model', 'Model field "%s" is missing.', req{i});
    end
end
switch model.type
    case 'first_order'
    case 'second_order'
        req2 = {'xIdx', 'zIdx', 'f1', 'G1'};
        for i = 1:numel(req2)
            if ~isfield(model, req2{i})
                error('waypoint_clf_cbf:model', 'Second-order model field "%s" is missing.', req2{i});
            end
        end
        if numel(model.xIdx) ~= model.n || numel(model.zIdx) ~= model.n
            error('waypoint_clf_cbf:model', 'xIdx and zIdx must have n = %d entries.', model.n);
        end
    otherwise
        error('waypoint_clf_cbf:model', 'model.type must be first_order or second_order.');
end
if model.rho < 0
    error('waypoint_clf_cbf:model', 'model.rho must be nonnegative.');
end
end

% =========================================================================
function [h, gh] = lse_cbf(x, A, b, kappa)
%LSE_CBF  h(x) = -(1/kappa) ln sum_i exp(kappa (a_i'x - b_i)); concave; h>=0 => Ax<=b.
g = A * x - b;
gmax = max(g);
ex = exp(kappa * (g - gmax));
S = sum(ex);
h = -(gmax + log(S) / kappa);
gh = -A.' * (ex / S);
end

% =========================================================================
function [v, aux] = top_level_qp(s, w, reg, model, P)
%TOP_LEVEL_QP  Closed-form CLF-CBF QP for ydot = f0(s) + G0(s) v.
y = model.output(s);
f = model.f0(s);
G = model.G0(s);
[hp, gh] = lse_cbf(y, reg.A, reg.b, P.kappa);
h = hp - model.rho;
gV = y - w;
V = 0.5 * (gV.' * gV);
LfV = gV.' * f;
LGV = gV.' * G;
Lfh = gh.' * f;
LGh = gh.' * G;
gammaV = 2 * P.gamma_bar * V;
alphah = P.alpha_bar * h;
[v, regime, feasible] = clf_cbf_qp_closed_form(LfV, LGV, gammaV, Lfh, LGh, alphah, P.tol_qp);
aux = struct('h', h, 'V', V, 'gV', gV, 'gh', gh, 'f', f, 'G', G, 'regime', regime, ...
    'feasible', feasible, ...
    'clfViol', LfV + LGV * v + gammaV, ...
    'cbfViol', -(Lfh + LGh * v + alphah));
end

% =========================================================================
function [u, aux] = control_law(s, w, reg, betaH, model, P)
%CONTROL_LAW  Input u for the active mode (first order or backstepping).
n = model.n;
[v, tq] = top_level_qp(s, w, reg, model, P);
if strcmp(model.type, 'first_order')
    u = v;
    aux = struct('h', tq.h, 'h1', tq.h, 'fallback', ~tq.feasible, ...
        'clfViol', tq.clfViol, 'cbfViol', tq.cbfViol, 'bsViol', -inf);
    return
end

xi = model.xIdx;
zi = model.zIdx;
z = s(zi);
% numerical Jacobian dk/dx of the (piecewise smooth) top-level QP controller
J = zeros(n, n);
for j = 1:n
    d = P.fd_step * (1 + abs(s(xi(j))));
    sp = s;
    sm = s;
    sp(xi(j)) = sp(xi(j)) + d;
    sm(xi(j)) = sm(xi(j)) - d;
    J(:, j) = (top_level_qp(sp, w, reg, model, P) - top_level_qp(sm, w, reg, model, P)) / (2 * d);
end
xdot = tq.f + tq.G * z;
kdot = J * xdot;
e = z - v;
ee = e.' * e;
V1 = tq.V + ee / (2 * P.beta_v);
h1 = tq.h - ee / (2 * betaH);
f1 = model.f1(s);
G1 = model.G1(s);
a = G1.' * e;
rr = -e.' * (f1 - kdot);
BV = -P.beta_v * (tq.gV.' * xdot + 2 * P.gamma_bar * P.gamma1_ratio * V1) + rr;
Bh = betaH * (tq.gh.' * xdot + P.alpha_bar * P.alpha1_ratio * h1) + rr;
B = min(BV, Bh);
if B >= 0 || (a.' * a) < P.e_tol^2
    u = zeros(model.m, 1);
else
    u = (B / (a.' * a)) * a;
end
aux = struct('h', tq.h, 'h1', h1, 'fallback', ~tq.feasible, ...
    'clfViol', -inf, 'cbfViol', -inf, 'bsViol', a.' * u - B);
end

% =========================================================================
function ds = closed_loop(s, w, reg, betaH, model, P)
u = control_law(s, w, reg, betaH, model, P);
ds = model.dyn(s, u);
ds = ds(:);
end

% =========================================================================
function res = final_residual(s, w, model)
y = model.output(s);
if strcmp(model.type, 'second_order')
    res = norm([y - w; s(model.zIdx)]);
else
    res = norm(y - w);
end
end

% =========================================================================
function [val, term, dir] = event_fcn(s, w, k, M, r, model, P)
%EVENT_FCN  Jump event (k < M) or convergence event (k = M).
if k < M
    val = norm(model.output(s) - w) - r(k);
else
    val = final_residual(s, w, model) - P.tol_final;
end
term = 1;
dir = -1;
end

% =========================================================================
function betaH = update_beta_h(s, w, reg, model, P)
%UPDATE_BETA_H  Eq. (84): guarantees h1 >= 0 at the mode update (second order).
if strcmp(model.type, 'first_order')
    betaH = P.beta_h_nom;
    return
end
[v, tq] = top_level_qp(s, w, reg, model, P);
if tq.h <= 0
    error('waypoint_clf_cbf:safety', 'Cannot set beta_h: h_eff = %.3e <= 0 at the mode update.', tq.h);
end
e = s(model.zIdx) - v;
betaH = max(P.beta_h_nom, 1.05 * (e.' * e) / (2 * tq.h));
end

% =========================================================================
function [k, betaH, jumps] = perform_jump(s, t, k, betaH, waypoints, regions, jumps, model, P)
%PERFORM_JUMP  Switch CLF target and CBF safe set from mode k to k+1.
hn = lse_cbf(model.output(s), regions(k + 1).A, regions(k + 1).b, P.kappa) - model.rho;
if hn < P.sigma - 1e-9    % tolerance: event location is accurate to ode45 tolerances
    error('waypoint_clf_cbf:sight', ...
        'Sight condition failed at jump %d -> %d (h_next = %.3e < sigma = %.3e).', k, k + 1, hn, P.sigma);
end
from = k;
k = k + 1;
betaH = update_beta_h(s, waypoints(:, k), regions(k), model, P);
jumps(end + 1, :) = [t, from, k, betaH]; %#ok<AGROW>
end
