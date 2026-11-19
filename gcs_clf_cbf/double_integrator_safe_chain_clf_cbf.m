function out = double_integrator_safe_chain_clf_cbf(sets, targets, x0, z0, params)
%DOUBLE_INTEGRATOR_SAFE_CHAIN_CLF_CBF  Hybrid CLF-CBF backstepping tracker of a
%   double integrator through a chain of bounded convex safe sets.
%
%   out = double_integrator_safe_chain_clf_cbf(sets, targets, x0, z0)
%   out = double_integrator_safe_chain_clf_cbf(sets, targets, x0, z0, params)
%
%   sets     1 x N cell of struct A, b: safe convex sets S_1..S_N {A x <= b}
%            (bounded; consecutive sets overlap; build_safe_chain)
%   targets  n x N: targets(:,i) in S_i ∩ S_{i+1} with h_i, h_{i+1} >= mu
%            (translate_waypoint); targets(:,N) = goal, in S_N
%   x0, z0   initial position (in S_1) and velocity
%   params   overrides of safe_chain_params
%
%   Adaptation of Matias & Silvestre (arXiv:2504.09760v2) from the unbounded
%   safe set outside one polytope to a chain of bounded convex safe sets:
%   MODE i = (safe set S_i, active target x_hat).
%     Flow: safe_chain_control with CLF about x_hat and the LSE CBF of S_i
%     (alpha_bar >= gamma_bar: CLF and CBF never conflict, x_hat in S_i).
%     Look-ahead ("after the next waypoint"): x_hat = targets(:,i+1) if it is
%     in sight, i.e. h_i(targets(:,i+1)) >= mu_los (S_i convex: the segment
%     x -> x_hat stays in S_i); otherwise x_hat = targets(:,i).
%     Jump (set-entry trigger, hysteresis as Eq. 36): when h_{i+1}(x) >=
%     sigma, i <- i+1 and the target is re-selected. Since x_hat lies in
%     S_{i+1} with margin mu > sigma and x -> x_hat exponentially, every jump
%     occurs in finite time; i only increases, so there are at most N-1 jumps
%     and the mode whose target is the goal converges to it (this can be
%     mode N-1 when the look-ahead already selects the goal).
%     beta_h+ = max(beta_h_nom, 1.05 |e+|^2 / (2 h_{i+1}(x)))  (Eq. 84)
%     keeps h1 >= 0 across jumps (safety of the second-order system).
%
%   out: t, s, x, z, u, mode, target (n x samples), h, h1, betaH, margin
%   (polytope margin of x in the active set), jumps [t from to betaH],
%   verification (minH1, minMargin, maxBackstepViolation, nQpFallback,
%   finalError, nJumps, converged), wallTime (simulation wall-clock time).

if nargin < 5, params = struct(); end
P = safe_chain_params(params);
t0 = tic;
x0 = x0(:);  z0 = z0(:);  n = numel(x0);
N = numel(sets);
assert(size(targets, 2) == N && size(targets, 1) == n, 'safe_chain:input', ...
    'targets must be %d x %d.', n, N);
sets = normalise_sets(sets);
check_chain(sets, targets, x0, P);

i = 1;  t = 0;  s = [x0; z0];
xh = select_target(i, sets, targets, P);
betaH = update_beta_h(s, xh, sets{i}, P);
jumps = zeros(0, 4);
Tc = {};  Sc = {};  Ic = {};  Bc = {};  Xc = {};
converged = false;
while true
    while i < N && lse_cbf(s(1:n), sets{i+1}.A, sets{i+1}.b, P.kappa) >= P.sigma - 1e-9
        [i, xh, betaH, jumps] = do_jump(s, t, i, sets, targets, jumps, P);
    end
    if isequal(xh, targets(:, N)) && norm([s(1:n) - targets(:, N); s(n+1:end)]) <= P.tol_final
        converged = true;  break
    end
    ii = i;  xx = xh;  bb = betaH;
    opts = odeset('RelTol', P.rel_tol, 'AbsTol', P.abs_tol, 'MaxStep', P.max_step, ...
        'Events', @(tt, ss) events(ss, ii, xx, sets, targets, N, n, P));
    [ts, ss, te, se, ie] = ode45(@(tt, ss) [ss(n+1:2*n); safe_chain_control(ss, xx, sets{ii}, bb, P)], ...
        [t, P.t_max], s', opts);
    Tc{end+1} = ts';  Sc{end+1} = ss';  Ic{end+1} = i * ones(1, numel(ts)); %#ok<AGROW>
    Bc{end+1} = betaH * ones(1, numel(ts));  Xc{end+1} = repmat(xh, 1, numel(ts)); %#ok<AGROW>
    if isempty(ie)
        warning('safe_chain:timeout', 'Reached t_max = %g in mode %d without a jump or convergence.', P.t_max, i);
        break
    end
    t = te(end);  s = se(end, :)';
    if ie(end) == 2, converged = true; break; end
end
out = postprocess(Tc, Sc, Ic, Bc, Xc, jumps, converged, sets, targets, n, P);
out.wallTime = toc(t0);
end

% =========================================================================
function sets = normalise_sets(sets)
for k = 1:numel(sets)
    nr = vecnorm(sets{k}.A, 2, 2);
    sets{k} = struct('A', sets{k}.A ./ nr, 'b', sets{k}.b(:) ./ nr);
end
end

function check_chain(sets, targets, x0, P)
N = numel(sets);
h0 = lse_cbf(x0, sets{1}.A, sets{1}.b, P.kappa);
if h0 <= 0
    error('safe_chain:input', 'x0 is not inside the first safe set (h = %.3g).', h0);
end
for k = 1:N
    hk = lse_cbf(targets(:, k), sets{k}.A, sets{k}.b, P.kappa);
    if hk < 0
        error('safe_chain:input', 'Target %d is outside safe set %d (h = %.3g).', k, k, hk);
    end
    if k < N
        hn = lse_cbf(targets(:, k), sets{k+1}.A, sets{k+1}.b, P.kappa);
        if hn <= P.sigma
            error('safe_chain:input', ['Target %d is not inside safe set %d with margin sigma ' ...
                '(h = %.3g <= %.3g): the jump %d -> %d could never trigger.'], k, k+1, hn, P.sigma, k, k+1);
        end
    end
end
end

function xh = select_target(i, sets, targets, P)
N = numel(sets);
xh = targets(:, i);
if i < N && lse_cbf(targets(:, i+1), sets{i}.A, sets{i}.b, P.kappa) >= P.mu_los
    xh = targets(:, i+1);      % after-next target in sight
end
end

function betaH = update_beta_h(s, xh, S, P)
n = numel(xh);
[~, aux] = safe_chain_control(s, xh, S, P.beta_h_nom, P);
if aux.h <= 0
    error('safe_chain:safety', 'Cannot set beta_h: h = %.3e <= 0 at the mode update.', aux.h);
end
e = s(n+1:end) - aux.k;
betaH = max(P.beta_h_nom, 1.05 * (e' * e) / (2 * aux.h));
end

function [i, xh, betaH, jumps] = do_jump(s, t, i, sets, targets, jumps, P)
from = i;  i = i + 1;
xh = select_target(i, sets, targets, P);
betaH = update_beta_h(s, xh, sets{i}, P);
jumps(end+1, :) = [t, from, i, betaH];
end

function [val, term, dir] = events(s, i, xh, sets, targets, N, n, P)
% 1: entry into S_{i+1} (jump);  2: convergence, active whenever the target
% is the goal (the look-ahead may select it before the last set is entered).
val = [1; 1];  dir = [1; -1];  term = [1; 1];
if i < N
    val(1) = lse_cbf(s(1:n), sets{i+1}.A, sets{i+1}.b, P.kappa) - P.sigma;
end
if isequal(xh, targets(:, N))
    val(2) = norm([s(1:n) - targets(:, N); s(n+1:2*n)]) - P.tol_final;
end
end

function out = postprocess(Tc, Sc, Ic, Bc, Xc, jumps, converged, sets, targets, n, P)
T = [Tc{:}];  S = [Sc{:}];  I = [Ic{:}];  Bh = [Bc{:}];  X = [Xc{:}];
K = numel(T);
u = zeros(n, K);  h = zeros(1, K);  h1 = zeros(1, K);  mg = zeros(1, K);
maxBs = -inf;  nFb = 0;
for k = 1:K
    [u(:, k), aux] = safe_chain_control(S(:, k), X(:, k), sets{I(k)}, Bh(k), P);
    h(k) = aux.h;  h1(k) = aux.h1;
    mg(k) = min(sets{I(k)}.b - sets{I(k)}.A * S(1:n, k));
    nFb = nFb + aux.fallback;  maxBs = max(maxBs, aux.bsViol);
end
ver = struct('minH1', min(h1), 'minMargin', min(mg), 'maxBackstepViolation', maxBs, ...
    'nQpFallback', nFb, 'finalError', norm(S(1:n, end) - targets(:, end)), ...
    'nJumps', size(jumps, 1), 'converged', converged);
if ver.minH1 < -1e-9
    warning('safe_chain:safety', 'h1 became negative (min %.3e).', ver.minH1);
end
if ver.minMargin < -1e-9
    warning('safe_chain:safety', 'The state left the active safe set (min margin %.3e).', ver.minMargin);
end
if ver.maxBackstepViolation > 1e-7
    warning('safe_chain:constraint', 'Backstepping constraint violated by %.3e.', ver.maxBackstepViolation);
end
if nFb > 0
    warning('safe_chain:qpFallback', '%d samples with incompatible CLF/CBF constraints.', nFb);
end
out = struct('t', T, 's', S, 'x', S(1:n, :), 'z', S(n+1:end, :), 'u', u, 'mode', I, ...
    'target', X, 'h', h, 'h1', h1, 'betaH', Bh, 'margin', mg, 'jumps', jumps, ...
    'verification', ver, 'params', P);
end
