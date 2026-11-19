function P = safe_chain_params(user)
%SAFE_CHAIN_PARAMS  Defaults and validation of the safe-set-chain CLF-CBF
%   trackers (main: double_integrator_safe_chain_clf_cbf; secondary:
%   online_firi/, gcs_bezier/).
%
%   P = safe_chain_params()                 defaults
%   P = safe_chain_params(struct('mu', 0.2)) overrides single fields
%
%   CBF of a safe convex set S = {A x <= b} (unit rows): log-sum-exp
%     h_S(x) = -(1/kappa) ln sum_q exp(kappa (a_q'x - b_q))  <=  min_q (b_q - a_q'x)
%   {h_S >= 0} is convex and contained in S (Matias & Silvestre, Sec. II).
%   CLF V = 1/2 |x - x_hat|^2, gamma(V) = 2 gamma_bar V, alpha(h) = alpha_bar h.
%   alpha_bar >= gamma_bar (Matias & Silvestre Eq. 59) makes the CLF and CBF
%   constraints compatible whenever h_S(x_hat) >= 0: h_S is concave, so its
%   tangent halfspace at x contains {h_S >= 0}, and Eq. 59 applies to that
%   halfspace CBF, whose condition at x equals the condition of h_S at x.

if nargin < 1, user = struct(); end
P = struct( ...
    'kappa', 40, ...          % LSE sharpness (gap to the polytope <= ln(m)/kappa)
    'gamma_bar', 1.0, ...     % CLF rate
    'alpha_bar', 2.0, ...     % CBF rate, >= gamma_bar
    'beta_v', 1.0, ...        % backstepping CLF weight
    'beta_h_nom', 1.0, ...    % nominal backstepping CBF weight
    'gamma1_ratio', 0.9, ...  % gamma1 = ratio * gamma, ratio in (0,1]
    'alpha1_ratio', 1.1, ...  % alpha1 = ratio * alpha, ratio >= 1
    'mu', 0.15, ...           % synergy gap: targets satisfy h >= mu in both adjacent sets
    'sigma', 0.05, ...        % hysteresis: jump into S_{i+1} once h_{i+1}(x) >= sigma (0 < sigma < mu)
    'mu_los', 0.05, ...       % look-ahead: after-next target used if h_i(target) >= mu_los
    'fd_step', 1e-6, ...      % central-difference step for dk/dx
    'tol_qp', 1e-12, ...
    'e_tol', 1e-12, ...
    'tol_final', 1e-3, ...    % |[x - goal; z]| at termination
    't_max', 400, ...
    'max_step', 0.05, ...
    'rel_tol', 1e-8, ...
    'abs_tol', 1e-10);
fn = fieldnames(user);
for i = 1:numel(fn)
    if ~isfield(P, fn{i})
        error('safe_chain:params', 'Unknown parameter "%s".', fn{i});
    end
    P.(fn{i}) = user.(fn{i});
end
if ~(P.gamma_bar > 0 && P.alpha_bar >= P.gamma_bar)
    error('safe_chain:params', 'Require 0 < gamma_bar <= alpha_bar (Matias & Silvestre Eq. 59).');
end
if ~(P.sigma > 0 && P.sigma < P.mu && P.mu_los > 0)
    error('safe_chain:params', 'Require 0 < sigma < mu and mu_los > 0.');
end
if ~(P.gamma1_ratio > 0 && P.gamma1_ratio <= 1 && P.alpha1_ratio >= 1)
    error('safe_chain:params', 'Require gamma1_ratio in (0,1] and alpha1_ratio >= 1.');
end
if ~(P.kappa > 0 && P.beta_v > 0 && P.beta_h_nom > 0 && P.fd_step > 0)
    error('safe_chain:params', 'kappa, beta_v, beta_h_nom and fd_step must be positive.');
end
end
