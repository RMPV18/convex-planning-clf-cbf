function [u, aux] = safe_chain_control(s, target, S, betaH, P)
%SAFE_CHAIN_CONTROL  Joint CLF-CBF backstepping controller of a double
%   integrator inside one bounded convex safe set (shared by the main and
%   the secondary safe-chain trackers).
%
%   [u, aux] = safe_chain_control(s, target, S, betaH, P)
%
%   s       [x; z] (2n x 1), plant xdot = z, zdot = u
%   target  n x 1 active target x_hat, h_S(x_hat) >= 0 (compatibility)
%   S       struct A, b (unit rows): safe convex set {A x <= b}
%   betaH   backstepping CBF weight of the current mode
%   P       safe_chain_params
%
%   TOP LEVEL (xdot = v): closed-form CLF-CBF QP (Matias & Silvestre Eq. 60,
%   closed form Eqs. 10-15 with p = inf):
%     k(x) = argmin 1/2|v|^2  s.t.  (x - x_hat)'v <= -2 gamma_bar V,
%                                    grad h_S'v   >= -alpha_bar h_S
%   with the log-sum-exp CBF h_S (safe_chain_params). Always feasible for
%   alpha_bar >= gamma_bar and h_S(x_hat) >= 0; aux.fallback flags a violation.
%   BACKSTEPPING (Matias & Silvestre Eqs. 72-77 = Taylor et al. Sec. V):
%     e = z - k(x),  V1 = V + |e|^2/(2 beta_v),  h1 = h_S - |e|^2/(2 betaH)
%     constraints e'u <= B_V and e'u <= B_h (same direction, always jointly
%     feasible); u = argmin 1/2|u|^2 s.t. e'u <= min(B_V, B_h).
%     kdot = (dk/dx) z by central differences (k is Lipschitz, piecewise
%     smooth; the C^1 requirement of the papers holds away from regime
%     switches of the top-level QP).

n = numel(target);
x = s(1:n);  z = s(n+1:2*n);
[kx, tq] = top_level(x, target, S, P);
J = zeros(n);
for j = 1:n
    d = P.fd_step * (1 + abs(x(j)));
    xp = x;  xm = x;  xp(j) = xp(j) + d;  xm(j) = xm(j) - d;
    J(:, j) = (top_level(xp, target, S, P) - top_level(xm, target, S, P)) / (2 * d);
end
kdot = J * z;
e = z - kx;  ee = e' * e;
V1 = tq.V + ee / (2 * P.beta_v);
h1 = tq.h - ee / (2 * betaH);
BV = -P.beta_v * (tq.gV' * z + 2 * P.gamma_bar * P.gamma1_ratio * V1) + e' * kdot;
Bh = betaH * (tq.gh' * z + P.alpha_bar * P.alpha1_ratio * h1) + e' * kdot;
B = min(BV, Bh);
if B >= 0 || ee < P.e_tol^2
    u = zeros(n, 1);
else
    u = (B / ee) * e;
end
aux = struct('h', tq.h, 'h1', h1, 'V', tq.V, 'V1', V1, 'k', kx, ...
    'fallback', ~tq.feasible, 'bsViol', e' * u - B);
end

% =========================================================================
function [v, aux] = top_level(x, target, S, P)
[h, gh] = lse_cbf(x, S.A, S.b, P.kappa);
gV = x - target;
V = 0.5 * (gV' * gV);
[v, feasible] = qp_closed_form(gV, 2 * P.gamma_bar * V, gh, P.alpha_bar * h, P.tol_qp);
aux = struct('h', h, 'V', V, 'gV', gV, 'gh', gh, 'feasible', feasible);
end

function [u, feasible] = qp_closed_form(a, FV, g, Fh, tol)
% min 1/2|u|^2 s.t. a'u <= -FV (CLF), g'u >= -Fh (CBF): active-set closed
% form S4, S1, S2, S3 (Matias & Silvestre Eqs. 10-15, p = inf). On
% incompatibility the CBF-only solution is returned (feasible = false).
feasible = true;
na = a' * a;  ng = g' * g;  q = numel(a);
if na <= tol^2, FV = min(FV, 0); end     % x at the target: V = 0, CLF row is 0 <= 0
if FV <= 0 && Fh >= 0, u = zeros(q, 1); return; end
u2 = [];
if na > tol^2 && FV > 0
    u1 = -(FV / na) * a;
    if g' * u1 >= -Fh - 1e-12, u = u1; return; end
end
if ng > tol^2 && Fh < 0
    u2 = -(Fh / ng) * g;
    if a' * u2 <= -FV + 1e-12, u = u2; return; end
end
L = a' * g;  Delta = L^2 - na * ng;
if abs(Delta) > 1e-12 * max(1, na * ng)
    l1 = (L * Fh - ng * FV) / Delta;
    l2 = (na * Fh - L * FV) / Delta;
    if l1 >= -1e-9 && l2 >= -1e-9
        u = -l1 * a + l2 * g;  return
    end
end
feasible = false;
if isempty(u2), u2 = zeros(q, 1); end
u = u2;
end
