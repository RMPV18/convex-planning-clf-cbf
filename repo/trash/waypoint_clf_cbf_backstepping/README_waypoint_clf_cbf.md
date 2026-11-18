# Waypoint-following closed-form CLF-CBF QP with backstepping

Files
- `clf_cbf_qp_closed_form.m` : exact closed-form QP for one CLF + one CBF (Matias & Silvestre Eq. 10-15, p = inf).
- `waypoint_models.m` : single integrator, double integrator, damped second order, unicycle; template for custom models.
- `waypoint_clf_cbf_backstepping.m` : hybrid controller + simulation (ode45 events).
- `demo_waypoint_clf_cbf.m` : corridor with 3 regions and 3 waypoints on three plants; quadprog cross-check.

Requirements: MATLAB R2018b+ (`yline`, `tiledlayout`). Core files need base MATLAB only; the demo cross-check uses the Optimization Toolbox if licensed.

## Call
`out = waypoint_clf_cbf_backstepping(W, regions, model, s0, params)`
- `W` n x M waypoints; `regions(k).A, .b` polytope active while targeting `W(:,k)`; `s0` full initial state.

## Top level (every model)
`ydot = f0(s) + G0(s) v`; `V = 1/2 |y - w_k|^2`; `h = h_LSE,k(y) - rho`.
`v = argmin 1/2|v|^2  s.t.  LfV + LGV v <= -2 gamma_bar V,  Lfh + LGh v >= -alpha_bar h`, solved by `clf_cbf_qp_closed_form` (regimes S1..S4; incompatible constraints are detected, CBF prioritised, warning raised).

## Model classes
- `first_order`: `u = v`. Covers single integrator and any fully actuated first-order system. Unicycle: control the offset point `y = p + l(cos th, sin th)`, `ydot = R(th) diag(1,l) u`, `rho = l` protects the body.
- `second_order`: `xdot = f0 + G0 z`, `zdot = f1 + G1 u`. `v` is the virtual control. Joint backstepping (`e = z - v`, `V1`, `h1`) gives `a'u <= min(B_V, B_h)`, `a = G1'e`; closed form `u = (B/|a|^2) a` if `B < 0`, else `0`. `dk/dx` by central differences.
- Custom model: supply the fields documented in `waypoint_models.m`.

## Switching
- Jump at `|y - w_k| <= r_k = min(r_sw, h_eff,k+1(w_k) - sigma)`.
- Sight: `h_eff,k+1(y) >= sigma`; convexity places the segment `y -> w_{k+1}` inside the new safe set. Guaranteed by 1-Lipschitz `h_LSE`, verified at runtime.
- Target and safe set change simultaneously; `beta_h` reset (Eq. 84) keeps `h1 >= 0`.

## Input requirements (validated)
`h_eff,k(w_k) > 0`, `h_eff,k+1(w_k) > sigma`, `h_eff,1(y0) > 0`, `0 < gamma_bar <= alpha_bar`.

## Verification (`out.verification`)
`minH1`, `minBodyMargin` (physical point inside active polytope), `maxCbfViolation`, `maxClfViolation`, `maxBackstepViolation`, `nQpFallback`, `finalOutputError`, `nJumps`. Violations raise explicit `warning()` calls.

## Limitations
- The QP controller is locally Lipschitz, not C^1; the backstepping theorems (Thm. 5-6) require C^1. Finite-difference `dk/dx` is exact inside a regime and a one-sided average at regime boundaries. The smooth Gaussian-centroid controller (HYB Eq. 79) removes this gap and is not implemented.
- Compatibility (Eq. 59, `alpha_bar >= gamma_bar`) is proven for `G0 = I`; for state-dependent `G0` it is monitored at runtime.
- Unicycle body converges to a ball of radius `l` around the waypoint.
- Second-order guarantees apply to the interior of each safe set; strict-feedback depth r = 1 only.
- Safety holds only with respect to the supplied convex regions; the CBF set is an inner approximation (gap about `ln(m)/kappa`).
