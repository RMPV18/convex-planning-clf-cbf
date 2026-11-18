# Waypoint-following closed-form CLF-CBF QP (unicycle / double integrator)

Files (each controller is self-contained: parameters, region preparation, region CBF, segment geometry and the closed-form QP are local functions at the end of the file)
- `unicycle_waypoint_clf_cbf.m` : first-order version (offset point, rho = l).
- `double_integrator_waypoint_clf_cbf.m` : second-order version (joint CLF-CBF backstepping).
- `run_waypoint_simulation.m` : script; edit `MODEL` and `SCENARIO` at the top. `CHECK_QP = true` compares the closed-form QP with `quadprog` (Optimization Toolbox) and checks that the QP copies in the controller files match the script's copy.

Requirements: MATLAB R2018b+ (`yline`, `tiledlayout`); base MATLAB only.

## Calls
- `unicycle_waypoint_clf_cbf(W, regions, [px; py; th], l, params)`
- `double_integrator_waypoint_clf_cbf(W, regions, x0, z0, params)`
- `regions(k)`: convex polytope `A` (m x n), `b` (m x 1); union of 2: `A = {A1, A2}`, `b = {b1, b2}`.

## Method
Mode k = (CLF target `w_k`, safe set `Q_k`). `V = 1/2|y - w_k|^2`, `h = h_region - rho`.
Closed-form QP (Matias & Silvestre Eq. 10-15, p = inf): `min 1/2|v|^2  s.t.  LGV v <= -2 gamma_bar V,  LGh v >= -alpha_bar h`.
- Unicycle: `u = v` with `G = R(th) diag(1,l)`, `rho = l`.
- Double integrator: `v = k(x)` is the virtual control; `e = z - k`, `V1`, `h1`; `u = (B/|e|^2) e` for `B = min(B_V, B_h) < 0`, else 0; `dk/dx` by central differences; `beta_h` reset at jumps (Eq. 84).
Switching at `|y - w_k| <= r_k`; target and safe set change together. Sight is guaranteed by construction (1-Lipschitz CBF; for a union region, tube test on segment `w_k -> w_{k+1}`) and verified at each jump.

## Limitations
- QP controller is Lipschitz, not C^1: backstepping theorems (Thm. 5-6) hold only approximately at regime boundaries.
- Union regions are non-convex: compatibility and deadlock-freedom are not guaranteed; stalls raise `wp:timeout`.
- Only unions of 2 polytopes; second-order depth r = 1.
- Safety holds only with respect to the supplied regions; CBF is an inner approximation (gap about `ln(m)/kappa`).
