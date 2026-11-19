# convex-planning-clf-cbf

MATLAB implementations of convex-optimization-based path planners and the
trackers that follow their output. The code is organised so that
planner + tracker **pipelines** can be compared side by side.

## Layout

| Folder | Content |
|---|---|
| `common/` | Shared by all pipelines: `ConvexSet`, `GraphOfConvexSets`, `ShortestPathProblem` (SPP in GCS), `free_space_decomposition`, `adjacent_regions`, `sort_regions`. `inflate_obstacles`, `convex_obstacle_pieces`, `firi_polytope`, `path_clearance`. Demos in `common/*/demo*/`. |
| `gcs_clf_cbf/` | Free-space decomposition → SPP in GCS on the intersection graph (`build_intersection_graph`, `path_filter`) → closed-form CLF-CBF waypoint following (`double_integrator_waypoint_clf_cbf`, `unicycle_waypoint_clf_cbf`); safe-set chain tracker (`build_safe_chain`, `translate_waypoint`, `chain_targets`, `double_integrator_safe_chain_clf_cbf`) with secondary versions in `online_firi/` and `gcs_bezier/`, compared by `demos/compare_safe_set_versions`. `wip/` holds work in progress, `tocheck/` review copies (both git-ignored or off the MATLAB path). |
| `gcs_pid/` | `GcsMotionPlanner` (Bézier trajectories in GCS, Marcucci et al. 2023) → joint PID (`track_gcs_pid`). |
| `dcp_pid/` | `dcp_polygon_planner` (penalty convex-concave procedure, CVX) → saturated discrete PID (`pid_trajectory_tracker`). |
| `library/` | Reusable code not called by the current pipelines: `line_of_sight/`, `clf_cbf_qp_closed_form` (relaxed-CLF variant), `lse_cbf`, `vert2hrep`. |

Each pipeline folder has a `demos/` subfolder. `demo_<pipeline>_planning`
runs the planner alone, `demo_<pipeline>_plan_track` / `_following` runs
planning and tracking together.

## Pipelines

**SPP in GCS** (`common/ShortestPathProblem.m`, Marcucci, Umenberger,
Parrilo, Tedrake, *Shortest Paths in Graphs of Convex Sets*, SIAM J. Optim.
2024, arXiv:2101.11565): convex relaxation of the perspective MICP →
y-weighted DFS rounding → convex restriction on the fixed path (feasible
upper bound; optimality gap vs. the relaxation). Polyhedral sets and
Euclidean edge length only. Not a line-by-line port of the Python repo;
reimplemented from the paper.

**GCS + CLF-CBF** (`gcs_clf_cbf/`): mode k = (CLF target `w_k`, safe set
`Q_k`), `V = 1/2|y - w_k|^2`, `h = h_region - rho`; closed-form QP
(Matias & Silvestre, Eq. 10-15, p = inf)
`min 1/2|v|^2 s.t. LGV v <= -2 gamma_bar V, LGh v >= -alpha_bar h`.
Unicycle: offset point, `rho = l`. Double integrator: CLF-CBF backstepping,
`beta_h` reset at jumps (Eq. 84). Switching at `|y - w_k| <= r_k`. Each
controller file is self-contained (all helpers are local functions);
`demo_clf_cbf_waypoint_simulation` with `CHECK_QP = true` checks the
closed-form QP against `quadprog` and that the local QP copies are identical.
Limitations: the QP controller is Lipschitz, not C^1; union regions (max 2
polytopes) are non-convex, so deadlock-freedom is not guaranteed
(`wp:timeout`); the CBF is an inner approximation (gap ≈ `ln(m)/kappa`).

**Safe-set chain** (`gcs_clf_cbf/`, `demo_gcs_clf_cbf_following`): offline,
once, after planning (the planner does not use these sets): for every facet
between consecutive path regions a FIRI polytope (Wang et al., T-RO 2025;
`common/regions/firi_polytope`) seeded on the facet and confined to the
bounding box of the two regions (Gonzalez-Garcia et al. heuristic) bridges
them, giving the chain region, bridge, region, ... Each crossing point is
translated along `n_qhat` into the intersection of the adjacent sets with
margin `mu` (Matias & Silvestre Eqs. 32, 42-44 adapted). Tracking: LSE CBF of
the active set, `alpha_bar >= gamma_bar` (no CLF/CBF conflict inside a convex
set), joint backstepping, after-next target when in sight, jump when the
state enters the next set (`h_{i+1} >= sigma`). Secondary versions:
`online_firi/` (FIRI set computed online from the current state and the
target) and `gcs_bezier/` (waypoints where a `GcsMotionPlanner` trajectory
crosses the `q_hat` hyperplanes).

**GCS + PID** (`gcs_pid/`): `GcsMotionPlanner` (Marcucci, Petersen, von
Wrangel, Tedrake, *Motion Planning around Obstacles with Convex
Optimization*, Sci. Robot. 8, eadf7843, 2023) extends `ShortestPathProblem`
with Bézier control points, time scaling, duration/length/energy costs,
edge pruning and randomized rounding; `track_gcs_pid` tracks the sampled
trajectory on a double integrator.

**DCP + PID** (`dcp_pid/`): `dcp_polygon_planner` assigns one obstacle face
per segment (penalty CCP, several homotopy starts); `pid_trajectory_tracker`
tracks it on a damped point mass with force saturation and anti-windup.

## Conventions

- Classes `UpperCamelCase`, functions and scripts `snake_case`.
- Loop indices: `i`, `j` pairs of regions/vertices; `k` sequences
  (regions in a list, waypoints, segments, time steps); `e` graph edges;
  `f` facets; `l` derivative order.
- `q0`, `qT` start/goal points; `regions` vertex arrays, `regions_cs`
  `ConvexSet` objects; `neighbors` adjacency lists; `obstacles` obstacle
  vertex arrays; `params` user option struct (missing fields → defaults).

## Dependencies (not included)

- [YALMIP](https://yalmip.github.io) and a conic solver (Mosek recommended;
  SeDuMi/SDPT3/SCS/ECOS work) — `common/`, `gcs_clf_cbf/` planning, `gcs_pid/`.
- [CVX](http://cvxr.com/cvx/) — `dcp_pid/`.
- MATLAB Optimization Toolbox (`linprog`, `quadprog`).

`startup.m` runs `../startup_local.m` (local, not in the repo), then adds
YALMIP and MOSEK. MOSEK is appended at the end of the path because its
toolbox ships its own `linprog`/`quadprog`. Run `cvx_setup` once from the
CVX folder; it saves the CVX path.
