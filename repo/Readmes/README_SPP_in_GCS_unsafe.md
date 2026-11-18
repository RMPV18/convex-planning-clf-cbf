# Shortest Paths in Graphs of Convex Sets — MATLAB

A MATLAB adaptation of the algorithm from:

> T. Marcucci, J. Umenberger, P. A. Parrilo, R. Tedrake, "Shortest Paths in
> Graphs of Convex Sets," SIAM J. Optim. 34(1), 2024.
> Python repo: https://github.com/TobiaMarcucci/shortest-paths-in-graphs-of-convex-sets

**Note on provenance:** GitHub blocks automated fetching of that repo's raw
file tree from this environment, so this is not a line-by-line port of the
Python source. It reimplements the same algorithmic pipeline described in
the paper (and mirrored in the repo's structure: `spp.graph`,
`spp.convex_sets`, `spp.shortest_path`) directly in MATLAB:

1. **Convex relaxation** of the perspective mixed-integer convex program
   (MICP): edge-selection variables `y_e` relaxed to `[0,1]`, with each
   edge carrying homogenized ("perspective") copies of its endpoints so
   the vertex-set and edge-length constraints stay convex.
2. **Rounding**: a `y`-weighted depth-first search recovers a single
   discrete source-to-target path from the (possibly fractional)
   relaxed solution.
3. **Convex restriction**: with the path fixed, a small convex program
   (no binaries) optimizes the continuous position within each visited
   set, giving a feasible upper bound and an optimality gap against the
   relaxation's lower bound.

Only polyhedral convex sets (`{x : A x <= b}` — boxes, polygons, points,
general halfspace intersections) are implemented, and edge cost is
Euclidean distance. Both are the common case in the repo's motion-planning
style examples; extending to other convex sets or edge costs means adding
a class with the same `constraints(x)` / `perspective_constraints(x,y)`
interface as `ConvexSet.m`.

## Files

- `ConvexSet.m` — polyhedral set + its perspective (homogenized) form.
- `GraphOfConvexSets.m` — named vertices (each a `ConvexSet`) + directed edges.
- `ShortestPathProblem.m` — relaxation, rounding, restriction.
- `demo_2d_example.m` — obstacle-avoidance-style 2D example with a plot.
- `demo_sanity_check.m` — shrinks all sets to points and cross-checks the
  result against MATLAB's built-in `shortestpath` (classic Dijkstra).

## Requirements

- [YALMIP](https://yalmip.github.io/) on the MATLAB path.
- Any conic solver YALMIP can call: Mosek, Gurobi, SCS, SDPT3, ECOS, etc.
  (Mosek is what the original repo uses; any of the others works equally
  well here since the problems are ordinary SOCPs.)

## Usage

```matlab
g = GraphOfConvexSets();
g.add_vertex('s', ConvexSet.point([0;0]));
g.add_vertex('A', ConvexSet.box([1;1], [3;3]));
g.add_vertex('t', ConvexSet.point([5;0]));
g.add_edge('s','A');
g.add_edge('A','t');

spp = ShortestPathProblem(g, 's', 't');
[path, positions, cost, info] = spp.solve();
```

Run `demo_2d_example.m` for a plotted example, and `demo_sanity_check.m`
to verify the point-set special case matches classic Dijkstra.
