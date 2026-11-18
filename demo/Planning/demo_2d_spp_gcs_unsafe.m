%% demo_2d_example.m
% Recreates, in spirit, the repo's 2d_example: a source point and a
% target point separated by an "obstacle", with two alternative convex
% regions (above/below) connecting them through an intermediate region.
% The solver picks not just which region to route through, but the
% optimal point within each region, to minimize total Euclidean length.
%
% Requires YALMIP + a conic solver (Mosek, Gurobi, SCS, SDPT3, ECOS, ...).

clear; clc;

g = GraphOfConvexSets();

% Source and target are singleton point sets.
g.add_vertex('s', ConvexSet.point([0; 0]));
g.add_vertex('t', ConvexSet.point([9; 0]));

% Two alternative "corridor" regions around an obstacle centered near x=4.
g.add_vertex('A', ConvexSet.box([2; 1],  [5; 4]));   % region above the obstacle
g.add_vertex('B', ConvexSet.box([2; -4], [5; -1]));  % region below the obstacle

% A convergence region both paths must pass through before the target.
g.add_vertex('C', ConvexSet.box([6; -2], [8; 2]));

g.add_edge('s', 'A');
g.add_edge('s', 'B');
g.add_edge('A', 'C');
g.add_edge('B', 'C');
g.add_edge('C', 't');

spp = ShortestPathProblem(g, 's', 't');
tic
[path, positions, cost, info] = spp.solve(false);
toc

fprintf('Path: %s\n', strjoin(path, ' -> '));
fprintf('Relaxation cost (lower bound): %.4f\n', info.relaxation_cost);
fprintf('Restriction cost (feasible):   %.4f\n', info.restriction_cost);
fprintf('Optimality gap:                %.4e\n', info.optimality_gap);

%% Plot
figure; hold on; axis equal; grid on;
colors = struct('s', [0 0 0], 't', [0 0 0], 'A', [0.2 0.6 1], ...
                 'B', [1 0.6 0.2], 'C', [0.6 0.8 0.4]);

for i = 1:g.num_vertices()
    name = g.vertex_names{i};
    cs = g.vertices(name);
    if isequal(cs.A, [eye(2); -eye(2)]) && all(cs.b(1:2) ~= -cs.b(3:4))
        % box: plot as a filled rectangle
        lb = [-cs.b(3); -cs.b(4)];
        ub = [cs.b(1); cs.b(2)];
        rectangle('Position', [lb(1), lb(2), ub(1)-lb(1), ub(2)-lb(2)], ...
            'FaceColor', [colors.(name), 0.3], 'EdgeColor', colors.(name));
    end
end

plot(0, 0, 'ko', 'MarkerFaceColor', 'k');
plot(9, 0, 'ko', 'MarkerFaceColor', 'k');
text(0, -0.5, 's');
text(9, -0.5, 't');

plot(positions(1,:), positions(2,:), '-o', 'LineWidth', 2, ...
    'Color', [0.85 0.1 0.1], 'MarkerFaceColor', [0.85 0.1 0.1]);
title(sprintf('Shortest path in GCS (path: %s, cost = %.3f)', ...
    strjoin(path, ' -> '), cost));
xlabel('x'); ylabel('y');
