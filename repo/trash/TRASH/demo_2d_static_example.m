%% demo_2d_static_example.m
% Recreates, in spirit, the repo's 2d_example: a source point and a
% target point separated by an "obstacle", with two alternative convex
% regions (above/below) connecting them through an intermediate region.
% The solver picks not just which region to route through, but the
% optimal point within each region, to minimize total Euclidean length.
%
% Requires YALMIP + a conic solver (Mosek, Gurobi, SCS, SDPT3, ECOS, ...).

clear; clc;

g = GraphOfConvexSets();

obstacles{1} = [2 2; 4 2; 4 3; 3 3; 3 5; 2 5];          % L-shape


% Source and target are singleton point sets.
g.add_vertex('s', ConvexSet.point([0,0]));
g.add_vertex('t', ConvexSet.point([10,10]));

% Obstacle vertices are singleton point sets. 
%   - Discover the respective obstacle convex hull vertices and edges
 
for i = 1:length(obstacles)
    g.add_vertex(['o',num2str(i), '-', 'v1'], ConvexSet.point([obstacles{i}(1,1); obstacles{i}(1,2)]));
    for j = 2:length(obstacles{i})
        g.add_vertex(['o',num2str(i), '-', 'v', num2str(j)], ConvexSet.point([obstacles{i}(j,1); obstacles{i}(j,2)]));
        g.add_edge(['o',num2str(i), '-', 'v', num2str(j)], ['o',num2str(i), '-', 'v', num2str(j-1)]);
        g.add_edge(['o',num2str(i), '-', 'v', num2str(j-1)], ['o',num2str(i), '-', 'v', num2str(j)]);
    end
    g.add_edge(['o',num2str(i), '-', 'v', num2str(j)], ['o',num2str(i), '-', 'v', num2str(1)]);
    g.add_edge(['o',num2str(i), '-', 'v', num2str(1)], ['o',num2str(i), '-', 'v', num2str(j)]);
end



% g.add_edge('A', 'B');
% g.add_edge('A', 'F');
% g.add_edge('B', 'E');
% 
% g.add_edge('F', 'C');
% 
% g.add_edge('E', 'D');
% 
% g.add_edge('C', 'D');
% 
% g.add_edge('s', 'A');
% g.add_edge('s', 'B');
% g.add_edge('s', 'F');
% 
% g.add_edge('B', 't');
% g.add_edge('F', 't');
% g.add_edge('E', 't');
% g.add_edge('C', 't');
% g.add_edge('D', 't');

spp = ShortestStaticPathProblem(g, 's', 't');
tic
[path, positions, cost, info] = spp.solve(false);
toc

fprintf('Path: %s\n', strjoin(path, ' -> '));
fprintf('Relaxation cost (lower bound): %.4f\n', info.relaxation_cost);
fprintf('Restriction cost (feasible):   %.4f\n', info.restriction_cost);
fprintf('Optimality gap:                %.4e\n', info.optimality_gap);

%% Plot
figure; hold on; axis equal; grid on;

for k = 1:numel(obstacles)
    fill(obstacles{k}(:,1), obstacles{k}(:,2), [0 0 0]);
end

plot(0, 0, 'ko', 'MarkerFaceColor', 'k');
plot(10, 10, 'ko', 'MarkerFaceColor', 'k');
text(0, -0.5, 's');
text(10, 9.5, 't');

plot(positions(1,:), positions(2,:), '-o', 'LineWidth', 2, ...
    'Color', [0.85 0.1 0.1], 'MarkerFaceColor', [0.85 0.1 0.1]);
title(sprintf('Shortest path in GCS (path: %s, cost = %.3f)', ...
    strjoin(path, ' -> '), cost));
xlabel('x'); ylabel('y');
