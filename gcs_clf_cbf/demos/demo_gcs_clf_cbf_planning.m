%% DEMO_GCS_CLF_CBF_PLANNING  Planning stage of the GCS + CLF-CBF pipeline.
% 1. Decompose the free space of [0,10]^2 into convex regions
%    (free_space_decomposition), sort them from q0 and find adjacencies.
% 2. Build the intersection graph (build_intersection_graph): one vertex per
%    pair of adjacent regions, so the SPP in GCS returns the crossing
%    points between consecutive regions of a minimum-length path.
% 3. Remove waypoints that do not touch an obstacle (path_filter).
% Requires: YALMIP + a conic solver, Optimization Toolbox (linprog).

clear; clc;

q0 = [0.1; 0.2];
qT = [9.9; 9.9];

bbox = [0 0; 10 0; 10 10; 0 10];
obstacles = { [2 2; 4 2; 4 3; 2 3], ...                    % L-shape (part 1)
              [3 3; 3 5; 2 5; 2 3], ...                    % L-shape (part 2)
              [6 6; 8 6; 8 7; 6 7], ...                    % L-shape (part 1)
              [7 7; 7 9; 6 9; 6 7], ...                    % L-shape (part 2)
              [1 7; 2 7; 2.5 8; 2 9; 1 9; 0.5 8], ...      % convex hexagon
              [6 1; 8 1; 8 2; 6 2], ...                    % L-shape (part 1)
              [8 2; 8 3; 7 3; 7 2] };                      % L-shape (part 2)

%% Free-space decomposition and graph
regions = free_space_decomposition(bbox, obstacles);
regions = sort_regions(regions, q0);
[neighbors, A_adj] = adjacent_regions(regions);

regions_cs = cell(1, numel(regions));
for k = 1:numel(regions), regions_cs{k} = ConvexSet.polygon(regions{k}); end
g = build_intersection_graph(q0, qT, regions_cs, neighbors);

%% Shortest path
spp = ShortestPathProblem(g, 's', 't');
tic
[path, positions, cost, info] = spp.solve(false);
toc

%% Waypoint filtering
% Convex pieces of the obstacles = decomposition of the complement of the
% free regions. NOTE: the loop runs over numel(obstacles), not
% numel(obstacles_dec); the two counts can differ (open issue).
obstacles_dec = free_space_decomposition(bbox, regions);
obstacles_cs = cell(1, numel(obstacles));
for k = 1:numel(obstacles), obstacles_cs{k} = ConvexSet.polygon(obstacles_dec{k}); end
[fil_path, fil_positions, rm_i] = path_filter(obstacles_cs, path, positions);

fprintf('Path: %s\n', strjoin(path, ' -> '));
fprintf('Filtered Path: %s\n', strjoin(fil_path, ' -> '));
fprintf('Relaxation cost (lower bound): %.4f\n', info.relaxation_cost);
fprintf('Restriction cost (feasible):   %.4f\n', info.restriction_cost);
fprintf('Optimality gap:                %.4e\n', info.optimality_gap);

%% Plot
figure; hold on; axis equal;
centroids = zeros(numel(regions), 2);
for k = 1:numel(regions)
    fill(regions{k}(:,1), regions{k}(:,2), rand(1,3), 'FaceAlpha',0.6,'EdgeColor','k');
    centroids(k,:) = mean(regions{k},1);
    text(centroids(k,1), centroids(k,2), num2str(k));
end
for k = 1:numel(obstacles)
    fill(obstacles{k}(:,1), obstacles{k}(:,2), [0 0 0]);
end
for i = 1:numel(regions)
    for j = neighbors{i}
        plot(centroids([i j],1), centroids([i j],2), 'white');
    end
end

plot(q0(1), q0(2), 'ko', 'MarkerFaceColor', 'k');
plot(qT(1), qT(2), 'ko', 'MarkerFaceColor', 'k');
text(q0(1)-0.1, q0(2)+0.1, 's','FontSize',15, 'FontWeight', 'bold');
text(qT(1)-0.1, qT(2)+0.1, 't','FontSize',15, 'FontWeight','bold');

plot(positions(1,:), positions(2,:), '-o', 'LineWidth', 2, ...
    'Color', [0.85 0.1 0.1], 'MarkerFaceColor', [0.85 0.1 0.1]);
plot(positions(1,rm_i), positions(2,rm_i), 'x', 'LineWidth', 2, ...
    'Color', [0.1 0.1 0.85], 'MarkerFaceColor', [0.1 0.1 0.85]);
xlabel('x'); ylabel('y');
title(sprintf('Free space decomposed into %d convex regions \n\n Shortest path in GCS\n (path: %s,\n cost = %.3f)', ...
    numel(regions), strjoin(path, ' -> '), cost));
