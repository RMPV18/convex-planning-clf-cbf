%% demo_2d_example.m
% Recreates, in spirit, the repo's 2d_example: a source point and a
% target point separated by an "obstacle", with two alternative convex
% regions (above/below) connecting them through an intermediate region.
% The solver picks not just which region to route through, but the
% optimal point within each region, to minimize total Euclidean length.
%
% Requires YALMIP + a conic solver (Mosek, Gurobi, SCS, SDPT3, ECOS, ...).

clear; clc;

% Source and target are singleton point sets.
s = [0.1;0.2];
t = [9.9;9.9];


bbox = [0 0; 10 0; 10 10; 0 10];

obstacles = {};
obstacles{1} = [2 2; 4 2; 4 3; 2 3];          % L-shape(1)
obstacles{2} = [3 3; 3 5; 2 5; 2 3];          % L-shape(2)
obstacles{3} = [6 6; 8 6; 8 7; 6 7];          % L-shape(1)
obstacles{4} = [7 7; 7 9; 6 9; 6 7];          % L-shape(2)
obstacles{5} = [1 7; 2 7; 2.5 8; 2 9; 1 9; 0.5 8];      % convex hexagon
obstacles{6} = [6 1; 8 1; 8 2; 6 2];          % L-shape(1)
obstacles{7} = [8 2; 8 3; 7 3; 7 2];          % L-shape(2)


freeRegions = free_space_decomposition(bbox, obstacles);

freeRegions = sortRegions(freeRegions, s);

[nb, Aadj]  = adjacentRegions(freeRegions); 

for k = 1:numel(freeRegions), freeRegions_cs{k} = ConvexSet.polygon(freeRegions{k}); end
g = GraphOfConvexSets.add_mygraph(s, t, freeRegions_cs, nb);

spp = ShortestPathProblem(g, 's', 't'); 
tic
[path, positions, cost, info] = spp.solve(false);
toc
obstacles_dec = free_space_decomposition(bbox, freeRegions);
for k = 1:numel(obstacles), obstacles_cs{k} = ConvexSet.polygon(obstacles_dec{k}); end
[fil_path, fil_positions, rm_i] = pathFilter(obstacles_cs, path, positions);

fprintf('Path: %s\n', strjoin(path, ' -> '));
fprintf('Filtered Path: %s\n', strjoin(fil_path, ' -> '));
fprintf('Relaxation cost (lower bound): %.4f\n', info.relaxation_cost);
fprintf('Restriction cost (feasible):   %.4f\n', info.restriction_cost);
fprintf('Optimality gap:                %.4e\n', info.optimality_gap);



%% ---------------- Plot ----------------
figure; hold on; axis equal;
for k = 1:numel(freeRegions)
    fill(freeRegions{k}(:,1), freeRegions{k}(:,2), rand(1,3), 'FaceAlpha',0.6,'EdgeColor','k');
    mfrv{k} = mean(freeRegions{k},1);
    text(mfrv{k}(1), mfrv{k}(2), num2str(k));
end
for k = 1:numel(obstacles)
    fill(obstacles{k}(:,1), obstacles{k}(:,2), [0 0 0]);
end
for i = 1:numel(freeRegions)
    for j = nb{i}
        plot([mfrv{i}(1), mfrv{j}(1)] , [mfrv{i}(2), mfrv{j}(2)], 'white');
    end
end
% legend(adj_net, 'net')

plot(s(1), s(2), 'ko', 'MarkerFaceColor', 'k');
plot(t(1), t(2), 'ko', 'MarkerFaceColor', 'k');
text(s(1)-0.1, s(2)+0.1, 's','FontSize',15, 'FontWeight', 'bold');
text(t(1)-0.1, t(2)+0.1, 't','FontSize',15, 'FontWeight','bold');

plot(positions(1,:), positions(2,:), '-o', 'LineWidth', 2, ...
    'Color', [0.85 0.1 0.1], 'MarkerFaceColor', [0.85 0.1 0.1]);
plot(positions(1,rm_i), positions(2,rm_i), 'x', 'LineWidth', 2, ...
    'Color', [0.1 0.1 0.85], 'MarkerFaceColor', [0.1 0.1 0.85]);
% title(sprintf('Shortest path in GCS\n (path: %s,\n cost = %.3f)', ...
%     strjoin(path, ' -> '), cost));
xlabel('x'); ylabel('y');


title(sprintf('Free space decomposed into %d convex regions \n\n Shortest path in GCS\n (path: %s,\n cost = %.3f)', ...
    numel(freeRegions), strjoin(path, ' -> '), cost));


% obstacles = {};
% obstacles{1} = [2 2; 4 2; 4 3; 3 3; 3 5; 2 5];          % L-shape
% obstacles{2} = [6 6; 8 6; 8 7; 7 7; 7 9; 6 9];          % L-shape
% obstacles{3} = [1 7; 2 7; 2.5 8; 2 9; 1 9; 0.5 8];      % non-convex hexagon
% obstacles{4} = [6 1; 8 1; 8 3; 7 3; 7 2; 6 2];          % L-shape

% freeRegions = free_space_decomposition(bbox, obstacles);

% obstacles = {};
% obstacles = free_space_decomposition(bbox, freeRegions);

