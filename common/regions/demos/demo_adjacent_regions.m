%% DEMO_ADJACENT_REGIONS  Adjacency graph of a convex free-space decomposition.
% Decomposes the workspace of demo_free_space_decomposition, computes the
% facet-adjacency of the regions with adjacent_regions and draws the
% adjacency graph between region centroids.
% Requires: Optimization Toolbox (linprog).

clear; clc;

bbox = [0 0; 10 0; 10 10; 0 10];
obstacles = { [2 2; 4 2; 4 3; 3 3; 3 5; 2 5], ...          % L-shape
              [6 6; 8 6; 8 7; 7 7; 7 9; 6 9], ...          % L-shape
              [1 7; 2 7; 2.5 8; 2 9; 1 9; 0.5 8], ...      % hexagon
              [6 1; 8 1; 8 3; 7 3; 7 2; 6 2] };            % L-shape

regions = free_space_decomposition(bbox, obstacles);
[neighbors, A_adj] = adjacent_regions(regions);

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
title('Adjacent Free Regions');
