%% DEMO_FREE_SPACE_DECOMPOSITION  Convex decomposition of a 2-D workspace.
% Decomposes [0,10]^2 minus four non-convex obstacles with
% free_space_decomposition, verifies that the regions do not overlap and
% that their union equals the free space, and plots the result.
% Requires: MATLAB (polyshape).

clear; clc;

bbox = [0 0; 10 0; 10 10; 0 10];
obstacles = { [2 2; 4 2; 4 3; 3 3; 3 5; 2 5], ...          % L-shape
              [6 6; 8 6; 8 7; 7 7; 7 9; 6 9], ...          % L-shape
              [1 7; 2 7; 2.5 8; 2 9; 1 9; 0.5 8], ...      % hexagon
              [6 1; 8 1; 8 3; 7 3; 7 2; 6 2] };            % L-shape

regions = free_space_decomposition(bbox, obstacles);
verify_decomposition(regions, bbox, obstacles);

%% Plot
figure; hold on; axis equal;
for k = 1:numel(regions)
    fill(regions{k}(:,1), regions{k}(:,2), rand(1,3), 'FaceAlpha',0.6,'EdgeColor','k');
end
for k = 1:numel(obstacles)
    fill(obstacles{k}(:,1), obstacles{k}(:,2), [0 0 0]);
end
title(sprintf('Free space decomposed into %d convex regions', numel(regions)));


function verify_decomposition(regions, bbox, obstacles)
% Warns on any pairwise overlap of positive area and on any area mismatch
% between the union of the regions and the true free space.
K = numel(regions);
polys = cell(K,1);
for k = 1:K
    polys{k} = polyshape(regions{k}(:,1), regions{k}(:,2));
end

for i = 1:K
    for j = i+1:K
        inter = intersect(polys{i}, polys{j});
        if area(inter) > 1e-8
            warning('demo_free_space_decomposition:overlap', ...
                'Overlap detected between region %d and %d (area=%g).', i, j, area(inter));
        end
    end
end

pobs = polyshape();
for k = 1:numel(obstacles)
    pobs = union(pobs, polyshape(obstacles{k}(:,1), obstacles{k}(:,2)));
end
pfree = subtract(polyshape(bbox(:,1), bbox(:,2)), pobs);

pUnion = polyshape();
for k = 1:K
    pUnion = union(pUnion, polys{k});
end

diffArea = area(subtract(pfree, pUnion)) + area(subtract(pUnion, pfree));
if diffArea > 1e-6
    warning('demo_free_space_decomposition:coverage', ...
        'Union of decomposed regions does not match free space (area mismatch=%g).', diffArea);
end
fprintf('Regions: %d | Free space area: %.4f | Covered area: %.4f\n', ...
    K, area(pfree), area(pUnion));
end
