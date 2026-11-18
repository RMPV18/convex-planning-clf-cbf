% FREE_SPACE_DECOMPOSITION_TEST
% Tests the decomposition of a 2D free space (bounded workspace minus non-convex obstacles)
% into a minimal set of non-overlapping convex polygons whose union
% equals the free space.

warning('on','all'); % warnings NOT suppressed

%% ---------------- Workspace and obstacle definition ----------------
bbox = [0 0; 10 0; 10 10; 0 10];

obstacles = {};
obstacles{1} = [2 2; 4 2; 4 3; 3 3; 3 5; 2 5];          % L-shape
obstacles{2} = [6 6; 8 6; 8 7; 7 7; 7 9; 6 9];          % L-shape
obstacles{3} = [1 7; 2 7; 2.5 8; 2 9; 1 9; 0.5 8];      % non-convex hexagon
obstacles{4} = [6 1; 8 1; 8 3; 7 3; 7 2; 6 2];          % L-shape


% obstacles{1} = [2 1; 2 3; 0 5; 3 4];
% obstacles{2} = [4 4.5; 4.5 5.5; 6.5 4; 4.5 3];


freeRegions = free_space_decomposition(bbox, obstacles);

%% ---------------- Verification ----------------
verifyDecomposition(freeRegions, bbox, obstacles);

%% ---------------- Plot ----------------
figure; hold on; axis equal;
for k = 1:numel(freeRegions)
    fill(freeRegions{k}(:,1), freeRegions{k}(:,2), rand(1,3), 'FaceAlpha',0.6,'EdgeColor','k');
end
for k = 1:numel(obstacles)
    fill(obstacles{k}(:,1), obstacles{k}(:,2), [0 0 0]);
end
title(sprintf('Free space decomposed into %d convex regions', numel(freeRegions)));


function verifyDecomposition(freeRegions, bbox, obstacles)
% Checks: no pairwise overlap (positive area intersection),
% union of regions matches free space area within tolerance.

nR = numel(freeRegions);
polys = cell(nR,1);
for k = 1:nR
    polys{k} = polyshape(freeRegions{k}(:,1), freeRegions{k}(:,2));
end

for i = 1:nR
    for j = i+1:nR
        inter = intersect(polys{i}, polys{j});
        if area(inter) > 1e-8
            warning('Overlap detected between region %d and %d (area=%g).', i, j, area(inter));
        end
    end
end

pbox = polyshape(bbox(:,1), bbox(:,2));
pobs = polyshape();
for k = 1:numel(obstacles)
    pk = polyshape(obstacles{k}(:,1), obstacles{k}(:,2));
    pobs = union(pobs, pk);
end
pfreeTrue = subtract(pbox, pobs);

pUnion = polyshape();
for k = 1:nR
    pUnion = union(pUnion, polys{k});
end

diffArea = area(subtract(pfreeTrue, pUnion)) + area(subtract(pUnion, pfreeTrue));
if diffArea > 1e-6
    warning('Union of decomposed regions does not exactly match free space (area mismatch=%g).', diffArea);
end

fprintf('Regions: %d | Free space area: %.4f | Covered area: %.4f\n', ...
    nR, area(pfreeTrue), area(pUnion));
end