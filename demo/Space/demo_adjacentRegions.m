clear; clc;

bbox = [0 0; 10 0; 10 10; 0 10];

obstacles = {};
obstacles{1} = [2 2; 4 2; 4 3; 3 3; 3 5; 2 5];          % L-shape
obstacles{2} = [6 6; 8 6; 8 7; 7 7; 7 9; 6 9];          % L-shape
obstacles{3} = [1 7; 2 7; 2.5 8; 2 9; 1 9; 0.5 8];      % non-convex hexagon
obstacles{4} = [6 1; 8 1; 8 3; 7 3; 7 2; 6 2];          % L-shape

freeRegions = free_space_decomposition(bbox, obstacles);

[nb, Aadj]   = adjacentRegions(freeRegions);  


%% ---------------- Plot ----------------
figure; hold on; axis equal;
for k = 1:numel(freeRegions)
    fill(freeRegions{k}(:,1), freeRegions{k}(:,2), rand(1,3), 'FaceAlpha',0.6,'EdgeColor','k');
    mfrv(k,:) = mean(freeRegions{k},1);
    text(mfrv(k,1), mfrv(k,2), num2str(k));
end
for k = 1:numel(obstacles)
    fill(obstacles{k}(:,1), obstacles{k}(:,2), [0 0 0]);
end
for i = 1:numel(freeRegions)
    for j = nb{i}
        plot([mfrv(i,1), mfrv(j,1)] , [mfrv(i,2), mfrv(j,2)], 'white');
    end
end

title('Adjacent Free Regions');
