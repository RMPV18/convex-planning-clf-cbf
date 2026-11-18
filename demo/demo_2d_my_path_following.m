%% demo_2d_followingPath_example.m
% Recreates, in spirit, the repo's 2d_example: a source point and a
% target point separated by an "obstacle", with two alternative convex
% regions (above/below) connecting them through an intermediate region.
% The solver picks not just which region to route through, but the
% optimal point within each region, to minimize total Euclidean length.
%
% Requires YALMIP + a conic solver (Mosek, Gurobi, SCS, SDPT3, ECOS, ...).

clear; clc;

%%---------------- Path Planning ----------------

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
[path, positions, cost, info] = spp.solve(false);

path = clear_path(path);
fprintf('---------Planning----------\n');
fprintf('Path: %s\n', strjoin(path, ' -> '));
fprintf('Relaxation cost (lower bound): %.4f\n', info.relaxation_cost);
fprintf('Restriction cost (feasible):   %.4f\n', info.restriction_cost);
fprintf('Optimality gap:                %.4e\n', info.optimality_gap);

%%---------------- Path Following ----------------

REGIONS = regions_selection(path, freeRegions_cs);  % Union
W = positions(:, 2:end);     % Previously obtained
X0 = s;                      % initial position
V0 = [0; 0];                 % initial velocity, 'double_integrator' only
PARAMS = struct();           % e.g. struct('kappa', 30, 'alpha_bar', 3)

out = double_integrator_my_waypoint_clf_cbf(W, REGIONS, X0, V0, PARAMS);

fprintf('\n---------Tracking----------\njumps [t from to (beta_h)]:\n');
disp(out.jumps);
fn = fieldnames(out.verification);
for i = 1:numel(fn)
    fprintf('%-22s = %g\n', fn{i}, out.verification.(fn{i}));
end
fprintf('expected jumps         = %d\n', size(W, 2) - 1);

% if PLOT
%     plot_results(out, boxes, W, X0, MODEL, SCENARIO);
% end



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




function clean_path = clear_path(path)
    % path{1} and path{end} are ignored, matching the original indexing
    p = path(2:end-1);
    e = cell2mat(cellfun(@(s) sscanf(s, '%d-%d')', p(:), 'UniformOutput', false));
    m = size(e, 1);
    assert(m >= 2, 'need at least 2 edges between path{1} and path{end}');
    nodes = zeros(1, m+1);
    nodes(1) = setdiff(e(1,:), e(2,:));
    for k = 1:m-1
        nodes(k+1) = intersect(e(k,:), e(k+1,:));
    end
    nodes(m+1) = setdiff(e(m,:), e(m-1,:));
    clean_path = cellstr([string(nodes)]);
end


function regions = regions_selection(path, freeR);
    M = numel(path);
    regions = struct('A', cell(1, M), 'b', cell(1, M));
    for i = 1:numel(path)-1
        A{i}  = freeR{str2num(path{i})}.A;
        A{i+1}  = freeR{str2num(path{i+1})}.A;
        b{i}  = freeR{str2num(path{i})}.b;
        b{i+1}  = freeR{str2num(path{i+1})}.b;
         Ao = [A{i}; A{i+1}];  bo = [b{i}; b{i+1}];
        [fr, fc] = find(abs(A{i}*A{i+1}' + 1) < 1e-6 & ...
                        abs(b{i} + b{i+1}')   < 1e-6*max(1,norm(bo,inf)));
        nrm = A{i}(fr,:)';
        pn  = vecnorm(Ao*(eye(size(nrm,1)) - nrm*nrm'), 2, 2);   % in-plane normal lengths
        keep = pn > 1e-9;
        regions(i).A = Ao(keep,:);
        regions(i).b = bo(keep,:);
    end
    regions(i+1).A = freeR{str2num(path{i+1})}.A;
    regions(i+1).b = freeR{str2num(path{i+1})}.b;
end