%% DEMO_GCS_CLF_CBF_FOLLOWING  Full GCS + CLF-CBF pipeline (plan, then follow).
% Planning as in demo_gcs_clf_cbf_planning (decomposition -> intersection
% graph -> SPP in GCS), then the crossing points are followed by the
% double-integrator CLF-CBF waypoint controller. The safe set while
% heading to waypoint k is built by regions_selection (below).
% Uses the work-in-progress tracker wip/double_integrator_my_waypoint_clf_cbf.
% Requires: YALMIP + a conic solver, Optimization Toolbox (linprog).

clear; clc;

%% Planning
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

regions = free_space_decomposition(bbox, obstacles);
regions = sort_regions(regions, q0);
[neighbors, A_adj] = adjacent_regions(regions);

regions_cs = cell(1, numel(regions));
for k = 1:numel(regions), regions_cs{k} = ConvexSet.polygon(regions{k}); end
g = build_intersection_graph(q0, qT, regions_cs, neighbors);

spp = ShortestPathProblem(g, 's', 't');
[path, positions, cost, info] = spp.solve(false);

path = clear_path(path);
fprintf('---------Planning----------\n');
fprintf('Path: %s\n', strjoin(path, ' -> '));
fprintf('Relaxation cost (lower bound): %.4f\n', info.relaxation_cost);
fprintf('Restriction cost (feasible):   %.4f\n', info.restriction_cost);
fprintf('Optimality gap:                %.4e\n', info.optimality_gap);

%% Following
safe_sets = regions_selection(path, regions_cs);
W = positions(:, 2:end);     % waypoints: crossing points and qT
V0 = [0; 0];                 % initial velocity
params = struct();           % e.g. struct('kappa', 30, 'alpha_bar', 3)

out = double_integrator_my_waypoint_clf_cbf(W, safe_sets, q0, V0, params);

fprintf('\n---------Tracking----------\njumps [t from to (beta_h)]:\n');
disp(out.jumps);
fn = fieldnames(out.verification);
for k = 1:numel(fn)
    fprintf('%-22s = %g\n', fn{k}, out.verification.(fn{k}));
end
fprintf('expected jumps         = %d\n', size(W, 2) - 1);

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
xlabel('x'); ylabel('y');
title(sprintf('Free space decomposed into %d convex regions \n\n Shortest path in GCS\n (path: %s,\n cost = %.3f)', ...
    numel(regions), strjoin(path, ' -> '), cost));


function region_path = clear_path(path)
% Converts the intersection-graph path {'s','i-j','j-k',...,'t'} into the
% sequence of traversed regions {'i','j','k',...} (needs >= 2 pair vertices).
pairs = path(2:end-1);
e = cell2mat(cellfun(@(s) sscanf(s, '%d-%d')', pairs(:), 'UniformOutput', false));
m = size(e, 1);
assert(m >= 2, 'need at least 2 edges between path{1} and path{end}');
nodes = zeros(1, m+1);
nodes(1) = setdiff(e(1,:), e(2,:));
for k = 1:m-1
    nodes(k+1) = intersect(e(k,:), e(k+1,:));
end
nodes(m+1) = setdiff(e(m,:), e(m-1,:));
region_path = cellstr(string(nodes));
end


function safe_sets = regions_selection(region_path, regions_cs)
% safe_sets(k), k < M: halfspaces of regions k and k+1 of the path stacked,
% with the rows of their shared facet removed. safe_sets(M): last region.
M = numel(region_path);
safe_sets = struct('A', cell(1, M), 'b', cell(1, M));
for k = 1:M-1
    Ri = regions_cs{str2double(region_path{k})};
    Rj = regions_cs{str2double(region_path{k+1})};
    Ao = [Ri.A; Rj.A];  bo = [Ri.b; Rj.b];
    [fi, ~] = find(abs(Ri.A*Rj.A' + 1) < 1e-6 & ...
                   abs(Ri.b + Rj.b')  < 1e-6*max(1,norm(bo,inf)));
    nrm = Ri.A(fi,:)';
    pn  = vecnorm(Ao*(eye(size(nrm,1)) - nrm*nrm'), 2, 2);   % in-plane normal lengths
    keep = pn > 1e-9;
    safe_sets(k).A = Ao(keep,:);
    safe_sets(k).b = bo(keep,:);
end
safe_sets(M).A = regions_cs{str2double(region_path{M})}.A;
safe_sets(M).b = regions_cs{str2double(region_path{M})}.b;
end
