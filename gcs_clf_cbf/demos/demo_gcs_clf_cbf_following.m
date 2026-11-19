%% DEMO_GCS_CLF_CBF_FOLLOWING  Full GCS + CLF-CBF pipeline (plan, then follow).
% OFFLINE (once, before tracking):
%   1. Inflate the obstacles by `margin` and decompose the free space
%      (free_space_decomposition); plan the crossing points with the SPP on
%      the intersection graph (build_intersection_graph). The safe sets below
%      are NOT used by the planner.
%   2. Safe-set chain (build_safe_chain): region, FIRI bridge, region, ...
%      Each bridge is a FIRI polytope (Wang et al.) seeded on the facet
%      between two consecutive regions and confined to their bounding box
%      (Gonzalez-Garcia et al. heuristic), so consecutive safe sets overlap.
%   3. Targets (translate_waypoint): each crossing point translated along
%      n_qhat into the intersection of the two adjacent safe sets with
%      margin mu (Matias & Silvestre Eqs. 32, 42-44 adapted).
% ONLINE:
%   4. Double-integrator CLF-CBF backstepping tracker with set-entry
%      switching and after-next look-ahead (double_integrator_safe_chain_clf_cbf).
% Verification: chain overlaps, target margins, safety signals, clearance of
% the executed path from the TRUE obstacles, convergence.
% Requires: YALMIP + a conic solver, Optimization Toolbox.

clear; clc;

q0 = [0.1; 0.2];
qT = [9.9; 9.9];
margin = 0.25;     % clearance from the true obstacles (planning and CBF safe sets)
bbox = [0 0; 10 0; 10 10; 0 10];
obstacles = { [2 2; 4 2; 4 3; 2 3], ...                    % L-shape (part 1)
              [3 3; 3 5; 2 5; 2 3], ...                    % L-shape (part 2)
              [6 6; 8 6; 8 7; 6 7], ...                    % L-shape (part 1)
              [7 7; 7 9; 6 9; 6 7], ...                    % L-shape (part 2)
              [1 7; 2 7; 2.5 8; 2 9; 1 9; 0.5 8], ...      % convex hexagon
              [6 1; 8 1; 8 2; 6 2], ...                    % L-shape (part 1)
              [8 2; 8 3; 7 3; 7 2] };                      % L-shape (part 2)
P = safe_chain_params();

%% Offline: planning
tOff = tic;
inflated = inflate_obstacles(obstacles, margin);
regions = sort_regions(free_space_decomposition(bbox, inflated), q0);
neighbors = adjacent_regions(regions);
regions_cs = cellfun(@(v) ConvexSet.polygon(v), regions, 'UniformOutput', false);
g = build_intersection_graph(q0, qT, regions_cs, neighbors);
[path, positions, cost, info] = ShortestPathProblem(g, 's', 't').solve(false);
seq = intersection_path_regions(path, regions_cs, q0);
W = positions(:, 2:end-1);                    % crossing points, W(:,j) on facet j
tPlan = toc(tOff);

%% Offline: safe-set chain and targets
chain = build_safe_chain(regions, seq, convex_obstacle_pieces(inflated), 'firi');
[targets, tinfo, tTargets] = chain_targets(chain, W, qT, P);
tOffline = toc(tOff);

%% Online: tracking
out = double_integrator_safe_chain_clf_cbf(chain.sets, targets, q0, [0; 0], P);

%% Report
clr = path_clearance(out.x, obstacles);
fprintf('---------Planning (offline)----------\n');
fprintf('Path: %s | regions %s | cost %.4f | gap %.2e\n', strjoin(path, ' -> '), ...
    mat2str(seq), cost, info.optimality_gap);
fprintf('Times: planning %.2f s, safe sets (FIRI) %.2f s, targets %.3f s, offline total %.2f s\n', ...
    tPlan, chain.time, tTargets, tOffline);
fr = [tinfo(1:end-1).facetRow];
fprintf(['Chain: %d safe sets | target margins min %.3f (mu %.2f), goal margin %.3f | ' ...
    'translations: n_qhat %d, bisector %d, Chebyshev fallback %d, none needed %d\n'], ...
    numel(chain.sets), min([tinfo(1:end-1).margin]), P.mu, tinfo(end).margin, ...
    sum(fr > 0), sum(fr == -1), sum([tinfo.fallback]), sum(fr == 0 & ~[tinfo(1:end-1).fallback]));
fprintf('---------Tracking (online)----------\n');
fprintf('Wall time %.2f s | arrival %.2f s | jumps %d/%d | converged %d | final error %.2e\n', ...
    out.wallTime, out.t(end), out.verification.nJumps, numel(chain.sets) - 1, ...
    out.verification.converged, out.verification.finalError);
fprintf('min h1 %.3e | min margin in active set %.3e | QP fallbacks %d | backstep viol %.2e\n', ...
    out.verification.minH1, out.verification.minMargin, out.verification.nQpFallback, ...
    out.verification.maxBackstepViolation);
fprintf('Min clearance of the executed path from the true obstacles: %.4f (margin %.2f)\n', clr, margin);
if clr < 0
    warning('demo_gcs_clf_cbf_following:collision', 'Executed path enters a true obstacle.');
end

%% Plot
figure('Color', 'w');  hold on;  axis equal;  box on;
for k = 1:numel(regions)
    fill(regions{k}(:,1), regions{k}(:,2), [0.9 0.93 1], 'EdgeColor', [0.7 0.75 0.9]);
end
for j = 1:numel(chain.bridgesV)
    B = chain.bridgesV{j};
    fill(B(:,1), B(:,2), [1 0.85 0.4], 'FaceAlpha', 0.35, 'EdgeColor', [0.85 0.55 0]);
end
for k = 1:numel(obstacles)
    fill(obstacles{k}(:,1), obstacles{k}(:,2), [0.2 0.2 0.2], 'EdgeColor', 'none');
end
h1 = plot(positions(1,:), positions(2,:), 'r--o', 'LineWidth', 1);
h2 = plot(targets(1,:), targets(2,:), 'kp', 'MarkerFaceColor', 'y', 'MarkerSize', 11);
h3 = plot(out.x(1,:), out.x(2,:), 'b-', 'LineWidth', 2);
J = out.jumps;
for k = 1:size(J, 1)
    [~, i] = min(abs(out.t - J(k, 1)));
    plot(out.x(1, i), out.x(2, i), 'mx', 'MarkerSize', 9, 'LineWidth', 2);
end
xlim([0 10]);  ylim([0 10]);  xlabel('x');  ylabel('y');
legend([h1 h2 h3], 'SPP crossing points', 'translated targets', 'executed', 'Location', 'southeast');
title(sprintf('Safe-set chain (%d sets: regions + FIRI bridges), jumps (x)', numel(chain.sets)));
