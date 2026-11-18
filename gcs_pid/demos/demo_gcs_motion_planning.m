%% DEMO_GCS_MOTION_PLANNING  Planning stage of the GCS + PID pipeline.
% GCS motion planning (Marcucci et al., arXiv:2205.04422) around a square
% obstacle [2,3]^2 in [0,5]^2, using ConvexSet / GraphOfConvexSets /
% ShortestPathProblem (via GcsMotionPlanner). Four overlapping boxes give an
% exact free-space cover with two homotopy classes (above / below).
%
% Solves two problems (min length; min time with C^1 and velocity limits).
% Requires: YALMIP + a conic solver (Mosek, Gurobi, ...), Optimization Toolbox.

clear; clc;

regions = { ConvexSet.box([0;0], [5;2]), ...   % bottom
            ConvexSet.box([0;3], [5;5]), ...   % top
            ConvexSet.box([0;0], [2;5]), ...   % left
            ConvexSet.box([3;0], [5;5]) };     % right
q0 = [0.2; 0.2];  qT = [4.8; 4.8];

% Problem 1: minimum Euclidean length, continuous path (eta=0, d=1).
o1 = struct('a', 0, 'b', 1, 'c', 0, 'eta', 0, 'd', 1);
% Problem 2: minimum time, |qdot|_inf <= 1, C^1, regularized acceleration.
o2 = struct('a', 1, 'b', 0, 'c', 0, 'eps', 1e-1, 'eta', 1, 'd', 3, 'vmax', 1, ...
            'hdotMin', 1e-1, 'v0', [0;0], 'vT', [0;0], 'Tmax', 100);
params = {o1, o2};  titles = {'min length', 'min time, C^1, |qdot|_\infty <= 1'};

figure('Color', 'w');
for k = 1:2
    planner = GcsMotionPlanner(regions, q0, qT, params{k});
    tic
    [path, positions, cost, info] = planner.solve(false);
    toc
    T = info.trajectory;

    fprintf('Problem %d: path %s\n', k, strjoin(path, ' -> '));
    fprintf('  Relaxation cost (lower bound): %.4f\n', info.relaxation_cost);
    fprintf('  Restriction cost (feasible):   %.4f\n', info.restriction_cost);
    fprintf('  Certified gap delta_relax:     %.3f %%\n', 100 * info.relative_gap);
    fprintf('  Distinct rounded paths: %d, costs: %s\n', numel(info.rounding_costs), ...
        mat2str(info.rounding_costs, 5));
    fprintf('  Verification: containment %.2e, continuity %.2e, pass %d\n', ...
        T.verification.maxRegionViolation, T.verification.maxContinuityGap, T.verification.pass);

    subplot(1, 2, k);  hold on;  axis equal;  grid on;
    for i = 1:numel(regions)
        cs = regions{i};  lb = -cs.b(3:4);  ub = cs.b(1:2);
        rectangle('Position', [lb(1) lb(2) ub(1)-lb(1) ub(2)-lb(2)], ...
            'FaceColor', [0.85 0.92 1], 'EdgeColor', [0.4 0.6 0.9]);
    end
    rectangle('Position', [2 2 1 1], 'FaceColor', [0.9 0.2 0.2], 'EdgeColor', 'none');
    plot(T.sample.q(:,1), T.sample.q(:,2), 'b-', 'LineWidth', 2);
    plot([q0(1) qT(1)], [q0(2) qT(2)], 'kx', 'MarkerSize', 10, 'LineWidth', 2);
    xlim([0 5]);  ylim([0 5]);
    title(sprintf('%s (gap %.2f %%)', titles{k}, 100 * info.relative_gap));
    xlabel('x'); ylabel('y');
    if k == 2, T2 = T; end
end

figure('Color', 'w');
plot(T2.sample.t, T2.sample.qdot, 'LineWidth', 1.5);
xlabel('t'); ylabel('qdot'); legend('qdot_x', 'qdot_y'); grid on;
