%% DEMO_GCS_PID_PLAN_TRACK  Full GCS + PID pipeline (plan, then track).
% Plan with GcsMotionPlanner (min time, C^1, |qdot|_inf <= 1) around the square
% obstacle [2,3]^2, then track the planned trajectory with track_gcs_pid
% from an offset initial state.
% Requires: YALMIP + a conic solver, Optimization Toolbox.

clear; clc;

regions = { ConvexSet.box([0;0], [5;2]), ConvexSet.box([0;3], [5;5]), ...
            ConvexSet.box([0;0], [2;5]), ConvexSet.box([3;0], [5;5]) };
q0 = [0.2; 0.2];  qT = [4.8; 4.8];

% ---- Plan -----------------------------------------------------------------
params = struct('a', 1, 'b', 0, 'c', 0, 'eps', 1e-1, 'eta', 1, 'd', 3, 'vmax', 1, ...
              'hdotMin', 1e-1, 'v0', [0;0], 'vT', [0;0], 'Tmax', 100);
planner = GcsMotionPlanner(regions, q0, qT, params);
[path, ~, ~, info] = planner.solve(false);
ref = info.trajectory.sample;
assert(info.trajectory.verification.pass, 'Planned trajectory failed verification.');
fprintf('Plan: path %s, duration %.3f, certified gap %.3f %%\n', ...
    strjoin(path, ' -> '), info.trajectory.verification.duration, 100 * info.relative_gap);

% ---- Track (initial offset of 0.1 in x and -0.1 in y) ---------------------
gains = struct('Kp', 25, 'Ki', 5, 'Kd', 10);
x0 = [q0 + [0.1; -0.1]; 0; 0];
out = track_gcs_pid(ref, gains, regions, x0);
fprintf('Tracking: rms error %.4f, max error %.4f, max region violation %.2e\n', ...
    out.rms_error, out.max_error, out.max_region_violation);

% ---- Plot -----------------------------------------------------------------
figure('Color', 'w');
subplot(1, 2, 1);  hold on;  axis equal;  grid on;
for i = 1:numel(regions)
    cs = regions{i};  lb = -cs.b(3:4);  ub = cs.b(1:2);
    rectangle('Position', [lb(1) lb(2) ub(1)-lb(1) ub(2)-lb(2)], ...
        'FaceColor', [0.85 0.92 1], 'EdgeColor', [0.4 0.6 0.9]);
end
rectangle('Position', [2 2 1 1], 'FaceColor', [0.9 0.2 0.2], 'EdgeColor', 'none');
h1 = plot(ref.q(:,1), ref.q(:,2), 'b-', 'LineWidth', 2);
h2 = plot(out.q(:,1), out.q(:,2), 'm--', 'LineWidth', 1.5);
h3 = plot([q0(1) qT(1)], [q0(2) qT(2)], 'kx', 'MarkerSize', 10, 'LineWidth', 2);
xlim([0 5]);  ylim([0 5]);  xlabel('x');  ylabel('y');
legend([h1 h2 h3], 'GCS reference', 'PID executed', 'start / goal', 'Location', 'best');
title('Reference vs executed path');

subplot(1, 2, 2);
plot(out.t, sqrt(sum(out.e.^2, 2)), 'LineWidth', 1.5);
grid on;  xlabel('t');  ylabel('||q_{ref} - q||');
title('Tracking error');
