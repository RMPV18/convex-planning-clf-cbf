%% DEMO_GCS_PID_PLAN_TRACK  GCS + PID with several non-convex obstacles.
% 1. Inflate the obstacles by a clearance margin and decompose the remaining
%    free space into convex polygons (free_space_decomposition). The margin
%    must exceed the PID tracking error: the GCS reference touches the
%    boundary of the planning regions at obstacle corners.
% 2. Plan with GcsMotionPlanner (min time, C^1, |qdot|_inf <= 1).
% 3. Track with track_gcs_pid from an offset initial state and check the
%    executed path against the TRUE (uninflated) free space.
% Requires: YALMIP + a conic solver, Optimization Toolbox.

clear; clc;

bbox = [0 0; 10 0; 10 10; 0 10];
obstacles = { [2 2; 4 2; 4 3; 3 3; 3 5; 2 5], ...          % L-shape
              [6 6; 8 6; 8 7; 7 7; 7 9; 6 9], ...          % L-shape
              [1 7; 2 7; 2.5 8; 2 9; 1 9; 0.5 8], ...      % hexagon
              [6 1; 8 1; 8 3; 7 3; 7 2; 6 2] };            % L-shape
q0 = [0.3; 0.3];  qT = [9.7; 9.7];
margin = 0.25;                              % clearance > expected tracking error

% ---- Regions --------------------------------------------------------------
V = free_space_decomposition(bbox, inflate_obstacles(obstacles, margin));   % planning regions
regions = cellfun(@(v) ConvexSet.polygon(v), V, 'UniformOutput', false);
Vtrue = free_space_decomposition(bbox, obstacles);     % for collision check
trueFree = cellfun(@(v) ConvexSet.polygon(v), Vtrue, 'UniformOutput', false);

% ---- Plan -----------------------------------------------------------------
params = struct('a', 1, 'b', 0, 'c', 0, 'eps', 1e-1, 'eta', 1, 'd', 3, 'vmax', 1, ...
              'hdotMin', 1e-1, 'v0', [0;0], 'vT', [0;0], 'Tmax', 100);
planner = GcsMotionPlanner(regions, q0, qT, params);
[path, ~, ~, info] = planner.solve(false);
ref = info.trajectory.sample;
assert(info.trajectory.verification.pass, 'Planned trajectory failed verification.');
fprintf('Plan: %d regions, path %s, duration %.3f, certified gap %.2f %%\n', numel(regions), ...
    strjoin(path, ' -> '), info.trajectory.verification.duration, 100 * info.relative_gap);

% ---- Track ----------------------------------------------------------------
gains = struct('Kp', 25, 'Ki', 5, 'Kd', 10);
x0 = [q0 + [0.1; -0.1]; 0; 0];
out = track_gcs_pid(ref, gains, trueFree, x0);
fprintf('Tracking: rms error %.4f, max error %.4f (margin %.2f)\n', ...
    out.rms_error, out.max_error, margin);
fprintf('Max violation of the true free space: %.2e (<= 0: collision-free)\n', ...
    out.max_region_violation);
if out.max_region_violation > 0
    warning('demo_gcs_pid_plan_track:collision', ...
        'Executed path enters an obstacle; increase margin or the PID gains.');
end

% ---- Plot -----------------------------------------------------------------
figure('Color', 'w');
subplot(1, 2, 1);  hold on;  axis equal;  grid on;
for k = 1:numel(V)
    fill(V{k}(:,1), V{k}(:,2), [0.85 0.92 1], 'EdgeColor', [0.4 0.6 0.9]);
end
for k = 1:numel(obstacles)
    fill(obstacles{k}(:,1), obstacles{k}(:,2), [0.9 0.2 0.2], 'EdgeColor', 'none');
end
h1 = plot(ref.q(:,1), ref.q(:,2), 'b-', 'LineWidth', 2);
h2 = plot(out.q(:,1), out.q(:,2), 'm--', 'LineWidth', 1.5);
h3 = plot([q0(1) qT(1)], [q0(2) qT(2)], 'kx', 'MarkerSize', 10, 'LineWidth', 2);
xlim([0 10]);  ylim([0 10]);  xlabel('x');  ylabel('y');
legend([h1 h2 h3], 'GCS reference', 'PID executed', 'start / goal', 'Location', 'best');
title(sprintf('%d regions, margin %.2f', numel(regions), margin));

subplot(1, 2, 2);
plot(out.t, sqrt(sum(out.e.^2, 2)), 'LineWidth', 1.5);
grid on;  xlabel('t');  ylabel('||q_{ref} - q||');
title('Tracking error');
