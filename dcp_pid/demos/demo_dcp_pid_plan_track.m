%% DEMO_DCP_PID_PLAN_TRACK  Full DCP + PID pipeline (plan, then track).
% Plans a collision-free reference among three convex polygons with
% dcp_polygon_planner (penalty convex-concave procedure), tracks it with
% pid_trajectory_tracker (point mass, saturated discrete PID) and reports
% tracking error and obstacle clearance.
% Requires: CVX (Grant & Boyd).
clear; clc; close all;
if ~exist('cvx_begin', 'file'), cvx_setup; end

q0 = [0; 0];
qT = [6; 4];
obstacles = { [1.5 0.5; 2.5 0.5; 2.5 1.5; 1.5 1.5], ...          % square
              [3.5 2.0; 4.8 2.2; 4.2 3.4], ...                   % triangle
              [3.0 -0.2; 4.0 -0.2; 4.5 0.6; 3.5 1.2; 2.8 0.6] }; % pentagon

plan_params.N = 60; plan_params.Ts = 0.1; plan_params.vmax = 2.0; plan_params.margin = 0.15;
plan = dcp_polygon_planner(q0, qT, obstacles, plan_params);
fprintf('plan: feasible=%d collisions=%d length=%.4f iters=%d time=%.2f s\n', ...
    plan.feasible, plan.collisions, plan.L, plan.iters, plan.time);
assert(plan.feasible, 'Planner returned a path that intersects an inflated obstacle.');

track_params.Kp = 25; track_params.Ki = 2; track_params.Kd = 10;
out = pid_trajectory_tracker(plan, track_params, obstacles);
fprintf('track: RMS %.4f | max %.4f | final %.4f | min clearance %.4f (margin %.4f)\n', ...
    sqrt(mean(out.err.^2)), max(out.err), out.err(end), min(out.clearance), plan.par.margin);
if min(out.clearance) <= 0
    warning('demo_dcp_pid_plan_track:collision', 'Tracked trajectory enters an obstacle.');
elseif min(out.clearance) < plan.par.margin
    warning('demo_dcp_pid_plan_track:margin', 'Tracked trajectory uses %.1f%% of the planner clearance margin.', ...
        100 * (plan.par.margin - min(out.clearance)) / plan.par.margin);
end

figure; hold on; axis equal; grid on;
for j = 1:numel(obstacles), fill(obstacles{j}(:,1), obstacles{j}(:,2), [0.85 0.85 0.85]); end
plot(out.pr(1,:), out.pr(2,:), 'k--', 'LineWidth', 1.2);
plot(out.p(1,:),  out.p(2,:),  'b-',  'LineWidth', 1.5);
legend('', '', '', 'reference', 'tracked'); title('Planned vs tracked');

figure;
subplot(2,1,1); plot(out.t, out.err); grid on; ylabel('|e|');
subplot(2,1,2); plot(out.t(1:end-1), out.u'); grid on; ylabel('u'); xlabel('t [s]');
