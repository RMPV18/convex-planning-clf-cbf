%% DEMO_PLAN_AND_TRACK  dcp_polygon_planner -> pid_trajectory_tracker, no files exchanged.
clear; clc; close all;
if ~exist('cvx_begin', 'file'), cvx_setup; end

x0 = [0; 0];
xf = [6; 4];
obs = { [1.5 0.5; 2.5 0.5; 2.5 1.5; 1.5 1.5], ...                 % square
        [3.5 2.0; 4.8 2.2; 4.2 3.4], ...                          % triangle
        [3.0 -0.2; 4.0 -0.2; 4.5 0.6; 3.5 1.2; 2.8 0.6] };        % pentagon

planPar.N = 60; planPar.Ts = 0.1; planPar.vmax = 2.0; planPar.margin = 0.15;
plan = dcp_polygon_planner(x0, xf, obs, planPar);
fprintf('plan: feasible=%d collisions=%d length=%.4f iters=%d time=%.2f s\n', ...
    plan.feasible, plan.collisions, plan.L, plan.iters, plan.time);
assert(plan.feasible, 'Planner returned a path that intersects an inflated obstacle.');

trkPar.Kp = 25; trkPar.Ki = 2; trkPar.Kd = 10;
out = pid_trajectory_tracker(plan, trkPar, obs);
fprintf('track: RMS %.4f | max %.4f | final %.4f | min clearance %.4f (margin %.4f)\n', ...
    sqrt(mean(out.err.^2)), max(out.err), out.err(end), min(out.clearance), plan.par.margin);
if min(out.clearance) <= 0
    warning('demo:collision', 'Tracked trajectory enters an obstacle.');
elseif min(out.clearance) < plan.par.margin
    warning('demo:margin', 'Tracked trajectory uses %.1f%% of the planner clearance margin.', ...
        100 * (plan.par.margin - min(out.clearance)) / plan.par.margin);
end

figure; hold on; axis equal; grid on;
for j = 1:numel(obs), fill(obs{j}(:,1), obs{j}(:,2), [0.85 0.85 0.85]); end
plot(out.pr(1,:), out.pr(2,:), 'k--', 'LineWidth', 1.2);
plot(out.p(1,:),  out.p(2,:),  'b-',  'LineWidth', 1.5);
legend('', '', '', 'reference', 'tracked'); title('Planned vs tracked');

figure;
subplot(2,1,1); plot(out.t, out.err); grid on; ylabel('|e|');
subplot(2,1,2); plot(out.t(1:end-1), out.u'); grid on; ylabel('u'); xlabel('t [s]');
