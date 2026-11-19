%% COMPARE_SAFE_SET_VERSIONS  Main safe-chain tracker vs its two secondary versions.
%   MAIN    intersection-graph SPP crossing points -> offline FIRI bridges
%           (build_safe_chain 'firi') -> targets translated along n_qhat
%           (chain_targets) -> double_integrator_safe_chain_clf_cbf
%   ONLINE  same SPP crossing points -> targets in the seed quadrilaterals
%           (build_safe_chain 'seed', no FIRI) -> FIRI safe set computed
%           online at every mode (online_firi/double_integrator_online_firi_clf_cbf)
%   GCS     GcsMotionPlanner (gcs_pid) trajectory -> offline FIRI bridges along
%           its region sequence -> targets where the trajectory crosses the
%           q_hat hyperplanes (gcs_bezier/gcs_trajectory_targets) ->
%           double_integrator_safe_chain_clf_cbf
%   Same map, start, goal, inflation margin and controller parameters.
%   OFFLINE time (planning + safe sets + targets) and ONLINE time (tracking
%   simulation wall-clock, including online FIRI) are measured separately.
% Requires: YALMIP + a conic solver, Optimization Toolbox.

clear; clc;

q0 = [0.3; 0.4];  qT = [9.7; 9.6];  z0 = [0; 0];   % start/goal off the diagonal facets
margin = 0.25;
bbox = [0 0; 10 0; 10 10; 0 10];  workspace = [0 10; 0 10];
obstacles = { [2 2; 4 2; 4 3; 2 3], [3 3; 3 5; 2 5; 2 3], ...      % L-shape
              [6 6; 8 6; 8 7; 6 7], [7 7; 7 9; 6 9; 6 7], ...      % L-shape
              [1 7; 2 7; 2.5 8; 2 9; 1 9; 0.5 8], ...              % hexagon
              [6 1; 8 1; 8 2; 6 2], [8 2; 8 3; 7 3; 7 2] };        % L-shape
P = safe_chain_params(struct('kappa', 80, 'mu', 0.15));
onlineOpts = struct('delta', 0.1, 'boxPad', 2.0, 'dtCheck', 0.05);
gcsParams = struct('a', 1, 'b', 0, 'c', 0, 'eps', 1e-1, 'eta', 1, 'd', 3, 'vmax', 1, ...
    'hdotMin', 1e-1, 'v0', z0, 'vT', [0; 0], 'Tmax', 100, 'nSample', 200);
names = {'main (offline FIRI)', 'online FIRI', 'GCS Bezier'};
R = struct('name', names, 'tPlan', 0, 'tSets', 0, 'tTargets', 0, 'tOffline', 0, ...
    'tOnline', 0, 'tFiri', 0, 'firiCalls', 0, 'nSets', 0, 'out', [], 'targets', [], 'setsV', []);

%% Warm-up (untimed): one-off first-call overheads of YALMIP, the solvers
% (quadprog, coneprog, linprog) and the JIT would otherwise be charged to
% whichever version runs first.
[regions, pieces, seq] = spp_plan(bbox, obstacles, margin, q0, qT);
build_safe_chain(regions, seq(1:2), pieces, 'firi');

%% MAIN
t0 = tic;
[regions, pieces, seq, W] = spp_plan(bbox, obstacles, margin, q0, qT);
R(1).tPlan = toc(t0);
chain = build_safe_chain(regions, seq, pieces, 'firi');
[R(1).targets, ~, R(1).tTargets] = chain_targets(chain, W, qT, P);
R(1).tSets = chain.time;  R(1).tOffline = toc(t0);  R(1).nSets = numel(chain.sets);
R(1).setsV = chain.bridgesV;
R(1).out = double_integrator_safe_chain_clf_cbf(chain.sets, R(1).targets, q0, z0, P);
R(1).tOnline = R(1).out.wallTime;

%% ONLINE
t0 = tic;
[regions, pieces, seq, W] = spp_plan(bbox, obstacles, margin, q0, qT);
R(2).tPlan = toc(t0);
chainSeed = build_safe_chain(regions, seq, pieces, 'seed');
[R(2).targets, ~, R(2).tTargets] = chain_targets(chainSeed, W, qT, P);
R(2).tSets = 0;  R(2).tOffline = toc(t0);
R(2).out = double_integrator_online_firi_clf_cbf(R(2).targets, pieces, workspace, q0, z0, P, onlineOpts);
R(2).tOnline = R(2).out.wallTime;  R(2).tFiri = R(2).out.firiTime;
R(2).firiCalls = R(2).out.firiCalls;  R(2).nSets = numel(R(2).out.sets);  R(2).setsV = R(2).out.setsV;

%% GCS
t0 = tic;
inflated = inflate_obstacles(obstacles, margin);
regions = free_space_decomposition(bbox, inflated);
pieces = convex_obstacle_pieces(inflated);
rcs = cellfun(@(v) ConvexSet.polygon(v), regions, 'UniformOutput', false);
[gpath, ~, ~, ginfo] = GcsMotionPlanner(rcs, q0, qT, gcsParams).solve(false);
seqG = cellfun(@str2double, gpath(2:end-1));
R(3).tPlan = toc(t0);
chainG = build_safe_chain(regions, seqG, pieces, 'firi');
[R(3).targets, ~, R(3).tTargets] = gcs_trajectory_targets(chainG, ginfo.trajectory.sample.q, qT, P);
R(3).tSets = chainG.time;  R(3).tOffline = toc(t0);  R(3).nSets = numel(chainG.sets);
R(3).setsV = chainG.bridgesV;
R(3).out = double_integrator_safe_chain_clf_cbf(chainG.sets, R(3).targets, q0, z0, P);
R(3).tOnline = R(3).out.wallTime;

%% Comparison table
fprintf('\n%-22s %9s %9s %9s %9s | %9s %9s %6s | %5s %5s %8s %8s %9s %9s %9s %4s\n', 'version', ...
    'plan[s]', 'sets[s]', 'targ[s]', 'OFFL[s]', 'ONL[s]', 'firi[s]', 'calls', ...
    'sets', 'jumps', 'arrive', 'length', 'clear', 'min h1', 'max|u|', 'conv');
for k = 1:3
    o = R(k).out;  v = o.verification;
    L = sum(vecnorm(diff(o.x, 1, 2)));
    clr = path_clearance(o.x, obstacles);
    fprintf('%-22s %9.3f %9.3f %9.3f %9.3f | %9.3f %9.3f %6d | %5d %5d %8.2f %8.3f %9.4f %9.2e %9.3f %4d\n', ...
        R(k).name, R(k).tPlan, R(k).tSets, R(k).tTargets, R(k).tOffline, R(k).tOnline, ...
        R(k).tFiri, R(k).firiCalls, R(k).nSets, v.nJumps, o.t(end), L, clr, v.minH1, ...
        max(vecnorm(o.u)), v.converged);
    if clr < margin - 1e-6
        warning('compare_safe_set_versions:clearance', '%s: clearance %.4f below the margin %.2f.', ...
            R(k).name, clr, margin);
    end
end
fprintf(['OFFL = plan + sets + targets (before tracking); ONL = tracking wall-clock ' ...
    '(online FIRI included, firi = its share); clear = distance to the TRUE obstacles.\n']);

%% Plot
figure('Color', 'w', 'Position', [100 100 1500 480]);
for k = 1:3
    subplot(1, 3, k);  hold on;  axis equal;  box on;
    for j = 1:numel(R(k).setsV)
        B = R(k).setsV{j};
        fill(B(:,1), B(:,2), [1 0.85 0.4], 'FaceAlpha', 0.2, 'EdgeColor', [0.85 0.55 0]);
    end
    for j = 1:numel(obstacles)
        fill(obstacles{j}(:,1), obstacles{j}(:,2), [0.2 0.2 0.2], 'EdgeColor', 'none');
    end
    if k == 3
        plot(ginfo.trajectory.sample.q(:,1), ginfo.trajectory.sample.q(:,2), 'r--', 'LineWidth', 1);
    end
    plot(R(k).targets(1,:), R(k).targets(2,:), 'kp', 'MarkerFaceColor', 'y', 'MarkerSize', 9);
    plot(R(k).out.x(1,:), R(k).out.x(2,:), 'b-', 'LineWidth', 2);
    xlim([0 10]);  ylim([0 10]);
    title(sprintf('%s: arrival %.1f s', R(k).name, R(k).out.t(end)));
end

% =========================================================================
function [regions, pieces, seq, W] = spp_plan(bbox, obstacles, margin, q0, qT)
% Planning of MAIN and ONLINE: inflated decomposition + SPP on the
% intersection graph; crossing points W and region sequence seq.
inflated = inflate_obstacles(obstacles, margin);
regions = sort_regions(free_space_decomposition(bbox, inflated), q0);
pieces = convex_obstacle_pieces(inflated);
rcs = cellfun(@(v) ConvexSet.polygon(v), regions, 'UniformOutput', false);
g = build_intersection_graph(q0, qT, rcs, adjacent_regions(regions));
[path, positions] = ShortestPathProblem(g, 's', 't').solve(false);
seq = intersection_path_regions(path, rcs, q0);
W = positions(:, 2:end-1);
end
