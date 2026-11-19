function [targets, info, time] = gcs_trajectory_targets(chain, q, qT, P)
%GCS_TRAJECTORY_TARGETS  Secondary version: waypoints of a safe-set chain
%   obtained by intersecting a GCS trajectory (GcsMotionPlanner) with the
%   q_hat hyperplanes, then turned into targets as in the main version.
%
%   [targets, info, time] = gcs_trajectory_targets(chain, q, qT, P)
%
%   chain   build_safe_chain output along the GCS region sequence
%   q       K x n samples of the GCS trajectory (info.trajectory.sample.q),
%           ordered in time, inside the union of the regions
%   qT      goal;  P  safe_chain_params (mu, kappa)
%
%   Waypoint of the pair (S_i, S_{i+1}): q_hat is the half-space of S_{i+1}
%   through which the trajectory enters S_{i+1} (the row most violated at the
%   sample before the entry: the forecast half-space of Matias & Silvestre
%   Eq. 32, which becomes active there); the waypoint w_i is the intersection
%   of the trajectory with the hyperplane a_qhat'x = b_qhat (linear
%   interpolation between samples; the entry sample itself if the
%   trajectory starts inside S_{i+1}). Waypoints are monotone in time.
%   Target: translate_waypoint(w_i, w_{i+1}, {S_i, S_{i+1}}, mu, kappa), the
%   same translation along n_qhat as the main version (the GCS trajectory
%   touches the obstacle corners, so w_i itself has no margin there).
%
%   info(i): sample (entry index), qhat (row of S_{i+1}), waypoint, and the
%   translate_waypoint fields (tau, direction, facetSet, facetRow, fallback,
%   margin).

t0 = tic;
N = numel(chain.sets);  [K, n] = size(q);
Wp = zeros(n, N-1);  ks = zeros(1, N-1);  qs = zeros(1, N-1);
k0 = 1;
for i = 1:N-1
    Sn = chain.sets{i+1};
    inN = all(q * Sn.A' - Sn.b' <= 1e-9, 2);
    ke = find(inN & (1:K)' >= k0, 1);
    if isempty(ke)
        error('gcs_trajectory_targets:entry', 'The trajectory never enters safe set %d.', i+1);
    end
    if ke > k0
        [~, qh] = max(Sn.A * q(ke-1, :)' - Sn.b);
        s1 = Sn.A(qh, :) * q(ke-1, :)' - Sn.b(qh);  s2 = Sn.A(qh, :) * q(ke, :)' - Sn.b(qh);
        lam = s1 / (s1 - s2);
        Wp(:, i) = q(ke-1, :)' + lam * (q(ke, :) - q(ke-1, :))';
    else
        [~, qh] = min(Sn.b - Sn.A * q(ke, :)');
        Wp(:, i) = q(ke, :)';
    end
    ks(i) = ke;  qs(i) = qh;  k0 = ke;
end
targets = zeros(n, N);
info = struct('sample', {}, 'qhat', {}, 'waypoint', {}, 'tau', {}, 'direction', {}, ...
    'facetSet', {}, 'facetRow', {}, 'fallback', {}, 'margin', {});
for i = 1:N-1
    if i < N-1, xbar = Wp(:, i+1); else, xbar = qT(:); end
    [targets(:, i), ti] = translate_waypoint(Wp(:, i), xbar, chain.sets([i i+1]), P.mu, P.kappa);
    info(i) = struct('sample', ks(i), 'qhat', qs(i), 'waypoint', Wp(:, i), 'tau', ti.tau, ...
        'direction', ti.direction, 'facetSet', ti.facetSet, 'facetRow', ti.facetRow, ...
        'fallback', ti.fallback, 'margin', ti.margin);
end
targets(:, N) = qT(:);
info(N) = struct('sample', K, 'qhat', 0, 'waypoint', qT(:), 'tau', 0, 'direction', zeros(n, 1), ...
    'facetSet', 0, 'facetRow', 0, 'fallback', false, ...
    'margin', lse_cbf(qT(:), chain.sets{N}.A, chain.sets{N}.b, P.kappa));
time = toc(t0);
end
