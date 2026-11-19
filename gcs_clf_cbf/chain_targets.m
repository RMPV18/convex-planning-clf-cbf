function [targets, info, time] = chain_targets(chain, W, qT, P)
%CHAIN_TARGETS  Targets of a safe-set chain from the planned crossing points.
%
%   [targets, info, time] = chain_targets(chain, W, qT, P)
%
%   chain    build_safe_chain output (N sets, pairs (S_i, S_{i+1}) belong
%            to crossing chain.crossing(i))
%   W        n x (M-1) planned crossing points, W(:,j) on facet j
%   qT       goal (last target)
%   P        safe_chain_params (mu, kappa)
%
%   targets(:,i), i < N: W(:,j) translated into S_i ∩ S_{i+1} with margin mu
%   (translate_waypoint), the reference direction pointing to the waypoint
%   after it (W(:,j+1), or qT for the last crossing); targets(:,N) = qT.
%   info(i): translate_waypoint info;  time: wall-clock time.

t0 = tic;
N = numel(chain.sets);
targets = zeros(size(W, 1), N);
info = struct('tau', {}, 'direction', {}, 'facetSet', {}, 'facetRow', {}, 'fallback', {}, 'margin', {});
for i = 1:N-1
    j = chain.crossing(i);
    if j < size(W, 2), xbar = W(:, j+1); else, xbar = qT(:); end
    [targets(:, i), info(i)] = translate_waypoint(W(:, j), xbar, chain.sets([i i+1]), P.mu, P.kappa);
end
targets(:, N) = qT(:);
hT = lse_cbf(qT(:), chain.sets{N}.A, chain.sets{N}.b, P.kappa);
info(N) = struct('tau', 0, 'direction', zeros(size(W, 1), 1), 'facetSet', 0, 'facetRow', 0, ...
    'fallback', false, 'margin', hT);
if hT < 0
    error('chain_targets:goal', 'The goal is outside the LSE safe set of the last region (h = %.3g).', hT);
end
time = toc(t0);
end
