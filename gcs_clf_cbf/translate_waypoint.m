function [g, info] = translate_waypoint(w, xbar, sets, mu, kappa)
%TRANSLATE_WAYPOINT  Move a planned waypoint into the interior of the
%   intersection of adjacent safe convex sets.
%
%   [g, info] = translate_waypoint(w, xbar, sets, mu, kappa)
%
%   w      n x 1 planned waypoint (typically on a facet of the sets)
%   xbar   n x 1 waypoint to be reached after w (the after-next waypoint)
%   sets   cell of struct A, b (unit rows), e.g. {S_i, S_{i+1}}
%   mu     required margin: h_S(g) >= mu for every S in sets (LSE CBF)
%   kappa  LSE sharpness (safe_chain_params)
%
%   g = w + tau * n_qhat, adapting Matias & Silvestre Eqs. 30-33 and 42-44:
%   - reference direction v = (xbar - w)/|xbar - w| (Eq. 30 uses the normal
%     of the half-space containing x_bar; here the target to reach is a point
%     inside a bounded set, so the direction to it is used);
%   - candidate half-spaces q: facets of the sets that w lies on or near
%     (slack < mu + ln(m)/kappa), with inward unit normal n_q; they are tried
%     in decreasing alignment v'n_q (Eqs. 32-33 select the facet whose
%     normal is better aligned with the reference direction);
%   - tau = smallest tau >= 0 with min_S h_S(w + tau n_q) >= mu (Eqs. 42-44
%     with the synergy gap mu; h_S is concave, so along the ray the feasible
%     tau form an interval, found by golden-section + bisection).
%   If no candidate normal admits a feasible tau (w at a corner where a
%   single normal leaves the sets before reaching the margin), the bisector
%   of the candidate normals is tried (info.facetRow = -1), then the
%   direction towards the Chebyshev centre of the intersection
%   (info.fallback = true).
%
%   info: tau, direction, facetSet / facetRow (the selected q_hat; 0 / -1 for
%   the bisector, 0 / 0 for the fallback), fallback, margin (min_S h_S(g)).

w = w(:);  xbar = xbar(:);  n = numel(w);
phi = @(p) min(cellfun(@(S) lse_cbf(p, S.A, S.b, kappa), sets));
if phi(w) >= mu
    g = w;
    info = struct('tau', 0, 'direction', zeros(n, 1), 'facetSet', 0, 'facetRow', 0, ...
        'fallback', false, 'margin', phi(w));
    return
end
v = xbar - w;
if norm(v) > 1e-12, v = v / norm(v); end

cand = zeros(0, 2);  N = zeros(n, 0);  score = zeros(1, 0);
for k = 1:numel(sets)
    S = sets{k};
    slack = S.b - S.A * w;
    thr = mu + log(size(S.A, 1)) / kappa;
    for r = find(slack < thr)'
        nq = -S.A(r, :)';
        cand(end+1, :) = [k r]; %#ok<AGROW>
        N(:, end+1) = nq; %#ok<AGROW>
        score(end+1) = v' * nq; %#ok<AGROW>
    end
end
[~, order] = sort(score, 'descend');
if size(N, 2) > 1                               % corner: bisector of the active normals
    nb = sum(N, 2);
    if norm(nb) > 1e-9
        N(:, end+1) = nb / norm(nb);  cand(end+1, :) = [0 -1];  order(end+1) = size(N, 2);
    end
end
for c = order
    [tau, ok] = ray_search(w, N(:, c), sets, phi, mu);
    if ok
        g = w + tau * N(:, c);
        info = struct('tau', tau, 'direction', N(:, c), 'facetSet', cand(c, 1), ...
            'facetRow', cand(c, 2), 'fallback', false, 'margin', phi(g));
        return
    end
end

% Fallback: direction to the Chebyshev centre of the intersection.
Aall = cell2mat(cellfun(@(S) S.A, sets(:), 'UniformOutput', false));
ball = cell2mat(cellfun(@(S) S.b, sets(:), 'UniformOutput', false));
lpo = optimoptions('linprog', 'Display', 'off');
[xc, ~, flag] = linprog([zeros(n, 1); -1], [Aall, ones(size(Aall, 1), 1)], ball, [], [], [], [], lpo);
if flag ~= 1
    error('translate_waypoint:empty', 'The safe sets have no common interior.');
end
d = xc(1:n) - w;
[tau, ok] = ray_search(w, d / norm(d), sets, phi, mu);
if ~ok
    error('translate_waypoint:margin', ...
        'No point of the intersection reaches the margin mu = %.3g (max %.3g): reduce mu or increase kappa.', ...
        mu, phi(xc(1:n)));
end
g = w + tau * d / norm(d);
info = struct('tau', tau, 'direction', d / norm(d), 'facetSet', 0, 'facetRow', 0, ...
    'fallback', true, 'margin', phi(g));
end

% =========================================================================
function [tau, ok] = ray_search(w, d, sets, phi, mu)
% Smallest tau >= 0 with phi(w + tau d) >= mu; phi concave along the ray.
tmax = inf;                                     % the ray leaves every set
for k = 1:numel(sets)
    ad = sets{k}.A * d;  sl = sets{k}.b - sets{k}.A * w;
    pos = ad > 1e-12;
    if any(pos), tmax = min(tmax, max(0, min(sl(pos) ./ ad(pos)))); end
end
if ~isfinite(tmax)
    error('translate_waypoint:unbounded', 'The safe sets must be bounded.');
end
f = @(t) phi(w + t * d);
a = 0;  b = tmax;  gr = (sqrt(5) - 1) / 2;      % golden section: argmax of f
c1 = b - gr * (b - a);  c2 = a + gr * (b - a);  f1 = f(c1);  f2 = f(c2);
for it = 1:80
    if f1 < f2, a = c1; c1 = c2; f1 = f2; c2 = a + gr * (b - a); f2 = f(c2);
    else,       b = c2; c2 = c1; f2 = f1; c1 = b - gr * (b - a); f1 = f(c1);
    end
end
tbest = (a + b) / 2;
ok = f(tbest) >= mu;
tau = NaN;
if ~ok, return; end
lo = 0;  hi = tbest;                            % f increasing on [0, tbest]
for it = 1:80
    mid = (lo + hi) / 2;
    if f(mid) >= mu, hi = mid; else, lo = mid; end
end
tau = hi;
end
