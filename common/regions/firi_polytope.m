function [A, b, info] = firi_polytope(seedV, obstacles, boxA, boxb, opts)
%FIRI_POLYTOPE  Fast Iterative Region Inflation: large obstacle-free convex
%   polytope that contains a convex seed (2-D / 3-D).
%
%   [A, b, info] = firi_polytope(seedV, obstacles, boxA, boxb)
%   [A, b, info] = firi_polytope(seedV, obstacles, boxA, boxb, opts)
%
%   seedV      k x n   vertices of the seed Q = conv(seedV)  (n = 2 or 3);
%                      Q must be obstacle-free (touching is allowed)
%   obstacles  cell of m_i x n vertex arrays, each a CONVEX obstacle conv(U_i)
%   boxA, boxb bounding box / region of interest {x : boxA x <= boxb}; its
%              facets are the initial hyperplanes (heuristic of Gonzalez-Garcia
%              et al., PiP Sec. "Heuristic for Task-Informed Polytope
%              Inflation"); it must contain the seed
%   opts       struct, optional fields (defaults):
%                maxIter (10), volTol (0.02 = stop when the MVIE volume grows
%                less than 2 %), tol (1e-9)
%
%   A, b       polytope P = {x : A x <= b}, unit-norm rows (box rows first)
%   info       iterations, volume (MVIE volume history, up to a constant),
%              ellipsoid (L, d: E = {L u + d : |u| <= 1}), nHalfspaces,
%              seedViolation (max over seed vertices of A v - b, <= tol),
%              obstacleOverlap (max Chebyshev radius of P ∩ conv(U_i), <= tol:
%              no obstacle interior inside P), time
%
%   Method: Wang et al., "Fast Iterative Region Inflation for Computing Large
%   2-D/3-D Convex Regions of Obstacle-Free Space", IEEE T-RO 41 (2025),
%   Algorithm 1, alternating
%     (a) restrictive inflation RsI (Sec. IV-A): in the coordinates of the
%         current ellipsoid (u = L^{-1}(x - d)) every obstacle gets the
%         halfspace {b'u <= 1} farthest from the origin that contains the
%         seed and excludes the obstacle, Eq. (15):
%             min |b|^2  s.t.  v'b <= 1 (seed vertices),  w'b >= 1 (obstacle)
%         (inequality directions reconstructed from the text, see KB §14);
%         halfspaces are added greedily nearest-first and every obstacle
%         already excluded by an added halfspace is dropped (as in IRIS);
%     (b) maximum-volume inscribed ellipsoid (Sec. V) as an SOCP in the
%         Cholesky factor (KB §15 cheat-sheet), solved with coneprog.
%   The SDMN solver of Sec. IV-B is replaced by quadprog (same QP, slower).
%   Every iterate contains the seed and is obstacle-free (Eq. 12); the MVIE
%   volume is non-decreasing (Eq. 13).

if nargin < 5, opts = struct(); end
def = struct('maxIter', 10, 'volTol', 0.02, 'tol', 1e-9);
fn = fieldnames(def);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}), opts.(fn{i}) = def.(fn{i}); end
end
t0 = tic;
n = size(seedV, 2);
assert(any(n == [2 3]), 'firi_polytope:dim', 'Only n = 2 or 3 is supported (FIRI is a 2-D/3-D method).');
nr = vecnorm(boxA, 2, 2);
boxA = boxA ./ nr;  boxb = boxb(:) ./ nr;
if any(max(seedV * boxA' - boxb', [], 1) > opts.tol)
    error('firi_polytope:seed', 'The bounding box does not contain the seed.');
end

% Obstacles entirely outside the box can never cut the polytope.
keep = true(1, numel(obstacles));
for i = 1:numel(obstacles)
    keep(i) = ~any(all(obstacles{i} * boxA' - boxb' >= -opts.tol, 1));
end
obstacles = obstacles(keep);
No = numel(obstacles);

qpo = optimoptions('quadprog', 'Display', 'off', 'OptimalityTolerance', 1e-12, ...
    'ConstraintTolerance', 1e-12, 'StepTolerance', 1e-14);
d = mean(seedV, 1)';
L = 1e-3 * max(1, max(abs(seedV(:)))) * eye(n);
vol = [];
A = boxA;  b = boxb;
for it = 1:opts.maxIter
    % ---- (a) restrictive inflation in ellipsoid coordinates -----------------
    Li = inv(L);
    Vt = (seedV - d') * Li';
    bh = zeros(n, No);  dist = inf(1, No);  Ut = cell(1, No);
    for i = 1:No
        Ut{i} = (obstacles{i} - d') * Li';
        [bi, ~, flag] = quadprog(2 * eye(n), zeros(n, 1), [Vt; -Ut{i}], ...
            [ones(size(Vt, 1), 1); -ones(size(Ut{i}, 1), 1)], [], [], [], [], [], qpo);
        if flag ~= 1
            error('firi_polytope:separation', ...
                'No hyperplane separates the seed from obstacle %d (quadprog flag %d): the seed is not obstacle-free.', i, flag);
        end
        bh(:, i) = bi;  dist(i) = 1 / norm(bi);
    end
    [~, order] = sort(dist);
    active = true(1, No);
    H = zeros(0, n);  hb = zeros(0, 1);
    for i = order
        if ~active(i), continue; end
        bi = bh(:, i);
        a = Li' * bi;                      % b'L^{-1}(x - d) <= 1
        H(end + 1, :) = a' / norm(a); %#ok<AGROW>
        hb(end + 1, 1) = (1 + a' * d) / norm(a); %#ok<AGROW>
        for j = find(active)
            if all(Ut{j} * bi >= 1 - 1e-9), active(j) = false; end
        end
    end
    A = [boxA; H];  b = [boxb; hb];

    % ---- (b) maximum-volume inscribed ellipsoid ------------------------------
    [L, d, v] = mvie(A, b, n);
    vol(end + 1) = v; %#ok<AGROW>
    if it > 1 && (vol(end) - vol(end - 1)) / vol(end - 1) < opts.volTol
        break
    end
end

% ---- verification ---------------------------------------------------------
seedViol = max(max(seedV * A' - b', [], 2));
overlap = -inf;
lpo = optimoptions('linprog', 'Display', 'off');
for i = 1:No
    [Ao, bo] = vertex_hull_hrep(obstacles{i});
    AA = [A; Ao];  bb = [b; bo];
    [~, fv, flag] = linprog([zeros(n, 1); -1], [AA, vecnorm(AA, 2, 2)], bb, [], [], [], [], lpo);
    if flag == 1, overlap = max(overlap, -fv); end
end
if seedViol > opts.tol
    warning('firi_polytope:seed', 'Seed not contained in the polytope (violation %.2e).', seedViol);
end
if overlap > 1e-7
    warning('firi_polytope:collision', 'Polytope overlaps an obstacle (Chebyshev radius %.2e).', overlap);
end
info = struct('iterations', it, 'volume', vol, 'ellipsoid', struct('L', L, 'd', d), ...
    'nHalfspaces', size(A, 1), 'seedViolation', seedViol, 'obstacleOverlap', overlap, ...
    'time', toc(t0));
end

% =========================================================================
function [L, d, vol] = mvie(A, b, n)
% max det L  s.t.  |(A L)_i| <= b_i - a_i'd,  L lower triangular,
% via the geometric mean of diag(L) (rotated second-order cones).
m = size(A, 1);
if n == 2
    % x = [l11 l21 l22 d1 d2 t]
    nv = 6;  iL = [1 0; 2 3];  id = [4 5];  it = 6;
else
    % x = [l11 l21 l22 l31 l32 l33 d1 d2 d3 s1 s2 t]
    nv = 12;  iL = [1 0 0; 2 3 0; 4 5 6];  id = [7 8 9];  it = 12;
end
soc = cell(1, m + n - 1);
for i = 1:m
    Ai = zeros(n, nv);
    for r = 1:n                % component r of (A L)_i = sum_c a_c L(c, r)
        for c = r:n
            Ai(r, iL(c, r)) = A(i, c);
        end
    end
    di = zeros(nv, 1);  di(id) = -A(i, :)';
    soc{i} = secondordercone(Ai, zeros(n, 1), di, -b(i));
end
geo = @(p, q, r) deal_cone(nv, p, q, r);     % r^2 <= p q
if n == 2
    soc{m + 1} = geo(iL(1, 1), iL(2, 2), it);
else
    soc{m + 1} = geo(iL(1, 1), iL(2, 2), 10);
    soc{m + 2} = geo(iL(3, 3), it, 11);
    soc{m + 3} = geo(10, 11, it);
end
f = zeros(nv, 1);  f(it) = -1;
co = optimoptions('coneprog', 'Display', 'off');
[x, ~, flag] = coneprog(f, soc, [], [], [], [], [], [], co);
if flag ~= 1
    error('firi_polytope:mvie', 'coneprog failed on the inscribed-ellipsoid SOCP (flag %d).', flag);
end
L = zeros(n);
for r = 1:n
    for c = 1:r, L(r, c) = x(iL(r, c)); end
end
d = x(id);
vol = prod(diag(L));
end

function s = deal_cone(nv, p, q, r)
% |[2 x_r; x_p - x_q]| <= x_p + x_q   <=>   x_r^2 <= x_p x_q, x_p, x_q >= 0
Ac = zeros(2, nv);  Ac(1, r) = 2;  Ac(2, p) = 1;  Ac(2, q) = -1;
dc = zeros(nv, 1);  dc([p q]) = 1;
s = secondordercone(Ac, zeros(2, 1), dc, 0);
end

function [A, b] = vertex_hull_hrep(U)
% H-representation of conv(U) with unit rows (2-D or 3-D).
cs = ConvexSet.polygon(U);
A = cs.A;  b = cs.b;
end
