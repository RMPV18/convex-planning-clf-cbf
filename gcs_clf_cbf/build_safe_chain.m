function chain = build_safe_chain(regionsV, seq, obstacles, mode, firiOpts)
%BUILD_SAFE_CHAIN  Offline chain of overlapping safe convex sets along a
%   region sequence: region, bridge, region, bridge, ..., region.
%
%   chain = build_safe_chain(regionsV, seq, obstacles)
%   chain = build_safe_chain(regionsV, seq, obstacles, mode, firiOpts)
%
%   regionsV   cell of m_k x 2 vertex arrays: convex free-space regions
%              (free_space_decomposition); consecutive regions of seq share
%              a facet
%   seq        1 x M region indices traversed by the planned path
%   obstacles  cell of CONVEX obstacle pieces (convex_obstacle_pieces of the
%              inflated obstacles)
%   mode       'firi' (default): bridge B_j = FIRI polytope (firi_polytope)
%              'seed'          : bridge B_j = its seed Q_j (no FIRI; used by
%                                the online version to place targets only)
%   firiOpts   options of firi_polytope
%
%   For the facet F_j = R_{seq(j)} ∩ R_{seq(j+1)} the seed is the
%   quadrilateral Q_j = conv(F_j, m - e1 n, m + e2 n) (m midpoint of F_j, n
%   unit normal, e1/e2 = 0.9 x the depth of each region along -n/+n): it is
%   convex, lies in R_{seq(j)} ∪ R_{seq(j+1)} (obstacle-free) and has interior
%   on both sides of F_j, so B_j ⊇ Q_j overlaps both regions with interior.
%   Bounding box of the FIRI inflation (Gonzalez-Garcia et al., PiP,
%   "Heuristic for Task-Informed Polytope Inflation"): the axis-aligned box
%   of the two adjacent regions, i.e. the box around their intersection F_j
%   spanned by the sets it joins (confines the inflation to the crossing).
%
%   chain.sets      1 x (2M-1) cell of struct A, b (unit rows)
%   chain.kind      1 x (2M-1) cell: 'region' | 'bridge'
%   chain.crossing  1 x (2M-1): index j of the crossing each pair (S_i,
%                   S_{i+1}) belongs to (pairs 2j-1 and 2j belong to F_j)
%   chain.facets    1 x (M-1) cell of 2 x 2 facet endpoints
%   chain.seeds     1 x (M-1) cell of seed vertices
%   chain.bridgesV  1 x (M-1) cell of bridge vertices (for plotting)
%   chain.firi      1 x (M-1) firi_polytope info ([] for 'seed')
%   chain.time      wall-clock time of the construction

if nargin < 4 || isempty(mode), mode = 'firi'; end
if nargin < 5, firiOpts = struct(); end
assert(any(strcmp(mode, {'firi', 'seed'})), 'build_safe_chain:mode', 'mode must be ''firi'' or ''seed''.');
t0 = tic;
M = numel(seq);
regH = cell(1, M);
for k = 1:M
    cs = ConvexSet.polygon(regionsV{seq(k)});
    regH{k} = struct('A', cs.A, 'b', cs.b);
end
chain = struct('sets', {{}}, 'kind', {{}}, 'crossing', [], 'facets', {cell(1, M-1)}, ...
    'seeds', {cell(1, M-1)}, 'bridgesV', {cell(1, M-1)}, 'firi', {cell(1, M-1)}, 'time', 0);
chain.sets{1} = regH{1};  chain.kind{1} = 'region';
for j = 1:M-1
    Va = regionsV{seq(j)};  Vb = regionsV{seq(j+1)};
    F = shared_facet(Va, Vb, regH{j}, regH{j+1});
    m = mean(F, 1)';  tdir = F(2, :)' - F(1, :)';
    nrm = [-tdir(2); tdir(1)] / norm(tdir);
    if any(regH{j+1}.A * (m + 1e-6 * nrm) - regH{j+1}.b > 1e-12), nrm = -nrm; end   % n points a -> b
    e1 = 0.9 * depth(regH{j}, m, -nrm);  e2 = 0.9 * depth(regH{j+1}, m, nrm);
    Q = [F(1, :); (m - e1 * nrm)'; F(2, :); (m + e2 * nrm)'];
    chain.facets{j} = F;  chain.seeds{j} = Q;
    if strcmp(mode, 'firi')
        V = [Va; Vb];
        boxA = [eye(2); -eye(2)];  boxb = [max(V, [], 1)'; -min(V, [], 1)'];
        [A, b, info] = firi_polytope(Q, obstacles, boxA, boxb, firiOpts);
        chain.firi{j} = info;
    else
        cs = ConvexSet.polygon(Q);  A = cs.A;  b = cs.b;
    end
    chain.bridgesV{j} = polytope_vertices_2d(A, b);
    chain.sets{end+1} = struct('A', A, 'b', b);  chain.kind{end+1} = 'bridge'; %#ok<AGROW>
    chain.sets{end+1} = regH{j+1};               chain.kind{end+1} = 'region'; %#ok<AGROW>
end
chain.crossing = ceil((1:numel(chain.sets)) / 2);
chain.time = toc(t0);
end

% =========================================================================
function F = shared_facet(Va, Vb, Ha, Hb)
% Endpoints of the common facet: vertices of either polygon lying in both.
tol = 1e-7;
P = [Va; Vb];
in = all(P * Ha.A' - Ha.b' <= tol, 2) & all(P * Hb.A' - Hb.b' <= tol, 2);
P = uniquetol(P(in, :), 1e-9, 'ByRows', true);
if size(P, 1) < 2
    error('build_safe_chain:facet', 'Consecutive regions do not share a facet.');
end
d = P(end, :) - P(1, :);
[~, i1] = min(P * d');  [~, i2] = max(P * d');
F = P([i1 i2], :);
if norm(F(2, :) - F(1, :)) < 1e-9
    error('build_safe_chain:facet', 'Consecutive regions touch at a single point.');
end
end

function t = depth(H, p, d)
% Largest t with p + t d in {A x <= b}.
ad = H.A * d;  sl = H.b - H.A * p;
t = min(sl(ad > 1e-12) ./ ad(ad > 1e-12));
end
