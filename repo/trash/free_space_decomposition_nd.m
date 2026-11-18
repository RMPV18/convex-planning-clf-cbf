function free_space_decomposition_nd()
% N-DIMENSIONAL FREE SPACE DECOMPOSITION
% Decomposes an axis-aligned bounding box minus a set of non-convex
% obstacles (each supplied as a union of convex polytopes) into a
% minimal-effort set of non-overlapping convex axis-aligned boxes whose
% union approximates the free space. Works for any dimension d (2D, 3D,
% 4D, ...). Demo below runs in 3D.

warning('on','all'); % warnings NOT suppressed

%% ---------------- Problem definition (3D demo) ----------------
bboxLow  = [0 0 0];
bboxHigh = [10 10 10];

obstacleParts = {};
obstacleParts{end+1} = boxVertices([2 2 2],[4 4 5]);          % box obstacle
obstacleParts{end+1} = boxVertices([6 6 0],[8 8 3]);          % L-shape part A
obstacleParts{end+1} = boxVertices([6 6 3],[7 7 6]);          % L-shape part B
obstacleParts{end+1} = [1 7 1; 3 7 1; 2 9 1; 2 8 4];          % tetrahedron
obstacleParts{end+1} = boxVertices([6 1 6],[9 3 9]);          % box obstacle

obstacles = cell(size(obstacleParts));
for k = 1:numel(obstacleParts)
    [A,b] = vert2conND(obstacleParts{k});
    obstacles{k} = struct('A',A,'b',b,'V',obstacleParts{k});
end

%% ---------------- Decomposition ----------------
maxDepth    = 6;
minCellSize = 0.15;

boxesLeaf = decomposeFreeSpaceND(bboxLow, bboxHigh, obstacles, ...
    maxDepth, minCellSize, 0, {});
fprintf('Leaf boxes before merge: %d\n', numel(boxesLeaf));

boxesMerged = mergeBoxesND(boxesLeaf);
fprintf('Boxes after merge: %d\n', numel(boxesMerged));

%% ---------------- Verification ----------------
verifyNoOverlap(boxesMerged);
checkCoverageMC(boxesMerged, bboxLow, bboxHigh, obstacles, 200000);

%% ---------------- Plot (3D) ----------------
plotResult3D(boxesMerged, obstacleParts, bboxLow, bboxHigh);

end


%% ========================= CORE ALGORITHM =========================

function boxesOut = decomposeFreeSpaceND(low, high, obstacles, maxDepth, minCellSize, depth, boxesOut)
d = numel(low);

overlapping = {};
for k = 1:numel(obstacles)
    if boxIntersectsPolytope(low, high, obstacles{k})
        overlapping{end+1} = obstacles{k}; %#ok<AGROW>
    end
end

if isempty(overlapping)
    boxesOut{end+1} = struct('low',low,'high',high);
    return;
end

for k = 1:numel(overlapping)
    if boxFullyInsidePolytope(low, high, overlapping{k})
        return; % fully occupied, discard
    end
end

extent = high - low;
if depth >= maxDepth || any(extent <= minCellSize)
    center = (low + high) / 2;
    if ~pointInsideAnyObstacle(center, overlapping)
        warning('Approximating boundary-adjacent cell at resolution limit (depth=%d, extent=%s).', ...
            depth, mat2str(extent,3));
        boxesOut{end+1} = struct('low',low,'high',high);
    end
    return;
end

mid = (low + high) / 2;
nChildren = 2^d;
for mask = 0:nChildren-1
    bits = bitget(mask, 1:d);
    childLow = low; childHigh = high;
    childLow(bits==1)  = mid(bits==1);
    childHigh(bits==0) = mid(bits==0);
    boxesOut = decomposeFreeSpaceND(childLow, childHigh, overlapping, ...
        maxDepth, minCellSize, depth+1, boxesOut);
end
end


function boxes = mergeBoxesND(boxes)
changed = true;
while changed
    changed = false;
    n = numel(boxes);
    for i = 1:n
        if isempty(boxes{i}), continue; end
        for j = i+1:n
            if isempty(boxes{j}), continue; end
            merged = tryMergeBoxPair(boxes{i}, boxes{j});
            if ~isempty(merged)
                boxes{i} = merged;
                boxes{j} = [];
                changed = true;
            end
        end
    end
    boxes = boxes(~cellfun(@isempty, boxes));
end
end


function merged = tryMergeBoxPair(b1, b2)
merged = [];
diffLow  = abs(b1.low  - b2.low)  > 1e-9;
diffHigh = abs(b1.high - b2.high) > 1e-9;
diffDims = find(diffLow | diffHigh);
if numel(diffDims) ~= 1
    return;
end
k = diffDims;
if abs(b1.high(k) - b2.low(k)) < 1e-9
    newLow = b1.low; newHigh = b1.high; newHigh(k) = b2.high(k);
elseif abs(b2.high(k) - b1.low(k)) < 1e-9
    newLow = b2.low; newHigh = b2.high; newHigh(k) = b1.high(k);
else
    return;
end
merged = struct('low', newLow, 'high', newHigh);
end


%% ===================== GEOMETRIC PRIMITIVES =====================

function [A,b] = vert2conND(V)
% Vertex-to-halfspace (H-representation) conversion for a convex
% polytope in arbitrary dimension d, via facet hyperplane fitting.
d  = size(V,2);
K  = convhulln(V);
nF = size(K,1);
A  = zeros(nF,d);
b  = zeros(nF,1);
c0 = mean(V,1);
for i = 1:nF
    P = V(K(i,:),:);
    M = [P, -ones(d,1)];
    [~,~,Vt] = svd(M);
    coeffs = Vt(:,end);
    a = coeffs(1:d)';
    c = coeffs(d+1);
    if a*c0' > c
        a = -a; c = -c;
    end
    nrm = norm(a);
    if nrm < 1e-12
        warning('Degenerate facet %d in vert2conND; constraint skipped (non-restrictive).', i);
        continue;
    end
    A(i,:) = a / nrm;
    b(i)   = c / nrm;
end
end


function C = boxCorners(low, high)
d = numel(low);
n = 2^d;
C = zeros(n,d);
for i = 0:n-1
    bits = bitget(i,1:d);
    pt = low;
    pt(bits==1) = high(bits==1);
    C(i+1,:) = pt;
end
end


function V = boxVertices(low, high)
V = boxCorners(low, high);
end


function tf = boxIntersectsPolytope(low, high, obs)
d = numel(low);
persistent hasOptim
if isempty(hasOptim)
    hasOptim = ~isempty(which('linprog'));
    if ~hasOptim
        warning('Optimization Toolbox (linprog) not found: using conservative SAT fallback for box-polytope intersection tests. This may over-subdivide but will not leak obstacle volume into free space.');
    end
end
if hasOptim
    opts = optimoptions('linprog','Display','none');
    [~,~,exitflag] = linprog(zeros(d,1), obs.A, obs.b, [], [], low(:), high(:), opts);
    tf = (exitflag == 1);
else
    tf = boxIntersectsPolytopeSAT(low, high, obs);
end
end


function tf = boxIntersectsPolytopeSAT(low, high, obs)
% Sufficient (not exhaustive) separating-axis test: box axes + polytope
% face normals. Never reports false separation (safe), may fail to
% detect some true separations (oblique edge-edge cases), which only
% costs efficiency, not correctness.
d = numel(low);
for dim = 1:d
    if max(obs.V(:,dim)) < low(dim) - 1e-9 || min(obs.V(:,dim)) > high(dim) + 1e-9
        tf = false; return;
    end
end
corners = boxCorners(low, high);
nF = size(obs.A,1);
for i = 1:nF
    proj = corners * obs.A(i,:)';
    if min(proj) > obs.b(i) + 1e-9
        tf = false; return;
    end
end
tf = true;
end


function tf = boxFullyInsidePolytope(low, high, obs)
corners = boxCorners(low, high);
tf = all(all(corners*obs.A' <= obs.b' + 1e-9));
end


function tf = pointInsideObstacle(pt, obs)
tf = all(obs.A*pt(:) <= obs.b + 1e-9);
end


function tf = pointInsideAnyObstacle(pt, obstacles)
tf = false;
for k = 1:numel(obstacles)
    if pointInsideObstacle(pt, obstacles{k})
        tf = true; return;
    end
end
end


%% ========================= VERIFICATION =========================

function verifyNoOverlap(boxes)
n = numel(boxes);
for i = 1:n
    for j = i+1:n
        if boxesOverlap(boxes{i}, boxes{j})
            warning('Boxes %d and %d overlap.', i, j);
        end
    end
end
end


function tf = boxesOverlap(b1,b2)
tf = all(b1.low < b2.high - 1e-9) && all(b2.low < b1.high - 1e-9);
end


function checkCoverageMC(boxes, bboxLow, bboxHigh, obstacles, nSamples)
d = numel(bboxLow);
pts = bboxLow + rand(nSamples,d).*(bboxHigh-bboxLow);

trueFree = true(nSamples,1);
for k = 1:numel(obstacles)
    inside = all(pts*obstacles{k}.A' <= obstacles{k}.b' + 1e-9, 2);
    trueFree(inside) = false;
end

decompFree = false(nSamples,1);
for i = 1:numel(boxes)
    inBox = all(pts >= boxes{i}.low & pts <= boxes{i}.high, 2);
    decompFree = decompFree | inBox;
end

occupiedLeak   = sum(~trueFree & decompFree);
freeUncovered  = sum(trueFree & ~decompFree);
bboxVol        = prod(bboxHigh-bboxLow);
estFreeVol     = mean(trueFree)*bboxVol;
decompVol      = sum(cellfun(@(bx) prod(bx.high-bx.low), boxes));

fprintf('Monte Carlo (%d samples): estimated free volume=%.4f, decomposed volume=%.4f, bbox volume=%.4f\n', ...
    nSamples, estFreeVol, decompVol, bboxVol);

if occupiedLeak > 0
    warning('%d/%d sampled obstacle points fall inside decomposed free-space boxes (unsafe leakage).', occupiedLeak, nSamples);
end
if freeUncovered > 0
    warning('%d/%d sampled free points not covered by decomposition (boundary approximation gap; controlled by minCellSize/maxDepth).', freeUncovered, nSamples);
end
end


%% =========================== PLOTTING ============================

function plotResult3D(boxes, obstacleParts, bboxLow, bboxHigh)
figure; hold on; axis equal; grid on; view(3);
cmap = lines(max(numel(boxes),1));
for i = 1:numel(boxes)
    plotBox3D(boxes{i}.low, boxes{i}.high, cmap(i,:), 0.35);
end
for k = 1:numel(obstacleParts)
    V = obstacleParts{k};
    K = convhulln(V);
    patch('Vertices',V,'Faces',K,'FaceColor','k','FaceAlpha',0.9,'EdgeColor','none');
end
xlim([bboxLow(1) bboxHigh(1)]); ylim([bboxLow(2) bboxHigh(2)]); zlim([bboxLow(3) bboxHigh(3)]);
title(sprintf('3D free space: %d convex boxes', numel(boxes)));
xlabel('x'); ylabel('y'); zlabel('z');
end


function plotBox3D(low, high, color, alpha)
C = boxCorners(low, high);
faces = [1 2 4 3; 5 6 8 7; 1 2 6 5; 3 4 8 7; 1 3 7 5; 2 4 8 6];
patch('Vertices', C, 'Faces', faces, 'FaceColor', color, 'FaceAlpha', alpha, 'EdgeColor','k');
end