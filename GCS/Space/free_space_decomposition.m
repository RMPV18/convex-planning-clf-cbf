function freeRegions = free_space_decomposition(bbox, obstacles)
% FREE_SPACE_DECOMPOSITION
% Decomposes 2D free space (bounded workspace minus non-convex obstacles)
% into a minimal set of non-overlapping convex polygons whose union
% equals the free space.
% Returns cell array of Nx2 vertex arrays, each a convex polygon.

pbox = polyshape(bbox(:,1), bbox(:,2));
pobs = polyshape();
for k = 1:numel(obstacles)
    pk = polyshape(obstacles{k}(:,1), obstacles{k}(:,2));
    if ~isempty(pk.Vertices)
        pobs = union(pobs, pk);
    end
end

pfree = subtract(pbox, pobs);
if isempty(pfree.Vertices)
    error('Free space is empty.');
end

% Triangulate each region (handles multiply-connected polygons)
T = triangulation(pfree);
tris = T.ConnectivityList;
pts  = T.Points;
nbrList = neighbors(T); % [nTri x 3], NaN where no neighbor

nTri = size(tris,1);
polys = cell(nTri,1);
for i = 1:nTri
    polys{i} = pts(tris(i,:), :);
end

% Adjacency as cell of neighbor index sets (initialized from triangulation)
adj = cell(nTri,1);
for i = 1:nTri
    n = nbrList(i,:);
    n = n(~isnan(n));
    adj{i} = n(:)';
end

active = true(nTri,1);

changed = true;
while changed
    changed = false;
    idxList = find(active)';
    for i = idxList
        if ~active(i), continue; end
        neigh = adj{i}(active(adj{i})); %Active (1-true) triangle (i) neighboors indexes.
        for j = neigh
            merged = tryMergeConvex(polys{i}, polys{j});
            if ~isempty(merged)
                polys{i} = merged;
                active(j) = false;

                % merge adjacency lists, drop self and dead refs
                newAdj = unique([adj{i}, adj{j}]);
                newAdj(newAdj==i | newAdj==j) = [];
                adj{i} = newAdj;

                % redirect any list pointing to j -> i
                for m = find(active)'
                    adj{m}(adj{m}==j) = i;
                    adj{m} = unique(adj{m});
                    adj{m}(adj{m}==m) = [];
                end

                changed = true;
                break; % restart neighbor scan for i
            end
        end
    end
end

freeRegions = polys(active);
end


function mergedVerts = tryMergeConvex(polyA, polyB)
% Attempts to merge two polygons (sharing an edge) into a single
% convex polygon. Returns [] if not possible (share no edge, or
% union not convex, or union not simply connected).

mergedVerts = [];

pa = polyshape(polyA(:,1), polyA(:,2));
pb = polyshape(polyB(:,1), polyB(:,2));

pu = union(pa, pb);

if pu.NumRegions ~= 1 || pu.NumHoles ~= 0
    return; % not a single simple polygon -> not adjacent along an edge
end

V = pu.Vertices;
V = V(~any(isnan(V),2), :);
if size(V,1) < 3
    return;
end


if ~isConvexPolygon(V)
    return;
end

mergedVerts = V;
end


function tf = isConvexPolygon(V)
% Checks convexity of a simple polygon given ordered vertices Vx2.
n = size(V,1);

signs = zeros(n,1);
for i = 1:n
    p0 = V(i,:);
    p1 = V(mod(i,n)+1,:);
    p2 = V(mod(i+1,n)+1,:);
    e1 = p1 - p0;
    e2 = p2 - p1;
    cr = e1(1)*e2(2) - e1(2)*e2(1);
    signs(i) = sign(cr);
end

signs(abs(signs) < 1e-9) = 0;
posCount = sum(signs > 0);
negCount = sum(signs < 0);
tf = (posCount == 0) || (negCount == 0);
end


