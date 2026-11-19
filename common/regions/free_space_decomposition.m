function regions = free_space_decomposition(bbox, obstacles)
%FREE_SPACE_DECOMPOSITION  Convex decomposition of 2-D free space.
%
%   regions = free_space_decomposition(bbox, obstacles)
%
%   bbox       m x 2 vertices of the workspace polygon
%   obstacles  cell of m_k x 2 vertex arrays (possibly non-convex)
%   regions    cell of m_k x 2 vertex arrays: non-overlapping convex
%              polygons whose union is bbox minus the obstacles
%
%   Method (Hertel-Mehlhorn style): constrained triangulation of the free
%   space, then greedy merging of edge-adjacent pieces while the union
%   stays convex. The count is small but not guaranteed minimal.

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
    error('free_space_decomposition:empty', 'Free space is empty.');
end

% Triangulation (handles multiply-connected polygons)
T = triangulation(pfree);
tris = T.ConnectivityList;
pts  = T.Points;
triNeighbors = neighbors(T);              % nTri x 3, NaN where no neighbor

nTri = size(tris,1);
polys = cell(nTri,1);
adj = cell(nTri,1);                       % neighbor index sets, updated while merging
for i = 1:nTri
    polys{i} = pts(tris(i,:), :);
    nb = triNeighbors(i,:);
    adj{i} = nb(~isnan(nb));
end

active = true(nTri,1);
changed = true;
while changed
    changed = false;
    for i = find(active)'
        if ~active(i), continue; end
        for j = adj{i}(active(adj{i}))    % active neighbors of piece i
            merged = try_merge_convex(polys{i}, polys{j});
            if ~isempty(merged)
                polys{i} = merged;
                active(j) = false;

                % merge adjacency lists, drop self and dead refs
                newAdj = unique([adj{i}, adj{j}]);
                newAdj(newAdj==i | newAdj==j) = [];
                adj{i} = newAdj;

                % redirect any list pointing to j -> i
                for k = find(active)'
                    adj{k}(adj{k}==j) = i;
                    adj{k} = unique(adj{k});
                    adj{k}(adj{k}==k) = [];
                end

                changed = true;
                break;                    % restart neighbor scan for i
            end
        end
    end
end

regions = polys(active);
end


function mergedVerts = try_merge_convex(polyA, polyB)
% Union of two polygons if it is a single convex polygon without holes,
% otherwise [].
mergedVerts = [];
pu = union(polyshape(polyA(:,1), polyA(:,2)), polyshape(polyB(:,1), polyB(:,2)));
if pu.NumRegions ~= 1 || pu.NumHoles ~= 0
    return;                               % not adjacent along an edge
end
V = pu.Vertices;
V = V(~any(isnan(V),2), :);
if size(V,1) < 3 || ~is_convex_polygon(V)
    return;
end
mergedVerts = V;
end


function tf = is_convex_polygon(V)
% Convexity of a simple polygon with ordered vertices V (m x 2).
m = size(V,1);
signs = zeros(m,1);
for k = 1:m
    e1 = V(mod(k,m)+1,:) - V(k,:);
    e2 = V(mod(k+1,m)+1,:) - V(mod(k,m)+1,:);
    signs(k) = sign(e1(1)*e2(2) - e1(2)*e2(1));
end
tf = all(signs >= 0) || all(signs <= 0);
end
