function c = path_clearance(W, obstacles)
%PATH_CLEARANCE  Distance from a 2-D polyline to the union of polygons.
%
%   c = path_clearance(W, obstacles)
%
%   W          2 x K polyline vertices (columns)
%   obstacles  cell of m_k x 2 vertex arrays (possibly non-convex)
%   c          minimum distance; -1 if the polyline passes through an
%              obstacle interior. Exact: for a polyline disjoint from the
%              polygons the minimum is attained at a polyline vertex or an
%              obstacle vertex.
pobs = polyshape();
for k = 1:numel(obstacles), pobs = union(pobs, polyshape(obstacles{k})); end
P = W';
inside = intersect(pobs, P);                 % parts of the polyline inside pobs
if ~isempty(inside) && sum(vecnorm(diff(inside(~any(isnan(inside), 2), :)), 2, 2)) > 1e-9
    c = -1;  return
end
[eA, eB] = boundary_edges(pobs);
c = Inf;
for k = 1:size(P, 1) - 1                     % obstacle vertices to path segments
    c = min(c, min(seg_dist(eA, P(k,:), P(k+1,:))));
end
for i = 1:size(eA, 1)                        % path vertices to obstacle edges
    c = min(c, min(seg_dist(P, eA(i,:), eB(i,:))));
end
end

function d = seg_dist(X, a, b)
% Distances from the rows of X to the segment [a, b].
ab = b - a;  t = max(0, min(1, ((X - a) * ab') / max(ab * ab', eps)));
d = vecnorm(X - (a + t .* ab), 2, 2);
end

function [A, B] = boundary_edges(p)
% Edges of every boundary loop of polyshape p.
A = [];  B = [];
for r = 1:p.NumRegions + p.NumHoles
    [x, y] = boundary(p, r);
    A = [A; x(1:end-1) y(1:end-1)]; %#ok<AGROW>
    B = [B; x(2:end) y(2:end)];     %#ok<AGROW>
end
end
