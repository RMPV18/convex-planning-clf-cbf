function pieces = convex_obstacle_pieces(obstacles)
%CONVEX_OBSTACLE_PIECES  Split 2-D polygonal obstacles into convex pieces.
%
%   pieces = convex_obstacle_pieces(obstacles)
%
%   obstacles  cell of m_k x 2 vertex arrays (possibly non-convex, possibly
%              overlapping)
%   pieces     cell of 3 x 2 triangles whose union equals the union of the
%              obstacles (constrained triangulation of the polyshape union),
%              as needed by methods that require convex obstacles (FIRI).

pu = polyshape();
for k = 1:numel(obstacles)
    pu = union(pu, polyshape(obstacles{k}));
end
T = triangulation(pu);
pieces = cell(1, size(T.ConnectivityList, 1));
for k = 1:numel(pieces)
    pieces{k} = T.Points(T.ConnectivityList(k, :), :);
end
end
