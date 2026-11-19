function inflated = inflate_obstacles(obstacles, margin)
%INFLATE_OBSTACLES  Grow 2-D polygonal obstacles by a clearance margin.
%
%   inflated = inflate_obstacles(obstacles, margin)
%
%   obstacles  cell of m_k x 2 vertex arrays (possibly non-convex)
%   margin     clearance >= 0
%   inflated   cell of vertex arrays; each contains the Minkowski sum of the
%              obstacle with a disc of radius margin (mitred corners), so any
%              point outside all inflated obstacles is >= margin from the
%              true ones. Overlapping inflated obstacles are left to the
%              union taken by free_space_decomposition.

assert(isscalar(margin) && margin >= 0, 'inflate_obstacles:margin', 'margin must be a scalar >= 0.');
inflated = obstacles;
if margin == 0, return; end
for k = 1:numel(obstacles)
    pb = polybuffer(polyshape(obstacles{k}), margin, 'JointType', 'miter');
    inflated{k} = pb.Vertices;
end
end
