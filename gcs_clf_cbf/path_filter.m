function [new_path, new_wp, rm_i, keep] = path_filter(obstacles_cs, path, wp)
%PATH_FILTER  Drop waypoints that do not touch any obstacle.
%
%   [new_path, new_wp, rm_i, keep] = path_filter(obstacles_cs, path, wp)
%
%   obstacles_cs  cell of ConvexSet, convex pieces of the obstacles
%   path          1 x M cell of vertex names (from ShortestPathProblem)
%   wp            n x M waypoints (positions of the path vertices)
%
%   An interior waypoint is kept iff it lies in some obstacle piece within
%   1e-3, i.e. it is a contact point (corner) that bends the path; the
%   first and last waypoints are always kept. A waypoint in the open free
%   space lies on a straight segment and is removed.
%   keep  1 x M logical mask;  rm_i  indices of removed waypoints.

M = size(wp, 2);
keep = true(1, M);
for k = 2:M-1
    keep(k) = any(cellfun(@(cs) cs.contains(wp(:,k), 1e-3), obstacles_cs));
end
rm_i = find(~keep);
new_wp = wp(:, keep);
new_path = path(keep);
end
