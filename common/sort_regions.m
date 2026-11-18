function regions = sort_regions(regions, ref, criteria)
%SORT_REGIONS  Reorder vertex-array regions by distance to a reference point.
%
%   regions = sort_regions(regions, ref)                 % by centroid distance
%   regions = sort_regions(regions, ref, 'closervertex') % by nearest-vertex
%                                                        % distance (ties keep
%                                                        % centroid order)
%   regions  cell of m_k x n vertex arrays;  ref  n x 1 point.

K = numel(regions);
dCentroid = zeros(1, K);
for k = 1:K, dCentroid(k) = vecnorm(mean(regions{k}, 1) - ref'); end
[~, order] = sort(dCentroid);
if nargin < 3
    regions = regions(order);
elseif strcmp(criteria, 'closervertex')
    dVertex = zeros(1, K);
    for k = 1:K, dVertex(k) = min(vecnorm(regions{k} - ref', 2, 2)); end
    [~, order2] = sort(dVertex(order));
    regions = regions(order(order2));
end
end
