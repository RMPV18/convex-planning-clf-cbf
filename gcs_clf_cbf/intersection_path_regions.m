function seq = intersection_path_regions(path, regions_cs, q0)
%INTERSECTION_PATH_REGIONS  Region sequence of an intersection-graph path.
%
%   seq = intersection_path_regions(path, regions_cs, q0)
%
%   path        {'s', 'i-j', 'j-k', ..., 't'} from ShortestPathProblem on the
%               graph of build_intersection_graph
%   regions_cs  cell of ConvexSet (the regions indexed in the vertex names)
%   q0          start point (resolves the first region when the path has a
%               single pair vertex)
%   seq         1 x (m+1) region indices traversed; the k-th pair vertex is
%               the facet between seq(k) and seq(k+1)

pairs = path(2:end-1);
e = cell2mat(cellfun(@(s) sscanf(s, '%d-%d')', pairs(:), 'UniformOutput', false));
m = size(e, 1);
assert(m >= 1, 'intersection_path_regions:path', 'The path has no pair vertex.');
seq = zeros(1, m + 1);
if m == 1
    if regions_cs{e(1)}.contains(q0, 1e-9), seq = e; else, seq = e([2 1]); end
    return
end
seq(1) = setdiff(e(1, :), e(2, :));
for k = 1:m-1
    seq(k+1) = intersect(e(k, :), e(k+1, :));
end
seq(m+1) = setdiff(e(m, :), e(m-1, :));
end
