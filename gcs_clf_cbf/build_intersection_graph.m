function g = build_intersection_graph(q0, qT, regions_cs, neighbors)
%BUILD_INTERSECTION_GRAPH  Graph of convex sets whose vertices are the
%   pairwise intersections of adjacent free-space regions.
%
%   g = build_intersection_graph(q0, qT, regions_cs, neighbors)
%
%   q0, qT      n x 1 start and goal points (vertices 's' and 't').
%   regions_cs  1 x K cell of ConvexSet, the convex free-space regions.
%   neighbors   1 x K cell, neighbors{i} = indices of regions adjacent to i
%               (from adjacent_regions).
%
%   Vertices: 's', 't' and one vertex 'i-j' (i < j) per adjacent pair,
%   holding the intersection Q_i ∩ Q_j. A shortest path through these
%   vertices is a sequence of crossing points between consecutive regions.
%   Edges: 'i-j' -> 'j-k' for every region j shared by the two pairs
%   (k ~= i); 's' -> every pair of a region containing q0; every pair of a
%   region containing qT -> 't'.

g = GraphOfConvexSets();
g.add_vertex('s', ConvexSet.point(q0));
g.add_vertex('t', ConvexSet.point(qT));

for i = 1:numel(regions_cs)
    for j = neighbors{i}
        if i < j
            g.add_vertex(pair_name(i, j), ConvexSet([regions_cs{i}.A; regions_cs{j}.A], ...
                                                    [regions_cs{i}.b; regions_cs{j}.b]));
        end
    end
end

for i = 1:numel(regions_cs)
    for j = neighbors{i}
        for k = neighbors{j}(neighbors{j} ~= i)
            g.add_edge(pair_name(i, j), pair_name(j, k));
        end
    end
end

for k = 1:numel(regions_cs)
    if regions_cs{k}.contains(q0, 1e-9)
        for j = neighbors{k}, g.add_edge('s', pair_name(k, j)); end
    end
    if regions_cs{k}.contains(qT, 1e-9)
        for j = neighbors{k}, g.add_edge(pair_name(k, j), 't'); end
    end
end
end

function name = pair_name(i, j)
name = [num2str(min(i, j)) '-' num2str(max(i, j))];
end
