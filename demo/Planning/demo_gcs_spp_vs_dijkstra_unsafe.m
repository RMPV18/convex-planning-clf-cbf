%% demo_sanity_check.m
% When every vertex set is a single point, the SPP in GCS collapses to
% an ordinary shortest-path problem with fixed edge weights (Euclidean
% distances between the points). This script builds such a graph, runs
% ShortestPathProblem, and cross-checks the result against MATLAB's
% built-in graph shortestpath (classic Dijkstra), which should agree
% up to solver tolerance.

clear; clc;

names = {'s','v1','v2','v3','v4','t'};
coords = [0 0; 2 3; 2 -3; 5 2; 5 -2; 8 0]'; % 2 x 6, one column per vertex

g = GraphOfConvexSets();
for i = 1:numel(names)
    g.add_vertex(names{i}, ConvexSet.point(coords(:,i)));
end

edge_list = {
    's','v1'; 's','v2'; ...
    'v1','v3'; 'v1','v4'; 'v2','v3'; 'v2','v4'; ...
    'v3','t'; 'v4','t'
};
for k = 1:size(edge_list,1)
    g.add_edge(edge_list{k,1}, edge_list{k,2});
end

spp = ShortestPathProblem(g, 's', 't');
[path, ~, cost, info] = spp.solve(false);

fprintf('GCS solver path: %s (cost = %.6f)\n', strjoin(path, ' -> '), cost);

%% Classic Dijkstra cross-check
src = zeros(size(edge_list,1),1);
dst = zeros(size(edge_list,1),1);
w   = zeros(size(edge_list,1),1);
idx = containers.Map(names, 1:numel(names));
for k = 1:size(edge_list,1)
    i = idx(edge_list{k,1});
    j = idx(edge_list{k,2});
    src(k) = i; dst(k) = j;
    w(k) = norm(coords(:,i) - coords(:,j));
end

G = digraph(src, dst, w, names);
[dijkstra_path, dijkstra_cost] = shortestpath(G, 's', 't');

fprintf('Dijkstra path:   %s (cost = %.6f)\n', strjoin(dijkstra_path, ' -> '), dijkstra_cost);

assert(abs(cost - dijkstra_cost) < 1e-4, ...
    'Mismatch between GCS relaxation+restriction and classic Dijkstra.');
disp('Sanity check passed: costs match.');
