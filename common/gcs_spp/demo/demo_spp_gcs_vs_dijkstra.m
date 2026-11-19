%% DEMO_SPP_GCS_VS_DIJKSTRA  Sanity check of ShortestPathProblem.
% When every vertex set is a single point, SPP in GCS reduces to a classic
% shortest-path problem with Euclidean edge weights. This script solves it
% with ShortestPathProblem (relaxation -> rounding -> restriction) and
% asserts that the cost matches MATLAB's shortestpath (Dijkstra).
% Requires: YALMIP + a conic solver.

clear; clc;

names = {'s','v1','v2','v3','v4','t'};
coords = [0 0; 2 3; 2 -3; 5 2; 5 -2; 8 0]';   % 2 x 6, one column per vertex

g = GraphOfConvexSets();
for k = 1:numel(names)
    g.add_vertex(names{k}, ConvexSet.point(coords(:,k)));
end

edge_list = {
    's','v1'; 's','v2'; ...
    'v1','v3'; 'v1','v4'; 'v2','v3'; 'v2','v4'; ...
    'v3','t'; 'v4','t'
};
for e = 1:size(edge_list,1)
    g.add_edge(edge_list{e,1}, edge_list{e,2});
end

spp = ShortestPathProblem(g, 's', 't');
[path, ~, cost, info] = spp.solve(false);

fprintf('GCS solver path: %s (cost = %.6f)\n', strjoin(path, ' -> '), cost);

%% Classic Dijkstra cross-check
nE = size(edge_list,1);
src = zeros(nE,1);  dst = zeros(nE,1);  w = zeros(nE,1);
idx = containers.Map(names, 1:numel(names));
for e = 1:nE
    src(e) = idx(edge_list{e,1});
    dst(e) = idx(edge_list{e,2});
    w(e) = norm(coords(:,src(e)) - coords(:,dst(e)));
end

G = digraph(src, dst, w, names);
[dijkstra_path, dijkstra_cost] = shortestpath(G, 's', 't');

fprintf('Dijkstra path:   %s (cost = %.6f)\n', strjoin(dijkstra_path, ' -> '), dijkstra_cost);

assert(abs(cost - dijkstra_cost) < 1e-4, ...
    'Mismatch between GCS relaxation+restriction and classic Dijkstra.');
disp('Sanity check passed: costs match.');
