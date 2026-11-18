classdef ShortestPathProblem < handle
    % SHORTESTPATHPROBLEM  Shortest Path Problem in a Graph of Convex
    % Sets (SPP in GCS), adapted from:
    %
    %   T. Marcucci, J. Umenberger, P. A. Parrilo, R. Tedrake,
    %   "Shortest Paths in Graphs of Convex Sets," SIAM J. Optim., 2024.
    %   https://github.com/TobiaMarcucci/shortest-paths-in-graphs-of-convex-sets
    %
    % ALGORITHM (mirrors the structure of the Python repo's
    % spp.shortest_path.ShortestPathProblem, using YALMIP in place of
    % Drake and any conic solver, e.g. Mosek/Gurobi/SCS, in place of
    % Mosek):
    %
    %   1. solve_relaxation(): build and solve the convex (SOCP)
    %      relaxation of the perspective MICP. Edge-selection variables
    %      y_e are relaxed to [0,1]; each edge carries "local" endpoint
    %      variables that are the homogenization (perspective) of the
    %      true vertex position, so the vertex-set and edge-length
    %      constraints stay convex for fractional y_e.
    %   2. round_solution(): recover a single source-target path from
    %      the (possibly fractional) relaxed y_e via a greedy,
    %      y-weighted depth-first search — analogous to the repo's
    %      rounding heuristics.
    %   3. solve_restriction(): with the discrete path fixed, solve the
    %      much smaller convex program that only optimizes the
    %      continuous vertex positions along that path. This gives a
    %      feasible upper bound; combined with the relaxation's optimal
    %      value it also gives an optimality gap.
    %
    % Requires YALMIP (https://yalmip.github.io) and a conic solver on
    % the MATLAB path (Mosek, Gurobi, SCS, SDPT3, ECOS, ...).

    properties
        graph GraphOfConvexSets
        source char
        target char
    end

    methods
        function obj = ShortestPathProblem(graph, source, target)
            if ~isKey(graph.vertices, source) || ~isKey(graph.vertices, target)
                error('ShortestPathProblem:novertex', ...
                    'Source and target must be vertices of the graph.');
            end
            obj.graph = graph;
            obj.source = source;
            obj.target = target;
        end

        function [path, positions, cost, info] = solve(obj, verbose)
            % Full pipeline: relax -> round -> restrict.
            if nargin < 2
                verbose = false;
            end
            [y_val, relax_cost] = obj.solve_relaxation(verbose);
            path = obj.round_solution(y_val); %greedy
            [positions, cost] = obj.solve_restriction(path, verbose);

            info = struct();
            info.relaxation_cost = relax_cost;   % lower bound
            info.restriction_cost = cost;        % feasible upper bound
            info.optimality_gap = cost - relax_cost;
            info.edge_weights = y_val;
        end

        function [y_val, cost] = solve_relaxation(obj, verbose)
            if nargin < 2
                verbose = false;
            end
            g = obj.graph;
            nE = g.num_edges();
            dim = g.vertices(g.vertex_names{1}).dim;

            y  = sdpvar(nE, 1, 'full');           % relaxed edge indicators, in [0,1]
            pu = sdpvar(dim, nE, 'full');         % perspective point of each edge's tail
            pv = sdpvar(dim, nE, 'full');         % perspective point of each edge's head
            l  = sdpvar(nE, 1, 'full');           % edge length slack

            Constraints = [0 <= y <= 1];
            Objective = 0;

            for e = 1:nE
                tail = g.edges{e,1};
                head = g.edges{e,2};
                Xu = g.vertices(tail);
                Xv = g.vertices(head);

                % perspective (homogenized) membership constraints
                Constraints = [Constraints, Xu.perspective_constraints(pu(:,e), y(e))]; %#ok<AGROW>
                Constraints = [Constraints, Xv.perspective_constraints(pv(:,e), y(e))]; %#ok<AGROW>

                % Euclidean edge length is positively homogeneous of
                % degree 1, so its perspective is just the norm of the
                % *already-scaled* points: y_e * ||x_v - x_u|| ==
                % ||pv(:,e) - pu(:,e)||. This second-order cone
                % constraint bounds it above by the slack l(e).
                Constraints = [Constraints, cone(pv(:,e) - pu(:,e), l(e))]; %#ok<AGROW>

                Objective = Objective + l(e);
            end

            % Degree / flow-conservation constraints (one unit of flow
            % from source to target; at most one path through each
            % intermediate vertex).
            for i = 1:g.num_vertices()
                v = g.vertex_names{i};
                out_idx = g.outgoing_edges(v);
                in_idx  = g.incoming_edges(v);
                out_sum = sum_or_zero(y, out_idx);
                in_sum  = sum_or_zero(y, in_idx);

                if strcmp(v, obj.source)
                    Constraints = [Constraints, out_sum - in_sum == 1]; %#ok<AGROW>
                elseif strcmp(v, obj.target)
                    Constraints = [Constraints, out_sum - in_sum == -1]; %#ok<AGROW>
                else
                    Constraints = [Constraints, out_sum - in_sum == 0]; %#ok<AGROW>
                    Constraints = [Constraints, in_sum <= 1]; %#ok<AGROW>
                end
            end

            % Spatial ("point") conservation: the position used for v
            % when arriving via one edge must match the position used
            % when leaving via the next edge. This is what stitches the
            % per-edge perspective points into a single continuous
            % trajectory through the chosen vertices.
            for i = 1:g.num_vertices()
                v = g.vertex_names{i};
                if strcmp(v, obj.source) || strcmp(v, obj.target)
                    continue
                end
                out_idx = g.outgoing_edges(v);
                in_idx  = g.incoming_edges(v);
                lhs = sum_cols_or_zero(pu, out_idx, dim);
                rhs = sum_cols_or_zero(pv, in_idx, dim);
                Constraints = [Constraints, lhs == rhs]; %#ok<AGROW>
            end

            ops = sdpsettings('verbose', double(verbose));
            diagnostics = optimize(Constraints, Objective, ops);
            if diagnostics.problem ~= 0
                warning('ShortestPathProblem:relaxation', ...
                    'Solver reported an issue solving the relaxation: %s', ...
                    diagnostics.info);
            end

            y_val = value(y);
            cost  = value(Objective);
        end

        function path = round_solution(obj, y_val)
            % Greedy, y-weighted depth-first search with backtracking,
            % analogous to the repo's rounding step: at each vertex,
            % try outgoing edges in decreasing order of relaxed y_e.
            g = obj.graph;
            visited = containers.Map('KeyType', 'char', 'ValueType', 'logical');
            path = obj.dfs_round(g, obj.source, y_val, visited, {});
            if isempty(path)
                error('ShortestPathProblem:rounding', ...
                    'Failed to recover a source-target path from the relaxation.');
            end
        end

        function [positions, cost] = solve_restriction(obj, path, verbose)
            % Convex restriction: path (sequence of vertex names) is
            % fixed, only the continuous positions are optimized.
            if nargin < 3
                verbose = false;
            end
            g = obj.graph;
            n = numel(path);
            dim = g.vertices(path{1}).dim;

            X = sdpvar(dim, n, 'full');
            Constraints = [];
            for i = 1:n
                Xv = g.vertices(path{i});
                Constraints = [Constraints, Xv.constraints(X(:,i))]; %#ok<AGROW>
            end

            Objective = 0;
            for i = 1:n-1
                Objective = Objective + norm(X(:,i+1) - X(:,i), 2);
            end

            ops = sdpsettings('verbose', double(verbose));
            diagnostics = optimize(Constraints, Objective, ops);
            if diagnostics.problem ~= 0
                warning('ShortestPathProblem:restriction', ...
                    'Solver reported an issue solving the restriction: %s', ...
                    diagnostics.info);
            end

            positions = value(X);
            cost = value(Objective);
        end
    end

    methods (Access = private)
        function path = dfs_round(obj, g, current, y_val, visited, prefix)
            prefix{end+1} = current; %#ok<AGROW>
            if strcmp(current, obj.target)
                path = prefix;
                return
            end
            visited(current) = true; %#ok<NASGU>

            out_idx = g.outgoing_edges(current);
            weights = y_val(out_idx);
            [~, order] = sort(weights, 'descend');
            out_idx = out_idx(order);

            path = {};
            for k = 1:numel(out_idx)
                e = out_idx(k);
                head = g.edges{e, 2};
                if isKey(visited, head) && visited(head)
                    continue
                end
                candidate = obj.dfs_round(g, head, y_val, visited, prefix);
                if ~isempty(candidate)
                    path = candidate;
                    return
                end
            end
        end
    end
end

%Auxiliary Functions
function s = sum_or_zero(y, idx)
    if isempty(idx)
        s = 0;
    else
        s = sum(y(idx));
    end
end

function s = sum_cols_or_zero(P, idx, dim)
    if isempty(idx)
        s = zeros(dim, 1);
    else
        s = sum(P(:, idx), 2);
    end
end
