classdef GraphOfConvexSets < handle
    % GRAPHOFCONVEXSETS  Directed graph whose vertices carry ConvexSet
    % objects. Mirrors spp.graph.GraphOfConvexSets of the Python repo of
    % Marcucci et al. (arXiv:2101.11565).
    %
    %   g = GraphOfConvexSets();
    %   g.add_vertex(name, convex_set);  g.add_edge(tail, head);
    %
    % Vertices are addressed by name (char); edges are stored in insertion
    % order as an E x 2 cell {tail, head}, so edge index e is stable.

    properties
        vertices containers.Map      % name (char) -> ConvexSet
        vertex_names cell            % names in insertion order
        edges cell = cell(0,2)       % E x 2 cell of {tail_name, head_name}
    end

    methods
        function obj = GraphOfConvexSets()
            obj.vertices = containers.Map('KeyType', 'char', 'ValueType', 'any');
            obj.vertex_names = {};
        end

        function add_vertex(obj, name, convex_set)
            if isKey(obj.vertices, name)
                error('GraphOfConvexSets:dup', 'Vertex "%s" already exists.', name);
            end
            obj.vertices(name) = convex_set;
            obj.vertex_names{end+1} = name;
        end

        function add_edge(obj, tail, head)
            if ~isKey(obj.vertices, tail) || ~isKey(obj.vertices, head)
                error('GraphOfConvexSets:novertex', ...
                    'Both endpoints must be added as vertices first.');
            end
            obj.edges(end+1, :) = {tail, head};
        end

        function idx = outgoing_edges(obj, v)
            idx = find(strcmp(obj.edges(:,1), v));
        end

        function idx = incoming_edges(obj, v)
            idx = find(strcmp(obj.edges(:,2), v));
        end

        function n = num_vertices(obj)
            n = numel(obj.vertex_names);
        end

        function n = num_edges(obj)
            n = size(obj.edges, 1);
        end
    end
end
