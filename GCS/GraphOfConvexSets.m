classdef GraphOfConvexSets < handle
    % GRAPHOFCONVEXSETS  Directed graph whose vertices carry ConvexSet
    % regions. Mirrors spp.graph.GraphOfConvexSets from the Python repo.

    properties
        vertices containers.Map      % name (char) -> ConvexSet
        vertex_names cell            % preserves insertion order
        edges cell = cell(0,2)       % Nx2 cell of {tail_name, head_name}
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
            idx = find(cellfun(@(t) strcmp(t, v), obj.edges(:,1)));
        end

        function idx = incoming_edges(obj, v)
            idx = find(cellfun(@(h) strcmp(h, v), obj.edges(:,2)));
        end

        function n = num_vertices(obj)
            n = numel(obj.vertex_names);
        end

        function n = num_edges(obj)
            n = size(obj.edges, 1);
        end
    end

    methods (Static)
        function g = add_mygraph(s,t, freecs, nb)
            g = GraphOfConvexSets();
            g.add_vertex('s', ConvexSet.point(s));
            g.add_vertex('t', ConvexSet.point(t));
            
            %Add Vertices as the intersection of the adjacent sets
            for i = 1:numel(freecs)
                for j = nb{i}
                    if i<j
                        g.add_vertex([num2str(i) '-' num2str(j)], ConvexSet([freecs{i}.A; freecs{j}.A],[freecs{i}.b; freecs{j}.b]));

                    end
                end
            end
            
            %Add Edges between the intersection sets and the other
            %intersection sets with a intersecting set in common 
            for i = 1:numel(freecs)
                for j = nb{i}
                    for k = nb{j}(find(nb{j} ~= i))
                        m = min(i,j);
                        n = max(i,j);
                        p = min(j,k);
                        q = max(j,k);
                        g.add_edge([num2str(m) '-' num2str(n)], [num2str(p) '-' num2str(q)]);
                    end
                end
            end
            
            % s & t edges
            for k = 1:numel(freecs)
                if freecs{k}.contains(s, 0)
                    for i = nb{k}, g.add_edge('s', [num2str(min(k,i)) '-' num2str(max(k,i))]); end
                end
                if freecs{k}.contains(t, 0)
                    for i = nb{k}, g.add_edge([num2str(min(k,i)) '-' num2str(max(k,i))], 't'); end
                end
            end
        end
    end
end






% function g = add_mygraph(s,t, freecs, nb)
%             g = GraphOfConvexSets();
%             g.add_vertex('s', ConvexSet.point(s));
%             g.add_vertex('t', ConvexSet.point(t));
% 
%             for i = 1:numel(freecs)
%                 if numel(nb{i}) > 3 
%                     g.add_vertex(num2str(i), freecs{i}); %Decreases the total edges
%                 end
%                 for j = nb{i}
%                     if i<j
%                         g.add_vertex([num2str(i) '-' num2str(j)], ConvexSet([freecs{i}.A; freecs{j}.A],[freecs{i}.b; freecs{j}.b]));
%                     end
%                 end
%             end
% 
%             for i = 1:numel(freecs)
%                 for j = nb{i}
%                     if any(strcmp(g.vertex_names, num2str(j))) %Edges Reduction
%                         m = min(i,j);
%                         n = max(i,j);
%                         g.add_edge([num2str(m) '-' num2str(n)], num2str(j));
%                         for k = nb{j}(find(nb{j} ~= i)) 
%                             p = min(j,k);
%                             q = max(j,k);
%                             g.add_edge(num2str(j),[num2str(p) '-' num2str(q)]);
%                         end
%                     else
%                         for k = nb{j}(find(nb{j} ~= i))
%                             m = min(i,j);
%                             n = max(i,j);
%                             p = min(j,k);
%                             q = max(j,k);
%                             g.add_edge([num2str(m) '-' num2str(n)], [num2str(p) '-' num2str(q)]);
%                         end
%                     end
%                 end
%             end
% 
%             % s & t edges
%             for k = 1:numel(freecs)
%                 if freecs{k}.contains(s, 0)
%                     if any(strcmp(g.vertex_names, num2str(k)))
%                         g.add_edge('s', num2str(k));
%                     else
%                         for i = nb{k}, g.add_edge('s', [num2str(min(k,i)) '-' num2str(max(k,i))]); end
%                     end
%                 end
%                 if freecs{k}.contains(t, 0)
%                     if any(strcmp(g.vertex_names, num2str(k)))
%                         g.add_edge(num2str(k),'t');
%                     else
%                         for i = nb{k}, g.add_edge([num2str(min(k,i)) '-' num2str(max(k,i))], 't'); end
%                     end
%                 end
%             end
%         end
% 
