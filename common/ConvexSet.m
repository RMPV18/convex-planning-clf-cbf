classdef ConvexSet
    % CONVEXSET  A polyhedral convex set {x : A*x <= b}.
    %
    % This plays the role of the Python repo's convex-set classes
    % (spp.convex_sets). Only polyhedra are implemented here (boxes,
    % polygons, points, halfspace intersections), which covers most
    % motion-planning-style examples and keeps every constraint an
    % ordinary linear/second-order-cone constraint that any conic
    % solver (Mosek, Gurobi, SCS, SDPT3, ECOS) can handle through
    % YALMIP.
    %
    % Interface used by ShortestPathProblem / GcsMotionPlanner (another
    % convex-set class, e.g. ellipsoids, only needs the same members):
    %   constraints(x)               -> x in X               (YALMIP)
    %   perspective_constraints(x,y) -> x/y in X, y >= 0      (YALMIP)
    %   contains(x_val, tol)         -> numeric membership test
    %   dim                          -> ambient dimension
    % Constructors: ConvexSet(A,b), ConvexSet.box(lb,ub),
    % ConvexSet.point(p), ConvexSet.polygon(V) (convex hull of rows of V).

    properties
        A double
        b double
        dim double
    end

    methods
        function obj = ConvexSet(A, b)
            % A: m x n, b: m x 1, set is {x in R^n : A*x <= b}
            if size(A,1) ~= numel(b)
                error('ConvexSet:dim', 'A and b row counts must match.');
            end
            obj.A = A;
            obj.b = b(:);
            obj.dim = size(A,2);
        end

        function cons = constraints(obj, x)
            % x: sdpvar (dim x 1). Returns YALMIP constraint x in X.
            cons = [obj.A * x <= obj.b];
        end

        function cons = perspective_constraints(obj, x, y)
            % Homogenization of X scaled by scalar y >= 0:
            %   x/y in X   <=>   A*x <= y*b,  y >= 0
            % This is exactly the perspective operator used in the
            % paper's MICP formulation (Section 4), and is what lets
            % the edge-selection binaries y_e "turn off" a set's
            % constraints when y_e = 0 (x forced to 0) while keeping
            % the constraint set convex for 0 <= y_e <= 1.
            cons = [obj.A * x <= obj.b * y, y >= 0];
        end

        function tf = contains(obj, x_val, tol)
            % Numeric membership check (for validation/plots), not YALMIP.
            if nargin < 3
                tol = 1e-6;
            end
            tf = all(obj.A * x_val(:) <= obj.b + tol);
        end
    end

    methods (Static)
        function cs = box(lb, ub)
            % Axis-aligned box [lb, ub] (column vectors of equal length).
            lb = lb(:); ub = ub(:);
            n = numel(lb);
            A = [eye(n); -eye(n)];
            b = [ub; -lb];
            cs = ConvexSet(A, b);
        end

        function cs = point(p)
            % Singleton {p}, encoded as a degenerate box.
            p = p(:);
            cs = ConvexSet.box(p, p);
        end

        function cs = polygon(V, tol)
            % Convex hull of the rows of V (m x n, any n) as {x : A x <= b}
            % with unit outward normals; coplanar facets merged within tol.
            if nargin < 2, tol = 1e-8; end
            n = size(V,2);  c = mean(V,1);
            T = convhulln(V, {'Qt'});                 % triangulated boundary facets
            A = zeros(size(T,1), n);  b = zeros(size(T,1),1);
            for f = 1:size(T,1)
                N = null([V(T(f,:),:), ones(n,1)]);   % [a; -off] with a'x = off on the facet
                N = N(:,1);  a = N(1:n);  off = -N(n+1);
                s = norm(a);  a = a/s;  off = off/s;
                if a'*c' > off, a = -a; off = -off; end   % orient outward
                A(f,:) = a';  b(f) = off;
            end
            [~, ia] = uniquetol([A b], max(tol,1e-8), 'ByRows', true);  % merge coplanar triangles
            A = A(ia,:);  b = b(ia);
            cs = ConvexSet(A, b);
        end
    end
end
