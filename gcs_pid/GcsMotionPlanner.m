classdef GcsMotionPlanner < ShortestPathProblem
    % GCSMOTIONPLANNER  Collision-free trajectory design with Graphs of
    % Convex Sets. Extends ShortestPathProblem (YALMIP relaxation ->
    % rounding -> restriction) to the motion-planning SPP of
    %
    %   T. Marcucci, M. Petersen, D. von Wrangel, R. Tedrake, "Motion
    %   Planning around Obstacles with Convex Optimization,"
    %   arXiv:2205.04422v1 / Sci. Robot. 8, eadf7843 (2023).
    %   (Section / equation numbers below refer to the arXiv text.)
    %
    % REUSED: ConvexSet, GraphOfConvexSets, ShortestPathProblem (common/).
    %
    % MAPPING TO THE BASE CLASSES
    %   - Vertex i (name num2str(i)) holds a ConvexSet X_i in R^m,
    %     m=(n+1)(d+1): the lifted polytope of Bezier control points
    %     x_i=[r_0;...;r_d;h_0;...;h_d], eq. (7a)-(7d). Its homogenization
    %     ConvexSet.perspective_constraints is the perspective used by the
    %     relaxation.
    %   - 's' and 't' are dim-0 ConvexSets (X_sigma = X_tau = empty, Sec. 5.2).
    %   - Edge sets X_e (boundary conditions, continuity (8)) and edge costs
    %     (9)-(10) are not representable by the base classes; they are added
    %     here in edge_constraints / edge_cost, shared by the relaxation
    %     (phi = y_e, perspective variables pu/pv) and the restriction
    %     (phi = 1, actual control points).
    %   - round_solution is replaced by the randomized DFS of Sec. 4.2
    %     (greedy rounding is explicitly discouraged there).
    %
    % ADDED vs. base ShortestPathProblem: Bezier trajectory variables,
    % duration/length/energy costs, acceleration regularizer (11), two-cycle
    % constraints (12), multiflow edge pruning (App. A.2), N-path rounding
    % with early global-optimality stop, trajectory sampling, verification.
    %
    % PARAMS (struct, all fields optional; defaults in brackets)
    %   a,b,c [0,1,0]   weights: duration, length, energy in (1a)
    %   eps [0]         acceleration regularizer weight (11), needs d>=2
    %   eta [0], d [1]  smoothness order and Bezier degree, d >= eta+1
    %   vmax [Inf]      box velocity limit D=[-vmax,vmax]^n
    %   hdotMin [1e-6]  lower bound on time-scaling derivative (7b)
    %   Tmin [0], Tmax [1e3]
    %   v0, vT []       boundary velocities (empty = free)
    %   nPaths [10], maxTrials [100]   rounding parameters N, M
    %   phiTol [1e-5], gapTol [1e-6], seed [0]
    %   preprocess [true]  multiflow LP edge pruning (App. A.2)
    %   minOverlap [0]  minimum Chebyshev overlap radius for an edge
    %   nSample [50], verifyTol [1e-6], solver ['' = YALMIP default]
    %
    % LIMITS: only the linear two-cycle constraints (12) are used (conic
    % tightening of Lemma 1(b), arXiv:2101.11565, is not in the project);
    % regions must be bounded polytopes; boundary conditions on derivatives
    % of order >= 2 are not imposed.
    %
    % Requires YALMIP + conic solver (Mosek recommended) and the
    % Optimization Toolbox (linprog, for overlap tests and edge pruning).

    properties
        regions cell      % configuration-space polytopes Q_i (ConvexSet)
        q0 double
        qT double
        params struct
        n double
        d double
        m double
    end

    methods
        function obj = GcsMotionPlanner(regions, q0, qT, params)
            if nargin < 4, params = struct(); end
            q0 = q0(:);  qT = qT(:);  n = numel(q0);
            params = GcsMotionPlanner.defaults(params, n);
            GcsMotionPlanner.validate(regions, q0, qT, params, n);
            g = GcsMotionPlanner.build_graph(regions, q0, qT, params, n);
            obj@ShortestPathProblem(g, 's', 't');
            obj.regions = regions;  obj.q0 = q0;  obj.qT = qT;  obj.params = params;
            obj.n = n;  obj.d = params.d;  obj.m = (n + 1) * (params.d + 1);
        end

        function [path, positions, cost, info] = solve(obj, verbose)
            % relax -> randomized rounding (N paths) -> restriction per path
            if nargin < 2, verbose = false; end
            o = obj.params;
            if o.preprocess, obj.prune_edges(); end
            [y_val, relax_cost] = obj.solve_relaxation(verbose);
            if any(isnan(y_val))
                error('GcsMotionPlanner:relaxation', 'Relaxation returned no solution.');
            end
            rs = RandStream('mt19937ar', 'Seed', o.seed);
            keys = {};  costs = [];  best = Inf;  positions = [];  path = {};
            trial = 0;
            while numel(keys) < o.nPaths && trial < o.maxTrials
                trial = trial + 1;
                p = obj.round_solution(y_val, rs);
                if isempty(p), continue; end
                key = strjoin(p, '>');
                if any(strcmp(keys, key)), continue; end
                keys{end+1} = key; %#ok<AGROW>
                [X, c] = obj.solve_restriction(p, verbose);
                costs(end+1) = c; %#ok<AGROW>
                if c < best, best = c;  positions = X;  path = p; end
                if best <= relax_cost + o.gapTol * max(1, abs(relax_cost)), break; end
            end
            if ~isfinite(best)
                error('GcsMotionPlanner:rounding', ...
                    'No rounded path with finite cost found in %d trials.', trial);
            end
            cost = best;
            info = struct();
            info.relaxation_cost = relax_cost;                    % lower bound
            info.restriction_cost = cost;                         % feasible upper bound
            info.optimality_gap = cost - relax_cost;
            info.relative_gap = (cost - relax_cost) / relax_cost; % delta_relax (Sec. 4.2)
            info.edge_weights = y_val;
            info.rounding_costs = costs;
            info.trajectory = obj.extract_trajectory(path, positions);
        end

        function [y_val, cost] = solve_relaxation(obj, verbose)
            % Perspective relaxation of the MICP; costs/constraints of Sec. 5.
            if nargin < 2, verbose = false; end
            g = obj.graph;  nE = g.num_edges();  m = obj.m;
            y  = sdpvar(nE, 1, 'full');
            pu = sdpvar(m, nE, 'full');   % perspective copy of tail control points
            pv = sdpvar(m, nE, 'full');   % perspective copy of head control points
            Constraints = [0 <= y <= 1];
            Objective = 0;
            for e = 1:nE
                tail = g.edges{e,1};  head = g.edges{e,2};
                if ~strcmp(tail, 's')
                    Constraints = [Constraints, g.vertices(tail).perspective_constraints(pu(:,e), y(e))]; %#ok<AGROW>
                    [oc, cc] = obj.edge_cost(pu(:,e), y(e));
                    Objective = Objective + oc;
                    Constraints = [Constraints, cc]; %#ok<AGROW>
                end
                if ~strcmp(head, 't')
                    Constraints = [Constraints, g.vertices(head).perspective_constraints(pv(:,e), y(e))]; %#ok<AGROW>
                end
                Constraints = [Constraints, obj.edge_constraints(tail, head, pu(:,e), pv(:,e), y(e))]; %#ok<AGROW>
            end
            % probability flow and spatial flow conservation
            for i = 1:g.num_vertices()
                v = g.vertex_names{i};
                out_idx = g.outgoing_edges(v);  in_idx = g.incoming_edges(v);
                out_sum = obj.sum_or_zero(y, out_idx);  in_sum = obj.sum_or_zero(y, in_idx);
                if strcmp(v, obj.source)
                    Constraints = [Constraints, out_sum - in_sum == 1]; %#ok<AGROW>
                elseif strcmp(v, obj.target)
                    Constraints = [Constraints, out_sum - in_sum == -1]; %#ok<AGROW>
                else
                    Constraints = [Constraints, out_sum - in_sum == 0, in_sum <= 1]; %#ok<AGROW>
                    Constraints = [Constraints, ...
                        obj.sum_cols_or_zero(pu, out_idx, m) == obj.sum_cols_or_zero(pv, in_idx, m)]; %#ok<AGROW>
                end
            end
            % two-cycle elimination (12): y_e + y_f <= flow through i and j
            for e = 1:nE
                tail = g.edges{e,1};  head = g.edges{e,2};
                if any(strcmp({tail, head}, 's')) || any(strcmp({tail, head}, 't')), continue; end
                f = find(strcmp(g.edges(:,1), head) & strcmp(g.edges(:,2), tail) & (1:nE)' > e, 1);
                if isempty(f), continue; end
                Constraints = [Constraints, ...
                    y(e) + y(f) <= obj.sum_or_zero(y, g.incoming_edges(tail)), ...
                    y(e) + y(f) <= obj.sum_or_zero(y, g.incoming_edges(head))]; %#ok<AGROW>
            end
            diagnostics = optimize(Constraints, Objective, obj.yalmip_ops(verbose));
            if diagnostics.problem ~= 0
                warning('GcsMotionPlanner:relaxation', ...
                    'Solver reported an issue solving the relaxation: %s', diagnostics.info);
            end
            y_val = value(y);
            cost = value(Objective);
        end

        function path = round_solution(obj, y_val, rs)
            % Randomized DFS with backtracking (Sec. 4.2): from each vertex
            % traverse an unvisited out-edge with probability proportional
            % to y_e. Returns {} if the search exhausts the graph.
            if nargin < 3, rs = RandStream('mt19937ar', 'Seed', obj.params.seed); end
            g = obj.graph;  tol = obj.params.phiTol;
            visited = containers.Map('KeyType', 'char', 'ValueType', 'logical');
            dead = containers.Map('KeyType', 'char', 'ValueType', 'logical');
            path = {obj.source};  visited(obj.source) = true;
            while ~isempty(path)
                u = path{end};
                if strcmp(u, obj.target), return; end
                cand = [];
                for e = g.outgoing_edges(u)'
                    h = g.edges{e,2};
                    if y_val(e) > tol && ~isKey(visited, h) && ~isKey(dead, h)
                        cand(end+1) = e; %#ok<AGROW>
                    end
                end
                if isempty(cand)
                    dead(u) = true;  path(end) = [];  % backtrack
                    continue;
                end
                w = y_val(cand);
                k = find(rand(rs) * sum(w) <= cumsum(w), 1);
                h = g.edges{cand(k), 2};
                path{end+1} = h;  visited(h) = true; %#ok<AGROW>
            end
        end

        function [positions, cost] = solve_restriction(obj, path, verbose)
            % Program (6): path fixed, optimize control points only.
            % positions: m x (numel(path)-2), column k = x_i of k-th region.
            if nargin < 3, verbose = false; end
            g = obj.graph;  nR = numel(path) - 2;
            X = sdpvar(obj.m, nR, 'full');
            Constraints = [];  Objective = 0;
            for k = 2:numel(path)-1
                Constraints = [Constraints, g.vertices(path{k}).constraints(X(:,k-1))]; %#ok<AGROW>
                [oc, cc] = obj.edge_cost(X(:,k-1), 1);
                Objective = Objective + oc;
                Constraints = [Constraints, cc]; %#ok<AGROW>
            end
            for k = 1:numel(path)-1
                xu = [];  xv = [];
                if k >= 2, xu = X(:,k-1); end
                if k+1 <= numel(path)-1, xv = X(:,k); end
                Constraints = [Constraints, ...
                    obj.edge_constraints(path{k}, path{k+1}, xu, xv, 1)]; %#ok<AGROW>
            end
            diagnostics = optimize(Constraints, Objective, obj.yalmip_ops(verbose));
            if diagnostics.problem ~= 0
                warning('GcsMotionPlanner:restriction', ...
                    'Restriction infeasible or solver issue (path %s): %s', ...
                    strjoin(path, '>'), diagnostics.info);
                positions = [];  cost = Inf;
                return
            end
            positions = value(X);
            cost = value(Objective);
        end

        function prune_edges(obj)
            % App. A.2: drop edge (u,v) when the LP (fractional multiflow)
            % relaxation of "vertex-disjoint paths s->u and v->t" is
            % infeasible (sufficient condition for redundancy).
            g = obj.graph;  names = g.vertex_names;  V = numel(names);  E = g.num_edges();
            eu = zeros(E,1);  ev = zeros(E,1);
            for e = 1:E
                eu(e) = find(strcmp(names, g.edges{e,1}));
                ev(e) = find(strcmp(names, g.edges{e,2}));
            end
            sig = find(strcmp(names, obj.source));  tau = find(strcmp(names, obj.target));
            Inc = sparse([eu; ev], [1:E, 1:E]', [ones(E,1); -ones(E,1)], V, E);
            In = sparse(ev, (1:E)', 1, V, E);
            Aeq = blkdiag(Inc, Inc);  Ain = [In In];
            keep = true(E, 1);
            lpo = optimoptions('linprog', 'Display', 'off');
            for e = 1:E
                b1 = zeros(V,1);  b1(sig) = b1(sig) + 1;  b1(eu(e)) = b1(eu(e)) - 1;
                b2 = zeros(V,1);  b2(ev(e)) = b2(ev(e)) + 1;  b2(tau) = b2(tau) - 1;
                cap = ones(V,1);  cap(sig) = cap(sig) - 1;  cap(ev(e)) = cap(ev(e)) - 1;
                ub = inf(2*E, 1);  ub([e, E+e]) = 0;
                [~, ~, fl] = linprog(zeros(2*E,1), Ain, cap, Aeq, [b1; b2], zeros(2*E,1), ub, lpo);
                if fl == -2, keep(e) = false; end
            end
            g2 = GraphOfConvexSets();
            for i = 1:V, g2.add_vertex(names{i}, g.vertices(names{i})); end
            for e = find(keep)', g2.add_edge(g.edges{e,1}, g.edges{e,2}); end
            obj.graph = g2;
        end

        function T = extract_trajectory(obj, path, positions)
            % Sec. 5.5: sequence the Bezier segments of the selected regions.
            n = obj.n;  d = obj.d;  nR = numel(path) - 2;
            segs = repmat(struct('r', [], 'h', [], 'region', 0), nR, 1);
            for k = 1:nR
                x = positions(:, k);
                segs(k).r = reshape(x(1:n*(d+1)), n, d+1)';
                segs(k).h = x(n*(d+1)+1:end);
                segs(k).region = str2double(path{k+1});
            end
            T.segments = segs;
            T.sample = GcsMotionPlanner.sample(segs, obj.params.nSample);
            T.verification = obj.verify(segs);
        end

        function V = verify(obj, segs)
            % Containment, continuity, velocity, boundary and time-scaling checks.
            o = obj.params;  s = linspace(0, 1, 401)';  tol = o.verifyTol;
            maxViol = -Inf;  minHdot = Inf;  maxVel = -Inf;  maxGap = 0;
            for k = 1:numel(segs)
                cs = obj.regions{segs(k).region};
                r = GcsMotionPlanner.bez(segs(k).r, s);
                maxViol = max(maxViol, max(max(r * cs.A' - cs.b', [], 2)));
                hd = GcsMotionPlanner.bez(GcsMotionPlanner.bez_diff(segs(k).h), s);
                minHdot = min(minHdot, min(hd));
                if ~isempty(o.D)
                    rd = GcsMotionPlanner.bez(GcsMotionPlanner.bez_diff(segs(k).r), s);
                    maxVel = max(maxVel, max(max(rd * o.D.A' - hd * o.D.b', [], 2)));
                end
            end
            for k = 1:numel(segs)-1
                Ra = segs(k).r;  Rb = segs(k+1).r;  Ha = segs(k).h;  Hb = segs(k+1).h;
                for l = 0:o.eta
                    maxGap = max([maxGap, norm(Ra(end,:) - Rb(1,:)), abs(Ha(end) - Hb(1))]);
                    if l < o.eta
                        Ra = GcsMotionPlanner.bez_diff(Ra);  Rb = GcsMotionPlanner.bez_diff(Rb);
                        Ha = GcsMotionPlanner.bez_diff(Ha);  Hb = GcsMotionPlanner.bez_diff(Hb);
                    end
                end
            end
            V.maxRegionViolation = maxViol;       % <= 0: inside union of Q_i
            V.minHdot = minHdot;                  % > 0: time scaling increasing
            V.maxVelocityViolation = maxVel;      % <= 0 (or -Inf if no limit)
            V.maxContinuityGap = maxGap;
            V.startResidual = norm(segs(1).r(1,:)' - obj.q0);
            V.endResidual = norm(segs(end).r(end,:)' - obj.qT);
            V.duration = segs(end).h(end) - segs(1).h(1);
            V.pass = maxViol <= tol && minHdot > 0 && maxGap <= tol && ...
                V.startResidual <= tol && V.endResidual <= tol && ...
                (isempty(o.D) || maxVel <= tol);
            if ~V.pass
                warning('GcsMotionPlanner:verification', ...
                    ['Verification beyond tol %.1e: containment %.2e, continuity %.2e, ' ...
                     'velocity %.2e, boundary %.2e/%.2e, min hdot %.2e.'], tol, maxViol, ...
                     maxGap, maxVel, V.startResidual, V.endResidual, minHdot);
            end
        end
    end

    methods (Access = private)
        function ops = yalmip_ops(obj, verbose)
            if isempty(obj.params.solver)
                ops = sdpsettings('verbose', double(verbose));
            else
                ops = sdpsettings('verbose', double(verbose), 'solver', obj.params.solver);
            end
        end

        function [oc, cc] = edge_cost(obj, x, phi)
            % Edge length (9) with bounds (10), perspective in phi, plus the
            % regularizer (11). x: local control-point vector (scaled by phi).
            n = obj.n;  d = obj.d;  o = obj.params;  oc = 0;  cc = [];
            if o.a > 0
                oc = oc + o.a * (GcsMotionPlanner.Hop(n,d,0,d) - GcsMotionPlanner.Hop(n,d,0,0)) * x;
            end
            for k = 0:d-1
                Dr = (GcsMotionPlanner.Rop(n,d,0,k+1) - GcsMotionPlanner.Rop(n,d,0,k)) * x;
                if o.b > 0
                    t = sdpvar(1);
                    cc = [cc, cone(Dr, t)]; %#ok<AGROW>
                    oc = oc + o.b * t;
                end
                if o.c > 0
                    Dh = (GcsMotionPlanner.Hop(n,d,0,k+1) - GcsMotionPlanner.Hop(n,d,0,k)) * x;
                    s = sdpvar(1);   % s >= ||Dr||^2/Dh  <=>  ||[2Dr; s-Dh]|| <= s+Dh
                    cc = [cc, cone([2*Dr; s - Dh], s + Dh)]; %#ok<AGROW>
                    oc = oc + o.c * s;
                end
            end
            if o.eps > 0
                for k = 0:d-2
                    wr = GcsMotionPlanner.Rop(n,d,2,k) * x;
                    wh = GcsMotionPlanner.Hop(n,d,2,k) * x;
                    gg = sdpvar(1);  % gg >= ||w||^2/phi
                    cc = [cc, cone([2*wr; 2*wh; gg - phi], gg + phi)]; %#ok<AGROW>
                    oc = oc + o.eps / (d - 1) * gg;
                end
            end
        end

        function cc = edge_constraints(obj, tail, head, xu, xv, phi)
            % Edge sets X_e (Sec. 5.3): boundary conditions on (s,i), (i,t);
            % continuity (8) on (i,j). Homogeneous in phi (perspective).
            n = obj.n;  d = obj.d;  o = obj.params;  cc = [];
            R = @(l,k) GcsMotionPlanner.Rop(n,d,l,k);
            H = @(l,k) GcsMotionPlanner.Hop(n,d,l,k);
            if strcmp(tail, 's')
                cc = [R(0,0)*xv == obj.q0*phi, H(0,0)*xv == 0];
                if ~isempty(o.v0), cc = [cc, R(1,0)*xv == o.v0*(H(1,0)*xv)]; end
            elseif strcmp(head, 't')
                cc = [R(0,d)*xu == obj.qT*phi, H(0,d)*xu >= o.Tmin*phi];
                if ~isempty(o.vT), cc = [cc, R(1,d-1)*xu == o.vT*(H(1,d-1)*xu)]; end
            else
                for l = 0:o.eta
                    cc = [cc, R(l,d-l)*xu == R(l,0)*xv, H(l,d-l)*xu == H(l,0)*xv]; %#ok<AGROW>
                end
            end
        end
    end

    methods (Static, Access = private)
        function o = defaults(o, n)
            def = struct('a',0,'b',1,'c',0,'eps',0,'eta',0,'d',1,'vmax',Inf, ...
                'hdotMin',1e-6,'Tmin',0,'Tmax',1e3,'v0',[],'vT',[],'nPaths',10, ...
                'maxTrials',100,'phiTol',1e-5,'gapTol',1e-6,'seed',0,'preprocess',true, ...
                'minOverlap',0,'nSample',50,'verifyTol',1e-6,'solver','');
            fn = fieldnames(def);
            for i = 1:numel(fn)
                if ~isfield(o, fn{i}), o.(fn{i}) = def.(fn{i}); end
            end
            o.D = [];
            vm = o.vmax(:);
            if any(isfinite(vm))
                if isscalar(vm), vm = repmat(vm, n, 1); end
                o.D = struct('A', [eye(n); -eye(n)], 'b', [vm; vm]);
            end
            if ~isempty(o.v0), o.v0 = o.v0(:); end
            if ~isempty(o.vT), o.vT = o.vT(:); end
        end

        function validate(regions, q0, qT, o, n)
            assert(iscell(regions) && ~isempty(regions), 'GcsMotionPlanner:input', ...
                'regions must be a nonempty cell array of ConvexSet.');
            for i = 1:numel(regions)
                assert(regions{i}.dim == n, 'GcsMotionPlanner:input', ...
                    'Region %d has dimension %d, expected %d.', i, regions{i}.dim, n);
            end
            assert(numel(qT) == n, 'GcsMotionPlanner:input', 'q0 and qT dimension mismatch.');
            assert(o.d >= o.eta + 1, 'GcsMotionPlanner:input', 'Need d >= eta+1.');
            assert(o.eps == 0 || o.d >= 2, 'GcsMotionPlanner:input', 'Regularizer needs d >= 2.');
            assert(all([o.a o.b o.c o.eps] >= 0) && any([o.a o.b o.c o.eps] > 0), ...
                'GcsMotionPlanner:input', 'Cost weights must be >= 0 and not all zero.');
            assert(o.hdotMin > 0 && o.Tmax > o.Tmin && o.Tmin >= 0, ...
                'GcsMotionPlanner:input', 'Invalid time bounds.');
            assert(isempty(o.v0) || numel(o.v0) == n, 'GcsMotionPlanner:input', 'v0 size.');
            assert(isempty(o.vT) || numel(o.vT) == n, 'GcsMotionPlanner:input', 'vT size.');
        end

        function g = build_graph(regions, q0, qT, o, n)
            % Sec. 5.1: vertex per region, edge (i,j) iff Q_i and Q_j overlap.
            N = numel(regions);  d = o.d;
            g = GraphOfConvexSets();
            none = ConvexSet(zeros(0,0), zeros(0,1));   % X_sigma = X_tau = empty
            g.add_vertex('s', none);
            g.add_vertex('t', none);
            for i = 1:N
                g.add_vertex(num2str(i), GcsMotionPlanner.lifted_set(regions{i}, n, d, o));
            end
            lpo = optimoptions('linprog', 'Display', 'off');
            for i = 1:N-1
                for j = i+1:N
                    A = [regions{i}.A; regions{j}.A];  b = [regions{i}.b; regions{j}.b];
                    nr = sqrt(sum(A.^2, 2));
                    [~, fv, fl] = linprog([zeros(n,1); -1], [A nr], b, [], [], [], [], lpo);
                    if fl == -3
                        error('GcsMotionPlanner:unbounded', ...
                            'Regions %d,%d: unbounded polytope; Q_i must be bounded.', i, j);
                    end
                    if fl ~= 1, continue; end
                    r = -fv;  tol = 1e-7;
                    % Touching regions (r ~ 0, e.g. from free_space_decomposition)
                    % are connected only through a shared facet: a vertex-only
                    % contact is a zero-clearance pinch between obstacle corners.
                    if r >= max(o.minOverlap, tol) || ...
                       (o.minOverlap <= tol && r >= -tol && ...
                        GcsMotionPlanner.facet_contact(A ./ nr, b ./ nr, n, tol))
                        g.add_edge(num2str(i), num2str(j));
                        g.add_edge(num2str(j), num2str(i));
                    end
                end
            end
            hasS = false;  hasT = false;
            for i = 1:N
                if regions{i}.contains(q0, 1e-9), g.add_edge('s', num2str(i)); hasS = true; end
                if regions{i}.contains(qT, 1e-9), g.add_edge(num2str(i), 't'); hasT = true; end
            end
            if ~hasS, error('GcsMotionPlanner:q0', 'q0 is not contained in any region.'); end
            if ~hasT, error('GcsMotionPlanner:qT', 'qT is not contained in any region.'); end
        end

        function tf = facet_contact(A, b, n, tol)
            % True iff the touching polytopes {A x <= b} (unit-norm rows,
            % stacked [Q_i; Q_j]) share an (n-1)-dimensional facet patch:
            % some opposing row pair spans a hyperplane in which the
            % intersection contains an (n-1)-ball of radius > tol.
            lpo = optimoptions('linprog', 'Display', 'off');
            [fi, fj] = find(abs(A * A' + 1) < 1e-6 & abs(b + b') < 1e-6 * max(1, norm(b, inf)));
            tf = false;
            for f = find(fi < fj)'
                nrm = A(fi(f), :)';
                pn = vecnorm(A * (eye(n) - nrm * nrm'), 2, 2);   % in-plane row norms
                keep = pn > 1e-9;
                [~, fv, fl] = linprog([zeros(n,1); -1], [A(keep,:), pn(keep)], b(keep), ...
                    [nrm; 0]', b(fi(f)), [], [], lpo);
                if fl == 1 && -fv > tol, tf = true; return; end
            end
        end

        function cs = lifted_set(region, n, d, o)
            % X_i of eq. (7) as {x : A x <= b}, x=[r_0..r_d; h_0..h_d].
            rows = {};  rhs = {};
            for k = 0:d
                rows{end+1} = region.A * GcsMotionPlanner.Rop(n,d,0,k); rhs{end+1} = region.b; %#ok<AGROW> (7a)
            end
            for k = 0:d-1
                rows{end+1} = -GcsMotionPlanner.Hop(n,d,1,k);  rhs{end+1} = -o.hdotMin; %#ok<AGROW> (7b)
                if ~isempty(o.D)                                                         % (7c)
                    rows{end+1} = o.D.A * GcsMotionPlanner.Rop(n,d,1,k) ...
                        - o.D.b * GcsMotionPlanner.Hop(n,d,1,k);                         %#ok<AGROW>
                    rhs{end+1} = zeros(size(o.D.A, 1), 1);                               %#ok<AGROW>
                end
            end
            rows{end+1} = -GcsMotionPlanner.Hop(n,d,0,0);  rhs{end+1} = 0;      %#ok<AGROW> (7d)
            rows{end+1} = GcsMotionPlanner.Hop(n,d,0,d);   rhs{end+1} = o.Tmax; %#ok<AGROW>
            cs = ConvexSet(vertcat(rows{:}), vertcat(rhs{:}));
        end

        function c = deriv_coef(d, l)
            % l-th derivative control point = c * [gamma_k ... gamma_{k+l}]
            j = 0:l;
            c = factorial(d) / factorial(d - l) * ((-1).^(l - j)) .* arrayfun(@(q) nchoosek(l, q), j);
        end

        function M = Rop(n, d, l, k)
            % n x m operator: r^{(l)}_k
            c = GcsMotionPlanner.deriv_coef(d, l);
            M = zeros(n, (n+1)*(d+1));
            for j = 0:l
                cols = (k + j) * n + (1:n);
                M(:, cols) = M(:, cols) + c(j+1) * eye(n);
            end
        end

        function M = Hop(n, d, l, k)
            % 1 x m operator: h^{(l)}_k
            c = GcsMotionPlanner.deriv_coef(d, l);
            M = zeros(1, (n+1)*(d+1));
            for j = 0:l
                M(1, n*(d+1) + k + j + 1) = c(j+1);
            end
        end
    end

    methods (Static)
        function out = bez(C, s)
            % Evaluate Bezier curve with (d+1) x n control points at s in [0,1].
            d = size(C, 1) - 1;  s = s(:);
            B = zeros(numel(s), d + 1);
            for k = 0:d
                B(:, k+1) = nchoosek(d, k) * s.^k .* (1 - s).^(d - k);
            end
            out = B * C;
        end

        function C = bez_diff(C)
            % Derivative control points: d*(gamma_{k+1}-gamma_k).
            C = (size(C, 1) - 1) * diff(C, 1, 1);
        end

        function S = sample(segs, nPer)
            % Sampled q(t), qdot(t) = rdot(s)/hdot(s) over all segments.
            t = [];  q = [];  qd = [];
            for k = 1:numel(segs)
                s = linspace(0, 1, nPer + 1)';
                if k < numel(segs), s(end) = []; end
                r  = GcsMotionPlanner.bez(segs(k).r, s);
                rd = GcsMotionPlanner.bez(GcsMotionPlanner.bez_diff(segs(k).r), s);
                h  = GcsMotionPlanner.bez(segs(k).h, s);
                hd = GcsMotionPlanner.bez(GcsMotionPlanner.bez_diff(segs(k).h), s);
                t = [t; h]; q = [q; r]; qd = [qd; rd ./ hd]; %#ok<AGROW>
            end
            S.t = t;  S.q = q;  S.qdot = qd;
        end
    end
end
