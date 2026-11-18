function [los, info] = line_of_sight(p, q, obstacles, opts)
%LINE_OF_SIGHT  Visibility test between two points among N-D obstacles.
%
%   [los, info] = line_of_sight(p, q, obstacles)
%   [los, info] = line_of_sight(p, q, obstacles, opts)
%
%   los = true  <=> the closed segment [p,q] does NOT cross the (eroded)
%                   interior of any obstacle.
%
% INPUTS
%   p, q       n-vector points (row or column).
%   obstacles  Either
%                - numeric m-by-n vertex matrix            : ONE convex obstacle
%                  (the convex hull of the rows is used),
%                - cell array whose elements are m_k-by-n vertex matrices
%                  or nested cell arrays of them,
%                  OR
%                - model struct array returned by prepare_los_model
%                  (precomputed H-representations; use for many queries).
%              A NON-CONVEX obstacle must be supplied as a cell array of
%              convex components (its union = the obstacle). A bare vertex
%              set of a non-convex body carries no connectivity, so its
%              interior is undefined; convex decomposition is a mandatory
%              pre-step. Segment-vs-union is exact: blocked iff blocked by
%              any component.
%   opts       struct, all fields optional:
%     .method      'auto' (default) | 'closed' | 'lp' | 'both'
%                  'closed' : H-representation (convhulln) + Cyrus-Beck
%                             clipping. Closed-form, O(#facets) per query.
%                  'lp'     : V-representation, one linprog feasibility
%                             problem per component. No facet enumeration.
%                  'both'   : runs both, warns on disagreement, returns
%                             the closed-form decision.
%                  'auto'   : 'closed' if n <= opts.maxHullDim, else 'lp'.
%     .maxHullDim  default 5. Facet count of a hull grows combinatorially
%                  with n; above this LP is the better method.
%     .shrink      default 1-1e-7. Each component is scaled about its
%                  centroid by this factor before testing, so grazing
%                  contact (touching a vertex/edge/face) does not block.
%                  Set to 1 to treat obstacles as closed sets
%                  (touching blocks).
%     .bboxPad     default 1e-9. Padding of the axis-aligned bounding-box
%                  prefilter (slab test).
%
% OUTPUT
%   los   logical.
%   info  struct: method, blockingObstacle (top-level index), 
%         blockingComponent (flattened index), tInterval ([t_lo t_hi] of
%         the segment inside the blocking component, closed-form only),
%         nTested (components surviving the bbox prefilter),
%         disagreement (method 'both'), degenerateFallback (a component
%         with empty interior was tested by LP).
%
% FORMULATION
%   Segment x(t) = p + t d, d = q - p, t in [0,1]. Convex component
%   P = {x : A x <= b}, outward unit normals A.
%   Closed form (Cyrus-Beck): each facet gives a half-line in t,
%       a_j'(p + t d) <= b_j   =>   t <= (b_j - a_j'p)/(a_j'd)  if a_j'd > 0
%                                   t >= (b_j - a_j'p)/(a_j'd)  if a_j'd < 0
%   t_lo = max(0, lower bounds), t_hi = min(1, upper bounds); the segment
%   crosses int(P) iff t_lo < t_hi. Parallel facets (a_j'd ~ 0) reject
%   when p is not strictly inside that facet.
%   LP: exists (w, t) :  V' w - t d = p,  sum(w) = 1,  w >= 0,  0 <= t <= 1
%   i.e. some point of the segment is a convex combination of the vertices.
%
% NOTES
%   - Endpoint inside an obstacle => blocked (t = 0 or 1 lies in P).
%   - Component with empty interior (rank(V - mean) < n, e.g. a flat
%     polygon in 3-D): no H-representation exists; tested by LP with the
%     closed set, so transversal crossings block.
%   - Convex-hull facets via Qhull (convhulln): requires n >= 2 for the
%     hull path; n = 1 is handled as an interval.

    if nargin < 4 || isempty(opts), opts = struct(); end
    opts = set_defaults(opts);

    p = p(:);
    q = q(:);
    n = numel(p);
    if numel(q) ~= n
        error('line_of_sight:dimMismatch', 'p and q must have the same dimension.');
    end
    if ~all(isfinite([p; q]))
        error('line_of_sight:nonFinite', 'p and q must be finite.');
    end

    if isstruct(obstacles)
        model = obstacles;
    else
        model = prepare_los_model(obstacles, n);
    end

    info = struct('method', '', 'blockingObstacle', [], 'blockingComponent', [], ...
                  'tInterval', [], 'nTested', 0, 'disagreement', false, ...
                  'degenerateFallback', false);

    method = opts.method;
    if strcmp(method, 'auto')
        if n <= opts.maxHullDim, method = 'closed'; else, method = 'lp'; end
    end
    info.method = method;

    los = true;
    if isempty(model), return; end
    if model(1).dim ~= n
        error('line_of_sight:dimMismatch', ...
              'Obstacle dimension (%d) differs from point dimension (%d).', model(1).dim, n);
    end

    d = q - p;
    for k = 1:numel(model)
        m = model(k);
        if slab_reject(p, d, m.lo, m.hi, opts.bboxPad)
            continue;
        end
        info.nTested = info.nTested + 1;
        tint = [];
        fb = false;
        switch method
            case 'closed'
                [hit, tint, fb] = test_closed(p, d, m, opts);
            case 'lp'
                hit = test_lp(p, d, m, opts);
            case 'both'
                [hc, tint, fb] = test_closed(p, d, m, opts);
                hl = test_lp(p, d, m, opts);
                if hc ~= hl && ~fb
                    warning('line_of_sight:methodDisagreement', ...
                        'Closed-form and LP disagree on component %d (closed=%d, lp=%d).', k, hc, hl);
                    info.disagreement = true;
                end
                hit = hc;
            otherwise
                error('line_of_sight:badMethod', 'Unknown method "%s".', opts.method);
        end
        if fb, info.degenerateFallback = true; end
        if hit
            los = false;
            info.blockingObstacle = m.owner;
            info.blockingComponent = k;
            info.tInterval = tint;
            return;
        end
    end
end

% ======================================================================
function [hit, tint, fallback] = test_closed(p, d, m, opts)
%   Cyrus-Beck clipping of [p, p+d] against the eroded polytope.
    tint = [];
    fallback = false;
    if ~m.fullDim
        hit = test_lp(p, d, m, opts);
        fallback = true;
        return;
    end
    s = opts.shrink;
    Ac = m.A * m.c;
    beff = Ac + s * (m.b - Ac);          % scaling about centroid
    num = beff - m.A * p;
    den = m.A * d;
    par = abs(den) <= 1e-12 * max(norm(d), realmin);

    if any(par & num <= 0)
        hit = false;
        return;
    end
    tlo = 0; thi = 1;
    pos = ~par & den > 0;
    neg = ~par & den < 0;
    if any(pos), thi = min(thi, min(num(pos) ./ den(pos))); end
    if any(neg), tlo = max(tlo, max(num(neg) ./ den(neg))); end
    hit = tlo < thi;
    if hit, tint = [tlo, thi]; end
end

% ======================================================================
function hit = test_lp(p, d, m, opts)
%   Feasibility of  V' w - t d = p, sum(w) = 1, w >= 0, 0 <= t <= 1.
    persistent lpOpt
    if isempty(lpOpt)
        lpOpt = optimoptions('linprog', 'Display', 'off');
    end
    V = m.V;
    Vs = m.c' + opts.shrink * (V - m.c');
    mv = size(V, 1);
    f = zeros(mv + 1, 1);
    Aeq = [Vs', -d; ones(1, mv), 0];
    beq = [p; 1];
    lb = [zeros(mv, 1); 0];
    ub = [inf(mv, 1); 1];
    [~, ~, flag] = linprog(f, [], [], Aeq, beq, lb, ub, lpOpt);
    hit = (flag == 1);
end

% ======================================================================
function rej = slab_reject(p, d, lo, hi, pad)
%   Segment vs axis-aligned box (slab method). true => segment misses box.
    pd = pad * max(1, max(hi - lo));
    lo = lo - pd; hi = hi + pd;
    t0 = 0; t1 = 1;
    rej = false;
    for i = 1:numel(p)
        if abs(d(i)) < 1e-15
            if p(i) < lo(i) || p(i) > hi(i)
                rej = true; return;
            end
        else
            ta = (lo(i) - p(i)) / d(i);
            tb = (hi(i) - p(i)) / d(i);
            if ta > tb, tmp = ta; ta = tb; tb = tmp; end
            t0 = max(t0, ta);
            t1 = min(t1, tb);
            if t0 > t1
                rej = true; return;
            end
        end
    end
end

% ======================================================================
function opts = set_defaults(opts)
    def = struct('method', 'auto', 'maxHullDim', 5, 'shrink', 1 - 1e-7, 'bboxPad', 1e-9);
    fn = fieldnames(def);
    for i = 1:numel(fn)
        if ~isfield(opts, fn{i}) || isempty(opts.(fn{i}))
            opts.(fn{i}) = def.(fn{i});
        end
    end
    if ~(opts.shrink > 0 && opts.shrink <= 1)
        error('line_of_sight:badShrink', 'opts.shrink must be in (0,1].');
    end
end
