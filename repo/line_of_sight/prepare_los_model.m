function model = prepare_los_model(obstacles, n)
%PREPARE_LOS_MODEL  Flatten obstacles into convex components and build
%   H-representations once (A x <= b, outward unit normals).
%   model(k): V (hull vertices), A, b, c (centroid), lo/hi (bbox), owner
%   (top-level obstacle index), fullDim, dim.

    if isnumeric(obstacles)
        obstacles = {obstacles};
    end
    if ~iscell(obstacles)
        error('prepare_los_model:badInput', 'obstacles must be numeric or cell.');
    end
    comps = {};
    owner = [];
    for k = 1:numel(obstacles)
        [cc, ~] = flatten(obstacles{k});
        comps = [comps, cc];                       %#ok<AGROW>
        owner = [owner, k * ones(1, numel(cc))];   %#ok<AGROW>
    end

    model = struct('V', {}, 'A', {}, 'b', {}, 'c', {}, 'lo', {}, 'hi', {}, ...
                   'owner', {}, 'fullDim', {}, 'dim', {});
    for i = 1:numel(comps)
        V = comps{i};
        if ~isnumeric(V) || ~ismatrix(V) || isempty(V) || ~all(isfinite(V(:)))
            error('prepare_los_model:badVertices', 'Component %d is not a finite numeric matrix.', i);
        end
        if nargin >= 2 && ~isempty(n) && size(V, 2) ~= n
            error('prepare_los_model:dimMismatch', ...
                  'Component %d has %d columns, expected %d.', i, size(V, 2), n);
        end
        V = unique(V, 'rows');
        dim = size(V, 2);
        m = struct('V', V, 'A', [], 'b', [], 'c', mean(V, 1)', ...
                   'lo', min(V, [], 1)', 'hi', max(V, [], 1)', ...
                   'owner', owner(i), 'fullDim', false, 'dim', dim);

        if dim == 1
            lo = m.lo; hi = m.hi;
            m.A = [1; -1];
            m.b = [hi; -lo];
            m.c = (lo + hi) / 2;
            m.V = [lo; hi];
            m.fullDim = hi > lo;
        else
            sv = svd(V - mean(V, 1), 0);
            tol = max(size(V)) * eps(max([sv(:); 1]));
            full = size(V, 1) >= dim + 1 && sum(sv > tol) == dim;
            if full
                try
                    K = convhulln(V);
                catch
                    K = [];
                    full = false;
                end
            end
            if full
                idx = unique(K(:));
                Vh = V(idx, :);
                c = mean(Vh, 1)';
                [A, b] = facets_from_hull(V, K, c);
                if isempty(A)
                    full = false;
                else
                    m.V = Vh; m.c = c; m.A = A; m.b = b;
                end
            end
            m.fullDim = full;
        end
        model(end + 1) = m; %#ok<AGROW>
    end
end

% ======================================================================
function [comps, depth] = flatten(x)
    if isnumeric(x)
        comps = {x}; depth = 0;
    elseif iscell(x)
        comps = {}; depth = 1;
        for i = 1:numel(x)
            [cc, dd] = flatten(x{i});
            comps = [comps, cc];   %#ok<AGROW>
            depth = max(depth, dd + 1);
        end
    else
        error('prepare_los_model:badInput', 'Unsupported obstacle element of class %s.', class(x));
    end
end

% ======================================================================
function [A, b] = facets_from_hull(V, K, c)
%   Hyperplane through each hull simplex, oriented outward, deduplicated
%   (coplanar facets are triangulated by Qhull).
    nf = size(K, 1);
    n = size(V, 2);
    A = zeros(nf, n);
    b = zeros(nf, 1);
    ok = false(nf, 1);
    for f = 1:nf
        F = V(K(f, :), :);
        M = F(2:end, :) - F(1, :);
        a = null(M);
        if size(a, 2) ~= 1, continue; end
        a = a(:)';
        bf = a * F(1, :)';
        if a * c > bf, a = -a; bf = -bf; end   % centroid must satisfy a'c <= b
        A(f, :) = a;
        b(f) = bf;
        ok(f) = true;
    end
    A = A(ok, :); b = b(ok);
    if isempty(A), return; end
    [~, ia] = unique(round([A, b] * 1e8), 'rows', 'stable');
    A = A(ia, :); b = b(ia);
end

