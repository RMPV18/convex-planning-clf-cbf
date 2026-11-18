function [neighbors, A_adj] = adjacent_regions(regions, mode, tol)
%ADJACENT_REGIONS  Neighbor lists of convex polytopes given by their vertices.
%
%   [neighbors, A_adj] = adjacent_regions(regions)
%   [neighbors, A_adj] = adjacent_regions(regions, mode, tol)
%
%   regions   1 x K cell, regions{k} is m_k x n, rows = vertices of a
%             full-dimensional convex polytope (vertices may be shared)
%   mode      'facet' (default): shared (n-1)-dimensional contact
%             'any'            : any contact (facet, edge or vertex)
%   tol       numerical tolerance (default 1e-7)
%
%   neighbors 1 x K cell, neighbors{k} = indices of regions adjacent to k
%   A_adj     K x K logical adjacency matrix
%
%   Method: regions i, j touch iff the signed Chebyshev radius of
%   Q_i ∩ Q_j is ~0 (LP). For 'facet', an opposing facet pair must also
%   carry an (n-1)-ball inside the common hyperplane (second LP).

if nargin < 2 || isempty(mode), mode = 'facet'; end
if nargin < 3, tol = 1e-7; end
K = numel(regions);  n = size(regions{1},2);
lpo = optimoptions('linprog','Display','none');

% H-representation of each region: A x <= b, unit-norm rows
A = cell(1,K); b = cell(1,K);
for k = 1:K, [A{k}, b{k}] = vert2hrep(regions{k}, tol); end

A_adj = false(K);
for i = 1:K-1
    for j = i+1:K
        % bounding-box rejection
        if any(max(regions{i},[],1) < min(regions{j},[],1) - tol) || ...
           any(max(regions{j},[],1) < min(regions{i},[],1) - tol), continue, end

        Ao = [A{i}; A{j}];  bo = [b{i}; b{j}];

        % Signed Chebyshev radius of the intersection:
        %   r > 0 overlap, r ~ 0 touching, r < 0 separated
        [~, fv] = linprog([zeros(n,1); -1], [Ao, ones(size(Ao,1),1)], bo, ...
                          [],[],[],[],lpo);
        r = -fv;
        if r > tol || r < -tol, continue, end      % overlap or separated

        if strcmp(mode,'any'), A_adj(i,j) = true; continue, end

        % Facet contact: opposing facet hyperplanes (n_i = -n_j, b_i = -b_j)
        [fi, ~] = find(abs(A{i}*A{j}' + 1) < 1e-6 & ...
                       abs(b{i} + b{j}')   < 1e-6*max(1,norm(bo,inf)));
        for f = 1:numel(fi)
            nrm = A{i}(fi(f),:)';  off = b{i}(fi(f));
            pn  = vecnorm(Ao*(eye(n) - nrm*nrm'), 2, 2);   % in-plane normal lengths
            keep = pn > 1e-9;
            % Chebyshev ball inside the common hyperplane: the contact has
            % (n-1)-dim extent iff its radius is > 0
            [~, rr] = linprog([zeros(n,1); -1], [Ao(keep,:), pn(keep)], bo(keep), ...
                              [nrm; 0]', off, [],[],lpo);
            if ~isempty(rr) && -rr > tol
                A_adj(i,j) = true;  break
            end
        end
    end
end
A_adj = A_adj | A_adj';
neighbors = arrayfun(@(k) find(A_adj(k,:)), 1:K, 'UniformOutput', false);
end

% ---------------------------------------------------------------------
function [A, b] = vert2hrep(V, tol)
% Convex hull of the rows of V as {x : A x <= b}, unit outward normals.
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
end
