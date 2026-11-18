function [adjFacet, adjAny, info] = tangentRegions(regions, tol)
% TANGENTREGIONS  Find convex regions that are tangent (touching, not overlapping).
%
%   regions : 1xK cell, regions{k} is an (m_k x n) matrix of vertices of a
%             full-dimensional convex polytope (rows = vertices).
%   tol     : numerical tolerance (default 1e-7)
%
%   adjFacet: KxK logical, true if regions i,j share a facet portion
%             (contact set is (n-1)-dimensional).
%   adjAny  : KxK logical, true if regions touch in any way (facet, edge,
%             vertex...) with disjoint interiors.
%   info    : KxK struct array with fields
%               .type  : 'disjoint' | 'overlap' | 'touch' | 'facet'
%               .gap   : signed Chebyshev radius of the intersection
%                        (<0 separated, ~0 touching, >0 overlapping)
%               .n, .b : contact hyperplane n'x = b (only for 'facet')

if nargin < 2, tol = 1e-7; end
K = numel(regions);
n = size(regions{1}, 2);

% ---- H-representation of every region: A x <= b (rows normalized) ----
A = cell(1,K); b = cell(1,K);
for k = 1:K
    [A{k}, b{k}] = vert2hrep(regions{k}, tol);
end

adjFacet = false(K); adjAny = false(K);
info = repmat(struct('type','disjoint','gap',-inf,'n',[],'b',[]), K, K);
opts = optimoptions('linprog','Display','none');

for i = 1:K-1
    for j = i+1:K
        % quick reject with bounding boxes
        if any(max(regions{i},[],1) < min(regions{j},[],1) - tol) || ...
           any(max(regions{j},[],1) < min(regions{i},[],1) - tol)
            continue
        end

        % ---- signed Chebyshev radius of the intersection ----
        % max r  s.t.  a_k'x + r <= b_k  (rows normalized), r free
        Aij = [A{i}; A{j}];  bij = [b{i}; b{j}];
        [~, r] = linprog([zeros(n,1); -1], [Aij, ones(size(Aij,1),1)], bij, ...
                         [], [], [], [], opts);
        r = -r;                                  % linprog minimizes -r
        res = struct('type','disjoint','gap',r,'n',[],'b',[]);

        if r > tol
            res.type = 'overlap';
        elseif r >= -tol
            res.type = 'touch';                  % interiors disjoint, sets meet
            adjAny(i,j) = true;

            % ---- is the contact (n-1)-dimensional? look for opposing facets ----
            [fi, fj] = find(abs(A{i}*A{j}'  + 1) < 1e-6 & ...
                            abs(b{i} + b{j}')   < 1e-6*max(1,norm(b{i},inf)));
            for c = 1:numel(fi)
                nrm = A{i}(fi(c),:)';  off = b{i}(fi(c));
                P   = eye(n) - nrm*nrm';
                Ao  = [A{i}; A{j}];  bo = [b{i}; b{j}];
                pn  = vecnorm((Ao*P)', 2, 1)';   % in-plane length of each normal
                keep = pn > 1e-9;                % drop constraints parallel to nrm
                % in-plane Chebyshev ball inside the common hyperplane
                [~, rr] = linprog([zeros(n,1); -1], ...
                                  [Ao(keep,:), pn(keep)], bo(keep), ...
                                  nrm', off, [], [], opts);
                if ~isempty(rr) && -rr > tol
                    res.type = 'facet';  res.n = nrm;  res.b = off;
                    adjFacet(i,j) = true;
                    break
                end
            end
        end
        info(i,j) = res;  info(j,i) = res;
    end
end
adjFacet = adjFacet | adjFacet';
adjAny   = adjAny   | adjAny';
end

% =====================================================================
function [A, b] = vert2hrep(V, tol)
% Vertices -> halfspaces A x <= b with unit-norm rows, duplicates removed.
n  = size(V,2);
c  = mean(V,1);
K  = convhulln(V, {'Qt','Pp'});          % triangulated facets (n vertices each)
A  = zeros(size(K,1), n);  b = zeros(size(K,1),1);
for f = 1:size(K,1)
    F = V(K(f,:), :);
    N = null([F, ones(n,1)]);            % [a; -b] with a'x - b = 0 on facet
    N = N(:,1);
    a = N(1:n);  off = -N(n+1);
    s = norm(a);  a = a/s;  off = off/s;
    if a'*c' > off, a = -a; off = -off; end    % orient outward
    A(f,:) = a';  b(f) = off;
end
[~, ia] = uniquetol([A b], max(tol,1e-8), 'ByRows', true);   % merge coplanar triangles
A = A(ia,:);  b = b(ia);
end