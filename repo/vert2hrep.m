function [A, b] = vert2hrep(V, tol)
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
end