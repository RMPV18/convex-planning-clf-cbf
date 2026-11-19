function V = polytope_vertices_2d(A, b)
%POLYTOPE_VERTICES_2D  Vertices of a bounded 2-D polytope {x : A x <= b},
%   counter-clockwise (pairwise row intersections that satisfy all rows).
V = zeros(0, 2);  m = size(A, 1);
for i = 1:m-1
    for k = i+1:m
        M2 = A([i k], :);
        if abs(det(M2)) < 1e-12, continue; end
        p = M2 \ b([i k]);
        if all(A * p - b <= 1e-9), V(end+1, :) = p'; end %#ok<AGROW>
    end
end
V = uniquetol(V, 1e-9, 'ByRows', true);
c = mean(V, 1);
[~, o] = sort(atan2(V(:, 2) - c(2), V(:, 1) - c(1)));
V = V(o, :);
end
