function [h, gh] = lse_cbf(x, A, b, kappa)
%LSE_CBF  Smooth CBF of a polytope {x : A x <= b} (log-sum-exp).
%   h(x) = -(1/kappa) ln sum_i exp(kappa (a_i'x - b_i)),  gh = grad h.
%   h is concave and h >= 0 implies A x <= b (under-approximation that
%   tightens as kappa grows). Same as the local copy in the controllers.
gi = A * x - b;          % a_i'x - b_i (<= 0 inside)
gmax = max(gi);
ex = exp(kappa * (gi - gmax));
S = sum(ex);
h = -(gmax + log(S) / kappa);
gh = -A.' * (ex / S);
end