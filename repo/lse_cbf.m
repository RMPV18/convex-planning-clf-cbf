function [h, gh] = lse_cbf(x, A, b, kappa)
%LSE_CBF  h(x) = -(1/kappa) ln sum_i exp(kappa (a_i'x - b_i)); concave; h>=0 => Ax<=b.
gi = A * x - b; % gi = -hi = a_i'x - b_i <= 0 
gmax = max(gi);
ex = exp(kappa * (gi - gmax));
S = sum(ex);
h = -(gmax + log(S) / kappa);
gh = -A.' * (ex / S);
end