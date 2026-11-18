function [u, regime, feasible] = clf_cbf_qp_closed_form(LfV, LGV, gammaV, Lfh, LGh, alphah, tol, p_nom, mode)

p = inf;
mode  = 2;

if nargin < 7 || isempty(tol)
    tol = 1e-12;
end
if nargin > 8 
    mode= 1;
    p = p_nom;
end

Fv = LfV + gammaV;
Fh = Lfh + alphah;

% Auxiliar Variables ===========================================================
L = LGV*LGh';
s1 = (1/p + norm(LGV)^2)*Fh - L*Fv;
s2 = norm(LGh)^2*Fv - L*Fh;
det = L^2 - (1/p + norm(LGV)^2)*norm(LGh)^2;
lambda1 = det^(-1)*(L*Fh - norm(LGh)^2*Fv);
lambda2 = det^(-1)*((1/p + norm(LGV)^2)*Fh - L*Fv);

% Case Selection ===============================================================
if Fv >= 0 && s1 > 0 
    u = -(1/p + norm(LGV)^2)^(-1)*Fv*LGV';
    regime = 1;
elseif det ~= 0 && lambda1 >= 0 && lambda2 >= 0 
    u = -lambda1*LGV' + lambda2*LGh';
     regime = 3;
elseif mode== 1
    clf_cbf_qp_closed_form(LfV, LGV, gammaV, Lfh, LGh, alphah, tol, p_nom, 2)
elseif Fv <= 0 && Fh >= 0
    u = zeros(size(LGV, 2), 1);
     regime = 4;
else
    u = -norm(LGh)^(-2)*Fh*LGh';
     regime = 2;
end
feasible = 1;
end