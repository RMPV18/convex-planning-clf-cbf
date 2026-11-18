function [u, regime, feasible] = clf_cbf_qp_closed_form(LfV, LGV, gammaV, Lfh, LGh, alphah, tol, p_nom, mode)
%CLF_CBF_QP_CLOSED_FORM  Closed-form CLF-CBF QP, optionally with relaxed CLF.
%
%   [u, regime, feasible] = clf_cbf_qp_closed_form(LfV, LGV, gammaV, Lfh, LGh, alphah)
%   [u, regime, feasible] = clf_cbf_qp_closed_form(..., tol, p_nom, mode)
%
%   Hard CLF (p = inf, default):  u = argmin 1/2|u|^2
%       s.t. LfV + LGV u <= -gammaV (CLF),  Lfh + LGh u >= -alphah (CBF)
%   Relaxed CLF (9 inputs): slack weight p = p_nom (Matias & Silvestre,
%   Eq. 9-15 with p^{-1} > 0); if neither the CLF-only nor the joint regime
%   applies, the hard-CLF solution is meant to be returned with feasible = 0.
%
%   regime 1: only CLF active   regime 3: both active
%   regime 4: u = 0 (both inactive)   regime 2: only CBF active (CBF prioritised)
%
%   Variant of the local copy inside the gcs_clf_cbf controllers (which have
%   only the first 7 inputs). Not used by current code.
%   KNOWN ISSUE: any 9-input call sets mode = 1, including the internal
%   fallback call (which passes mode = 2), so that fallback recurses
%   without end when regimes 1 and 3 do not apply. tol is unused.

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

% Auxiliary quantities =========================================================
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
    feasible = 1;
elseif det ~= 0 && lambda1 >= 0 && lambda2 >= 0 
    u = -lambda1*LGV' + lambda2*LGh';
    regime = 3;
    feasible = 1;
elseif mode== 1
    [u, regime, feasible] = clf_cbf_qp_closed_form(LfV, LGV, gammaV, Lfh, LGh, alphah, tol, p_nom, 2);
    feasible = 0;
elseif Fv <= 0 && Fh >= 0
    u = zeros(size(LGV, 2), 1);
    regime = 4;
    feasible = 0;
else
    u = -norm(LGh)^(-2)*Fh*LGh';
    regime = 2;
    feasible = 0;
end
end