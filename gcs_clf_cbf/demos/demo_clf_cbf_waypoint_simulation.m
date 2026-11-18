%% DEMO_CLF_CBF_WAYPOINT_SIMULATION  CLF-CBF waypoint following on box regions.
%   Tracking stage of the GCS + CLF-CBF pipeline, run on hand-made regions.
%   Edit the SELECTION block only.
%     MODEL     'unicycle'          -> unicycle_waypoint_clf_cbf (first order)
%               'double_integrator' -> double_integrator_waypoint_clf_cbf (second order)
%     SCENARIO  'convex'            three convex regions
%               'union'             region 2 is a union of two boxes (L-shape)
%   CHECK_QP compares the closed-form QP with quadprog (Optimization Toolbox) on
%   random instances and verifies that the QP copies in the controller files are
%   identical to the copy in this script.
%   PARAMS overrides any field of WP_PARAMS (local function at the end of each
%   controller file). Each controller file is self-contained.

clear; clc;

% ----------------------------- SELECTION ---------------------------------
MODEL = 'double_integrator';          % 'unicycle' | 'double_integrator'
SCENARIO = 'union';          % 'convex' | 'union'
UNICYCLE_OFFSET = 0.3;       % l (m), 'unicycle' only
X0 = [1; 1];                 % initial position
HEADING0 = 0;                % initial heading (rad), 'unicycle' only
V0 = [0; 0];                 % initial velocity, 'double_integrator' only
PARAMS = struct();           % e.g. struct('kappa', 30, 'alpha_bar', 3)
PLOT = true;
CHECK_QP = false;
% -------------------------------------------------------------------------

[boxes, W] = scenario_data(SCENARIO);
regions = boxes_to_regions(boxes);

if CHECK_QP
    check_closed_form_vs_quadprog();
    check_qp_copies_identical();
end

switch MODEL
    case 'unicycle'
        out = unicycle_waypoint_clf_cbf(W, regions, [X0; HEADING0], UNICYCLE_OFFSET, PARAMS);
    case 'double_integrator'
        out = double_integrator_waypoint_clf_cbf(W, regions, X0, V0, PARAMS);
    otherwise
        error('demo_clf_cbf_waypoint_simulation:model', 'Unknown MODEL "%s".', MODEL);
end

fprintf('\n=== %s / %s ===\njumps [t from to (beta_h)]:\n', MODEL, SCENARIO);
disp(out.jumps);
fn = fieldnames(out.verification);
for i = 1:numel(fn)
    fprintf('%-22s = %g\n', fn{i}, out.verification.(fn{i}));
end
fprintf('expected jumps         = %d\n', size(W, 2) - 1);

if PLOT
    plot_results(out, boxes, W, X0, MODEL, SCENARIO);
end

% ============================= LOCAL FUNCTIONS ============================
function [boxes, W] = scenario_data(name)
%SCENARIO_DATA  boxes{k}: rows [xmin xmax ymin ymax] (1 row convex, 2 rows union).
switch name
    case 'convex'
        boxes = {[0 7 0 2], [5 10 0 8], [2 10 3 9]};
        W = [6 1; 8 5; 3 7.5].';
    case 'union'
        boxes = {[0 7 0 2], [5 7.5 0 8; 5 10 3 8], [2 10 3 9]};
        W = [6 1; 8 5; 3 7.5].';
    otherwise
        error('demo_clf_cbf_waypoint_simulation:scenario', 'Unknown SCENARIO "%s".', name);
end
end

function regions = boxes_to_regions(boxes)
M = numel(boxes);
regions = struct('A', cell(1, M), 'b', cell(1, M));
for k = 1:M
    nb = size(boxes{k}, 1);
    Ac = cell(1, nb);
    bc = cell(1, nb);
    for i = 1:nb
        bx = boxes{k}(i, :);
        Ac{i} = [1 0; -1 0; 0 1; 0 -1];
        bc{i} = [bx(2); -bx(1); bx(4); -bx(3)];
    end
    if nb == 1
        regions(k).A = Ac{1};
        regions(k).b = bc{1};
    else
        regions(k).A = Ac;
        regions(k).b = bc;
    end
end
end

function plot_results(out, boxes, W, X0, MODEL, SCENARIO)
M = size(W, 2);
colors = lines(M);
figure('Name', sprintf('%s / %s', MODEL, SCENARIO));
tiledlayout(1, 2);
nexttile;
hold on; axis equal; grid on;
for k = 1:M
    for i = 1:size(boxes{k}, 1)
        b = boxes{k}(i, :);
        patch([b(1) b(2) b(2) b(1)], [b(3) b(3) b(4) b(4)], colors(k, :), ...
            'FaceAlpha', 0.15, 'EdgeColor', colors(k, :));
    end
    idx = out.mode == k;
    plot(out.p(1, idx), out.p(2, idx), '-', 'Color', colors(k, :), 'LineWidth', 2);
end
plot(W(1, :), W(2, :), 'kp', 'MarkerFaceColor', 'y', 'MarkerSize', 12);
plot(X0(1), X0(2), 'ko', 'MarkerFaceColor', 'k');
for j = 1:size(out.jumps, 1)
    [~, i] = min(abs(out.t - out.jumps(j, 1)));
    plot(out.p(1, i), out.p(2, i), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
end
xlabel('x_1'); ylabel('x_2');
title('Body trajectory, regions (colour = active mode), jumps (x)');

nexttile;
hold on; grid on;
plot(out.t, out.h, '--', 'LineWidth', 1.5);
names = {'h_{eff}'};
if isfield(out, 'h1')
    plot(out.t, out.h1, 'LineWidth', 1.5);
    names{end + 1} = 'h_1';
end
plot(out.t, out.bodyMargin, ':', 'LineWidth', 1.5);
names{end + 1} = 'body margin';
yline(0, 'k-');
xlabel('t');
legend(names, 'Location', 'best');
title('Safety signals');
end

function check_closed_form_vs_quadprog()
%CHECK_CLOSED_FORM_VS_QUADPROG  Random-instance comparison with quadprog.
if ~license('test', 'Optimization_Toolbox')
    warning('demo_clf_cbf_waypoint_simulation:noOptTbx', 'Optimization Toolbox unavailable: QP check skipped.');
    return
end
rng(1);
q = 2;
nT = 2000;
worst = 0;
nInfeasible = 0;
reg = zeros(1, 4);
opts = optimoptions('quadprog', 'Display', 'off', 'OptimalityTolerance', 1e-12, 'ConstraintTolerance', 1e-12);
for i = 1:nT
    LfV = randn;
    LGV = randn(1, q);
    gam = abs(randn);
    Lfh = randn;
    LGh = randn(1, q);
    alp = abs(randn);
    [u, rg, feas] = clf_cbf_qp_closed_form(LfV, LGV, gam, Lfh, LGh, alp);
    [uq, ~, flag] = quadprog(eye(q), zeros(q, 1), [LGV; -LGh], [-(LfV + gam); Lfh + alp], [], [], [], [], [], opts);
    if flag <= 0 || ~feas
        nInfeasible = nInfeasible + 1;
        continue
    end
    reg(rg) = reg(rg) + 1;
    worst = max(worst, norm(u - uq));
end
fprintf('Closed-form QP vs quadprog: %d feasible instances, max |u - u_qp| = %.3e, regimes S1..S4 = %s, skipped = %d\n', ...
    nT - nInfeasible, worst, mat2str(reg), nInfeasible);
if worst > 1e-6
    warning('demo_clf_cbf_waypoint_simulation:qpMismatch', 'Closed-form QP differs from quadprog by %.3e.', worst);
end
end

function check_qp_copies_identical()
%CHECK_QP_COPIES_IDENTICAL  The closed-form QP is a local function in each
%   controller file and in this script. Compare the source text (whitespace
%   normalised) so that the copy tested above is the copy the controllers use.
header = ['function [u, regime, ', 'feasible] = clf_cbf_qp_closed_form('];
ref = extract_function_text([mfilename('fullpath') '.m'], header);
names = {'unicycle_waypoint_clf_cbf', 'double_integrator_waypoint_clf_cbf'};
for i = 1:numel(names)
    f = which(names{i});
    if isempty(f)
        warning('demo_clf_cbf_waypoint_simulation:notOnPath', '%s is not on the MATLAB path: copy check skipped.', names{i});
        continue
    end
    if strcmp(extract_function_text(f, header), ref)
        fprintf('QP copy in %s: identical.\n', names{i});
    else
        warning('demo_clf_cbf_waypoint_simulation:qpCopyDiffers', 'QP copy in %s differs from the one tested here.', names{i});
    end
end
end

function txt = extract_function_text(file, header)
%EXTRACT_FUNCTION_TEXT  Source of the function starting at "header" up to the next
%   function or end of file, whitespace normalised.
src = fileread(file);
i0 = strfind(src, header);
if isempty(i0)
    error('demo_clf_cbf_waypoint_simulation:extract', 'Function header not found in %s.', file);
end
rest = src(i0(1):end);
nxt = regexp(rest(2:end), [char(10) 'function '], 'once');
if ~isempty(nxt)
    rest = rest(1:nxt);
end
txt = strtrim(regexprep(rest, '\s+', ' '));
end

function [u, regime, feasible] = clf_cbf_qp_closed_form(LfV, LGV, gammaV, Lfh, LGh, alphah, tol)
%CLF_CBF_QP_CLOSED_FORM  Exact solution of the hard-constrained CLF-CBF QP
%
%       u = argmin 1/2 |u|^2
%           s.t.  LfV + LGV u <= -gammaV        (CLF)
%                 Lfh + LGh u >= -alphah        (CBF)
%
%   for a single CLF/CBF pair (Matias & Silvestre, Eq. 60 = Eq. 9 with p = inf;
%   closed form Eq. 10-15 with p^{-1} = 0).
%
%   [u, regime, feasible] = clf_cbf_qp_closed_form(LfV, LGV, gammaV, Lfh, LGh, alphah)
%
% INPUTS
%   LfV, Lfh   scalars  Lie derivatives along the drift f
%   LGV, LGh   1 x q    Lie derivatives along G (row vectors)
%   gammaV     scalar   gamma(V(x))
%   alphah     scalar   alpha(h(x))
%   tol        scalar   optional numerical tolerance (default 1e-12)
%
% OUTPUTS
%   u         q x 1   minimiser
%   regime    1..4    active set (Eq. 14): 1 CLF only (S1), 2 CBF only (S2),
%                     3 both (S3), 4 none (S4, u = 0)
%   feasible  logical false if the two constraints are incompatible
%                     (aligned gradients, Eq. 48-53) or the CLF gradient
%                     vanishes while the CLF constraint is violated. In that
%                     case the CBF-only solution (safety priority) is returned.
%
% With FV = LfV + gammaV, Fh = Lfh + alphah, L = LGV LGh',
% Delta = L^2 - |LGV|^2 |LGh|^2:
%   S4: FV <= 0 and Fh >= 0                         u = 0
%   S1: u = -FV LGV' / |LGV|^2        (CBF satisfied)
%   S2: u = -Fh LGh' / |LGh|^2        (CLF satisfied)
%   S3: u = -l1 LGV' + l2 LGh',
%       l1 = (L Fh - |LGh|^2 FV)/Delta, l2 = (|LGV|^2 Fh - L FV)/Delta, l1,l2 >= 0.
% The candidates are tested in the order S4, S1, S2, S3, which yields the
% unique KKT point of the strictly convex QP.

if nargin < 7 || isempty(tol)
    tol = 1e-12;
end
a = LGV(:);
g = LGh(:);
q = numel(a);
if numel(g) ~= q
    error('clf_cbf_qp:dim', 'LGV and LGh must have the same length.');
end
FV = LfV + gammaV;
Fh = Lfh + alphah;
na = a.' * a;
ng = g.' * g;
feasible = true;

% S4: both constraints inactive at u = 0
if FV <= 0 && Fh >= 0
    u = zeros(q, 1);
    regime = 4;
    return
end

u1 = [];
u2 = [];
% S1: CLF only
if na > tol^2 && FV > 0
    u1 = -(FV / na) * a;
    if g.' * u1 >= -Fh - 1e-12
        u = u1;
        regime = 1;
        return
    end
end
% S2: CBF only
if ng > tol^2 && Fh < 0
    u2 = -(Fh / ng) * g;
    if a.' * u2 <= -FV + 1e-12
        u = u2;
        regime = 2;
        return
    end
end
% S3: both active
L = a.' * g;
Delta = L^2 - na * ng;
if abs(Delta) > 1e-12 * max(1, na * ng)
    l1 = (L * Fh - ng * FV) / Delta;
    l2 = (na * Fh - L * FV) / Delta;
    if l1 >= -1e-9 && l2 >= -1e-9
        u = -l1 * a + l2 * g;
        regime = 3;
        return
    end
end

% incompatible constraints: keep safety, report infeasibility
feasible = false;
regime = 2;
if isempty(u2)
    if ng > tol^2 && Fh < 0
        u2 = -(Fh / ng) * g;
    else
        u2 = zeros(q, 1);
    end
end
u = u2;
end
