function results = demo_waypoint_clf_cbf()
%DEMO_WAYPOINT_CLF_CBF  Same corridor and waypoints, three plant models:
%   (1) double integrator, (2) damped second-order system, (3) kinematic unicycle.
%   Also verifies the closed-form CLF-CBF QP against quadprog (if available).

boxes = [0 7 0 2; 5 10 0 8; 2 10 3 9];
W = [6 1; 8 5; 3 7.5].';
M = size(W, 2);
regions = struct('A', cell(1, M), 'b', cell(1, M));
for k = 1:M
    [regions(k).A, regions(k).b] = box_to_halfspaces(boxes(k, :));
end

check_closed_form_vs_quadprog();

cases = {
    'double integrator',   waypoint_models('double_integrator', 2),                [1; 1; 0; 0]
    'damped second order', waypoint_models('damped_second_order', 2, 2, 0.8, 0.5), [1; 1; 0.5; 0]
    'unicycle (l = 0.3)',  waypoint_models('unicycle', 0.3),                       [1; 1; 0]
    };
nC = size(cases, 1);
results = cell(1, nC);
colors = lines(M);
figure('Name', 'Waypoint CLF-CBF backstepping');
tiledlayout(2, nC);
for c = 1:nC
    fprintf('\n=== %s ===\n', cases{c, 1});
    out = waypoint_clf_cbf_backstepping(W, regions, cases{c, 2}, cases{c, 3});
    results{c} = out;
    ver = out.verification;
    fprintf('jumps [t from to beta_h]:\n');
    disp(out.jumps);
    fprintf('min h1 / h_eff          = %.4e\n', ver.minH1);
    fprintf('min body margin         = %.4e\n', ver.minBodyMargin);
    fprintf('max CBF violation       = %.4e\n', ver.maxCbfViolation);
    fprintf('max CLF violation       = %.4e\n', ver.maxClfViolation);
    fprintf('max backstep violation  = %.4e\n', ver.maxBackstepViolation);
    fprintf('QP fallback samples     = %d\n', ver.nQpFallback);
    fprintf('final output error      = %.4e\n', ver.finalOutputError);
    fprintf('jumps                   = %d (expected %d)\n', ver.nJumps, M - 1);

    nexttile(c);
    hold on; axis equal; grid on;
    for k = 1:M
        b = boxes(k, :);
        patch([b(1) b(2) b(2) b(1)], [b(3) b(3) b(4) b(4)], colors(k, :), ...
            'FaceAlpha', 0.15, 'EdgeColor', colors(k, :));
    end
    P = zeros(2, numel(out.t));
    for i = 1:numel(out.t)
        P(:, i) = cases{c, 2}.position(out.s(:, i));
    end
    for k = 1:M
        idx = out.mode == k;
        plot(P(1, idx), P(2, idx), '-', 'Color', colors(k, :), 'LineWidth', 2);
    end
    plot(W(1, :), W(2, :), 'kp', 'MarkerFaceColor', 'y', 'MarkerSize', 12);
    for j = 1:size(out.jumps, 1)
        [~, i] = min(abs(out.t - out.jumps(j, 1)));
        plot(P(1, i), P(2, i), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    end
    title(cases{c, 1});
    xlabel('x_1'); ylabel('x_2');

    nexttile(nC + c);
    hold on; grid on;
    plot(out.t, out.h1, 'LineWidth', 1.5);
    plot(out.t, out.h, '--', 'LineWidth', 1.5);
    plot(out.t, out.bodyMargin, ':', 'LineWidth', 1.5);
    yline(0, 'k-');
    xlabel('t');
    legend('h_1', 'h_{eff}', 'body margin', 'Location', 'best');
end
end

% =========================================================================
function check_closed_form_vs_quadprog()
%CHECK_CLOSED_FORM_VS_QUADPROG  Random-instance comparison with quadprog.
if ~license('test', 'Optimization_Toolbox')
    warning('demo:noOptTbx', 'Optimization Toolbox unavailable: closed-form QP not compared with quadprog.');
    return
end
rng(1);
q = 2;
nT = 2000;
worst = 0;
nInfeasible = 0;
opts = optimoptions('quadprog', 'Display', 'off', 'OptimalityTolerance', 1e-12, 'ConstraintTolerance', 1e-12);
reg = zeros(1, 4);
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
fprintf('Closed-form QP vs quadprog: %d feasible instances, max |u - u_qp| = %.3e, regimes S1..S4 = %s, skipped (infeasible) = %d\n', ...
    nT - nInfeasible, worst, mat2str(reg), nInfeasible);
if worst > 1e-6
    warning('demo:qpMismatch', 'Closed-form QP differs from quadprog by %.3e.', worst);
end
end

function [A, b] = box_to_halfspaces(bx)
%BOX_TO_HALFSPACES  [xmin xmax ymin ymax] -> A x <= b.
A = [1 0; -1 0; 0 1; 0 -1];
b = [bx(2); -bx(1); bx(4); -bx(3)];
end
