%% PID_TRAJECTORY_TRACKER  Tracks the dcp_polygon_planner reference with a discrete PID.
% Run dcp_polygon_planner.m first (it writes plan.mat). Local functions are at the end.
%
% Plant (replace with your robot): 2-D point mass, m*a = u - c*v, force saturation |u| <= umax.
% Controller (every Ts, zero-order hold), per axis vector form:
%   u = Kp*e + Ki*int(e) + Kd*(v_ref - v) + m*a_ref + c*v_ref
%   e = p_ref - p;  derivative acts on the known reference velocity, so no noisy differentiation.
%   Anti-windup: conditional integration (integrator frozen while the force is saturated).

clear; clc; close all;
assert(isfile('plan.mat'), 'plan.mat not found. Run dcp_polygon_planner.m first.');
S = load('plan.mat');
plan = S.plan; obs = S.obs; planPar = S.par;

%% Parameters
par.m       = 1.0;     % mass
par.c       = 0.5;     % viscous damping
par.umax    = 10;      % force-norm saturation
par.Kp      = 25;
par.Ki      = 2;
par.Kd      = 10;
par.ff      = true;    % model-based feedforward from reference velocity/acceleration
par.nSub    = 10;      % RK4 substeps per control period
par.tSettle = 2.0;     % hold final reference after the plan ends [s]
par.p0Off   = [0; 0];  % initial position offset from the reference start

%% Simulate
out = simulatePid(plan, par);

%% Metrics and checks
clr = min(cell2mat(cellfun(@(V) polyDist(out.p, V), obs', 'UniformOutput', false)), [], 1);
fprintf('RMS error %.4f | max error %.4f | final error %.4f\n', ...
    sqrt(mean(out.err.^2)), max(out.err), out.err(end));
fprintf('min clearance to obstacles %.4f (planner margin %.4f)\n', min(clr), planPar.margin);
if min(clr) <= 0
    warning('pidTracker:collision', 'Tracked trajectory enters an obstacle.');
elseif min(clr) < planPar.margin
    warning('pidTracker:margin', 'Tracked trajectory uses %.1f%% of the planner clearance margin.', ...
        100 * (planPar.margin - min(clr)) / planPar.margin);
end
assert(all(isfinite(out.p(:))), 'Simulation diverged.');

%% Plots
figure; hold on; axis equal; grid on;
for j = 1:numel(obs), fill(obs{j}(:,1), obs{j}(:,2), [0.85 0.85 0.85]); end
plot(out.pr(1,:), out.pr(2,:), 'k--', 'LineWidth', 1.2);
plot(out.p(1,:),  out.p(2,:),  'b-',  'LineWidth', 1.5);
legend('', '', '', 'reference', 'tracked'); title('PID tracking');

figure;
subplot(2,1,1); plot(out.t, out.err); grid on; ylabel('|e|');
subplot(2,1,2); plot(out.t(1:end-1), out.u'); grid on; ylabel('u'); xlabel('t [s]');

%% ------------------------------------------------------------------------------------------
function out = simulatePid(plan, par)
Ts = plan.t(2) - plan.t(1);
N  = numel(plan.t) - 1;
K  = N + round(par.tSettle / Ts);

pr = [plan.x, repmat(plan.x(:,end), 1, K - N)];      % 2 x (K+1) position reference
vr = [plan.v, zeros(2, K - N)];                      % 2 x (K+1) velocity reference
ar = [diff(vr, 1, 2) / Ts, zeros(2, 1)];             % acceleration reference

p = zeros(2, K+1); v = zeros(2, K+1); u = zeros(2, K);
p(:,1) = plan.x(:,1) + par.p0Off;
I = zeros(2, 1);
h = Ts / par.nSub;

for k = 1:K
    e  = pr(:,k) - p(:,k);
    ed = vr(:,k) - v(:,k);
    uu = par.Kp*e + par.Ki*I + par.Kd*ed;
    if par.ff, uu = uu + par.m*ar(:,k) + par.c*vr(:,k); end

    nu = norm(uu);
    if nu > par.umax
        uu = uu * par.umax / nu;
    else
        I = I + e * Ts;                              % conditional integration (anti-windup)
    end
    u(:,k) = uu;

    s = [p(:,k); v(:,k)];
    for q = 1:par.nSub                               % RK4, input held constant over Ts
        k1 = plant(s, uu, par);
        k2 = plant(s + h/2*k1, uu, par);
        k3 = plant(s + h/2*k2, uu, par);
        k4 = plant(s + h*k3, uu, par);
        s = s + h/6 * (k1 + 2*k2 + 2*k3 + k4);
    end
    p(:,k+1) = s(1:2); v(:,k+1) = s(3:4);
end

out.t = (0:K) * Ts;
out.p = p; out.v = v; out.u = u; out.pr = pr;
out.err = sqrt(sum((pr - p).^2, 1));
end

function ds = plant(s, u, par)
ds = [s(3:4); (u - par.c * s(3:4)) / par.m];
end

function d = polyDist(p, P)
% Signed distance (negative inside) from points p (2xM) to the convex hull of vertices P (Nx2).
k = convhull(P(:,1), P(:,2));
V = P(k(1:end-1), :);
n = size(V, 1);
d = inf(1, size(p, 2));
for i = 1:n
    a = V(i,:)'; b = V(mod(i, n) + 1, :)'; ab = b - a;
    t = min(max(((p - a)' * ab) / (ab' * ab), 0), 1)';
    q = a + ab * t;
    d = min(d, sqrt(sum((p - q).^2, 1)));
end
in = inpolygon(p(1,:), p(2,:), V(:,1), V(:,2));
d(in) = -d(in);
end
