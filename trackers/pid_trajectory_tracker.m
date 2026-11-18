function out = pid_trajectory_tracker(plan, par, obs)
%PID_TRAJECTORY_TRACKER  Tracks a dcp_polygon_planner trajectory with a discrete PID; plant via ode45.
%   out = pid_trajectory_tracker(plan)
%   out = pid_trajectory_tracker(plan, par)
%   out = pid_trajectory_tracker(plan, par, obs)   % obs: cell of Nx2 convex polygons -> out.clearance
%
%   plan : struct from dcp_polygon_planner (fields x, t, v)
%   par  : optional struct, missing fields take the defaults below
%       m 1 | c 0.5 | umax 10 | Kp 25 | Ki 2 | Kd 10 | ff true | tSettle 2 | p0Off [0;0]
%       relTol 1e-8 | absTol 1e-10 (ode45)
%
%   Plant (replace plantDyn with your robot): 2-D point mass, m*a = u - c*v, |u| <= umax.
%   Controller, every Ts with zero-order hold:
%       u = Kp*e + Ki*int(e) + Kd*(v_ref - v) + m*a_ref + c*v_ref,   e = p_ref - p
%   Derivative acts on the known reference velocity (no differentiation of measurements).
%   Anti-windup: conditional integration (integrator frozen while the force is saturated).
%
%   out fields: t, p, v, u, pr, err (tracking error norm), clearance (signed distance to the
%   nearest obstacle per sample, only when obs is given).

if nargin < 2, par = struct(); end
par = fillDefaults(par, struct('m', 1.0, 'c', 0.5, 'umax', 10, 'Kp', 25, 'Ki', 2, 'Kd', 10, ...
    'ff', true, 'tSettle', 2.0, 'p0Off', [0; 0], 'relTol', 1e-8, 'absTol', 1e-10));

Ts = plan.t(2) - plan.t(1);
N  = numel(plan.t) - 1;
K  = N + round(par.tSettle / Ts);

pr = [plan.x, repmat(plan.x(:,end), 1, K - N)];      % 2 x (K+1) position reference
vr = [plan.v, zeros(2, K - N)];                      % 2 x (K+1) velocity reference
ar = [diff(vr, 1, 2) / Ts, zeros(2, 1)];             % acceleration reference

p = zeros(2, K+1); v = zeros(2, K+1); u = zeros(2, K);
p(:,1) = plan.x(:,1) + par.p0Off(:);
I = zeros(2, 1);
opts = odeset('RelTol', par.relTol, 'AbsTol', par.absTol);

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

    [~, Y] = ode45(@(~, s) plantDyn(s, uu, par), [0 Ts], [p(:,k); v(:,k)], opts);
    p(:,k+1) = Y(end, 1:2)';
    v(:,k+1) = Y(end, 3:4)';
end

assert(all(isfinite(p(:))), 'Simulation diverged.');
out.t = (0:K) * Ts;
out.p = p; out.v = v; out.u = u; out.pr = pr;
out.err = sqrt(sum((pr - p).^2, 1));
if nargin >= 3 && ~isempty(obs)
    d = inf(1, K+1);
    for j = 1:numel(obs), d = min(d, polyDist(p, obs{j})); end
    out.clearance = d;
end
end

%% ------------------------------------------------------------------------------------------
function par = fillDefaults(par, def)
f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(par, f{i}), par.(f{i}) = def.(f{i}); end
end
end

function ds = plantDyn(s, u, par)
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
