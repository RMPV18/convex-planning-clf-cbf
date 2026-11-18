function out = pid_trajectory_tracker(plan, params, obstacles)
%PID_TRAJECTORY_TRACKER  Tracks a dcp_polygon_planner trajectory with a discrete PID; plant via ode45.
%   out = pid_trajectory_tracker(plan)
%   out = pid_trajectory_tracker(plan, params)
%   out = pid_trajectory_tracker(plan, params, obstacles)
%
%   plan      : struct from dcp_polygon_planner (fields x, t, v)
%   params    : optional struct, missing fields take the defaults below
%   obstacles : optional cell of Nx2 convex polygons -> out.clearance
%       m 1 | c 0.5 | umax 10 | Kp 25 | Ki 2 | Kd 10 | ff true | tSettle 2 | p0Off [0;0]
%       relTol 1e-8 | absTol 1e-10 (ode45)
%
%   Plant (replace plant_dyn with your robot): 2-D point mass, m*a = u - c*v, |u| <= umax.
%   Controller, every Ts with zero-order hold:
%       u = Kp*e + Ki*int(e) + Kd*(v_ref - v) + m*a_ref + c*v_ref,   e = p_ref - p
%   Derivative acts on the known reference velocity (no differentiation of measurements).
%   Anti-windup: conditional integration (integrator frozen while the force is saturated).
%
%   out fields: t, p, v, u, pr, err (tracking error norm), clearance (signed distance to the
%   nearest obstacle per sample, only when obstacles is given).

if nargin < 2, params = struct(); end
params = fill_defaults(params, struct('m', 1.0, 'c', 0.5, 'umax', 10, 'Kp', 25, 'Ki', 2, 'Kd', 10, ...
    'ff', true, 'tSettle', 2.0, 'p0Off', [0; 0], 'relTol', 1e-8, 'absTol', 1e-10));

Ts = plan.t(2) - plan.t(1);
N  = numel(plan.t) - 1;
K  = N + round(params.tSettle / Ts);

pr = [plan.x, repmat(plan.x(:,end), 1, K - N)];      % 2 x (K+1) position reference
vr = [plan.v, zeros(2, K - N)];                      % 2 x (K+1) velocity reference
ar = [diff(vr, 1, 2) / Ts, zeros(2, 1)];             % acceleration reference

p = zeros(2, K+1); v = zeros(2, K+1); u = zeros(2, K);
p(:,1) = plan.x(:,1) + params.p0Off(:);
I = zeros(2, 1);
opts = odeset('RelTol', params.relTol, 'AbsTol', params.absTol);

for k = 1:K
    e  = pr(:,k) - p(:,k);
    ed = vr(:,k) - v(:,k);
    uu = params.Kp*e + params.Ki*I + params.Kd*ed;
    if params.ff, uu = uu + params.m*ar(:,k) + params.c*vr(:,k); end

    nu = norm(uu);
    if nu > params.umax
        uu = uu * params.umax / nu;
    else
        I = I + e * Ts;                              % conditional integration (anti-windup)
    end
    u(:,k) = uu;

    [~, Y] = ode45(@(~, s) plant_dyn(s, uu, params), [0 Ts], [p(:,k); v(:,k)], opts);
    p(:,k+1) = Y(end, 1:2)';
    v(:,k+1) = Y(end, 3:4)';
end

assert(all(isfinite(p(:))), 'Simulation diverged.');
out.t = (0:K) * Ts;
out.p = p; out.v = v; out.u = u; out.pr = pr;
out.err = sqrt(sum((pr - p).^2, 1));
if nargin >= 3 && ~isempty(obstacles)
    d = inf(1, K+1);
    for j = 1:numel(obstacles), d = min(d, poly_dist(p, obstacles{j})); end
    out.clearance = d;
end
end

%% ------------------------------------------------------------------------------------------
function params = fill_defaults(params, def)
f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(params, f{i}), params.(f{i}) = def.(f{i}); end
end
end

function ds = plant_dyn(s, u, params)
ds = [s(3:4); (u - params.c * s(3:4)) / params.m];
end

function d = poly_dist(p, P)
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
