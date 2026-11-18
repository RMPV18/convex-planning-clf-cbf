function out = track_gcs_pid(traj, gains, regions, x0)
%TRACK_GCS_PID  Track a GcsMotionPlanner trajectory with a joint-wise PID.
%
%   out = track_gcs_pid(traj, gains, regions, x0)
%
%   traj     info.trajectory.sample from GcsMotionPlanner.solve:
%            t (N x 1), q (N x n), qdot (N x n)
%   gains    struct Kp, Ki, Kd (scalar or n-vector)
%   regions  optional cell of ConvexSet (the planner's regions); [] to skip
%   x0       optional initial state [q; qdot] (2n x 1); default = reference
%
%   Plant (assumed, fully actuated double integrator): qddot = u
%   Control: u = Kp*e + Ki*int(e) + Kd*(qdot_ref - qdot),  e = q_ref - q
%   Stable iff Kd*Kp > Ki (characteristic s^3 + Kd s^2 + Kp s + Ki).
%
%   out: t, q, qdot, u, e, rms_error, max_error, max_region_violation
%        (<= 0: executed samples lie in the union of the regions).
%
%   shortcut: no actuator limits, delay or noise; containment is checked at
%   the trajectory samples only. Upgrade when targeting a real robot.
%   shortcut: no acceleration feedforward; tracking error is transient after
%   acceleration jumps (eta = 0 trajectories).
%
%   track_gcs_pid() with no arguments runs a self-check.

if nargin == 0, selfCheck(); return; end
if nargin < 3, regions = []; end
t = traj.t(:);  n = size(traj.q, 2);
assert(numel(t) > 2 && all(diff(t) > 0), 'track_gcs_pid:time', 'traj.t must be strictly increasing.');
Kp = gains.Kp(:) .* ones(n, 1);  Ki = gains.Ki(:) .* ones(n, 1);  Kd = gains.Kd(:) .* ones(n, 1);
assert(all(Kd .* Kp > Ki), 'track_gcs_pid:gains', 'Unstable gains: need Kd*Kp > Ki.');
if nargin < 4 || isempty(x0), x0 = [traj.q(1, :)'; traj.qdot(1, :)']; end

f = @(tt, z) dyn(tt, z, traj, Kp, Ki, Kd, n);
[~, Z] = ode45(f, t, [x0(:); zeros(n, 1)], odeset('RelTol', 1e-8, 'AbsTol', 1e-10));

q = Z(:, 1:n);  v = Z(:, n+1:2*n);  I = Z(:, 2*n+1:end);
e = traj.q - q;
u = e .* Kp' + I .* Ki' + (traj.qdot - v) .* Kd';
out = struct('t', t, 'q', q, 'qdot', v, 'u', u, 'e', e);
err = sqrt(sum(e.^2, 2));
out.rms_error = sqrt(mean(err.^2));
out.max_error = max(err);
out.max_region_violation = NaN;
if ~isempty(regions)
    viol = zeros(numel(t), 1);
    for k = 1:numel(t)
        viol(k) = min(cellfun(@(cs) max(cs.A * q(k, :)' - cs.b), regions));
    end
    out.max_region_violation = max(viol);
end
end

function dz = dyn(tt, z, traj, Kp, Ki, Kd, n)
q = z(1:n);  v = z(n+1:2*n);  I = z(2*n+1:end);
qr = interp1(traj.t, traj.q, tt, 'pchip')';
vr = interp1(traj.t, traj.qdot, tt, 'pchip')';
e = qr - q;
dz = [v; Kp .* e + Ki .* I + Kd .* (vr - v); e];
end

function selfCheck()
t = linspace(0, 30, 301)';
traj = struct('t', t, 'q', [t, 0.5*t], 'qdot', repmat([1 0.5], numel(t), 1));
g = struct('Kp', 25, 'Ki', 5, 'Kd', 10);
box = struct('A', [eye(2); -eye(2)], 'b', [40; 40; 1; 1]);   % contains the line
L = track_gcs_pid(traj, g, {box}, [0.2; -0.1; 1; 0.5]);
assert(norm(L.e(end, :)) < 1e-2, 'Self-check: error did not converge.');
assert(L.max_error < 0.3, 'Self-check: transient too large.');
assert(L.max_region_violation <= 0, 'Self-check: containment report wrong.');
disp('track_gcs_pid self-check passed.');
end
