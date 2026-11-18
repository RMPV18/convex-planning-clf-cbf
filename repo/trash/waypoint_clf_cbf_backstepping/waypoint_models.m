function model = waypoint_models(name, varargin)
%WAYPOINT_MODELS  System models for WAYPOINT_CLF_CBF_BACKSTEPPING.
%
%   model = waypoint_models('single_integrator', n)
%   model = waypoint_models('double_integrator', n)
%   model = waypoint_models('damped_second_order', n, mass, damping, stiffness)
%   model = waypoint_models('unicycle', l)
%
% MODEL STRUCT (also the template for custom models)
%   Common fields
%     type      'first_order' | 'second_order'
%     n         dimension of the controlled output y (top-level state x)
%     m         input dimension
%     output    @(s) y         point driven to the waypoints and used by the
%                              CLF/CBF (n x 1)
%     position  @(s) p         physical reference point used ONLY for the
%                              a-posteriori safety check (n x 1)
%     f0, G0    @(s) n x 1, n x q   top-level dynamics  ydot = f0(s) + G0(s) v
%                              first_order : v = u   (q = m)
%                              second_order: v = z   (q = n, G0 full row rank)
%     dyn       @(s,u) sdot    full state derivative (column)
%     rho       safety radius: the CBF is h(y) = h_poly(y) - rho
%   Additional fields for 'second_order' (Problem 2 of Matias & Silvestre)
%       xdot = f0(s) + G0(s) z,   zdot = f1(s) + G1(s) u
%     xIdx, zIdx   indices of x and z in the state vector s (y = s(xIdx))
%     f1, G1       @(s) n x 1, n x m   (G1 full row rank)
%
%   Mixed relative degree (unicycle): kinematic unicycle
%       pdot = v (cos th, sin th),  thdot = w,  u = (v, w)
%   is reduced to a first-order, state-dependent fully actuated system by
%   controlling the offset point y = p + l (cos th, sin th) (Glotfelter et
%   al.; Matias & Silvestre Sec. V-D):  ydot = R(th) diag(1, l) u.
%   The body is protected by rho = l (|p - y| = l): h_eff = h_poly(y) - l.
%   The body position p converges to a ball of radius l around the waypoint.

switch lower(name)
    case 'single_integrator'
        n = varargin{1};
        model = struct('type', 'first_order', 'n', n, 'm', n, ...
            'output', @(s) s(1:n), 'position', @(s) s(1:n), ...
            'f0', @(s) zeros(n, 1), 'G0', @(s) eye(n), ...
            'dyn', @(s, u) u(:), 'rho', 0);

    case 'double_integrator'
        n = varargin{1};
        model = struct('type', 'second_order', 'n', n, 'm', n, ...
            'output', @(s) s(1:n), 'position', @(s) s(1:n), ...
            'xIdx', 1:n, 'zIdx', n + 1:2 * n, ...
            'f0', @(s) zeros(n, 1), 'G0', @(s) eye(n), ...
            'f1', @(s) zeros(n, 1), 'G1', @(s) eye(n), ...
            'dyn', @(s, u) [s(n + 1:2 * n); u(:)], 'rho', 0);

    case 'damped_second_order'
        % m xddot = -d xdot - kx + u
        n = varargin{1};
        mass = varargin{2};
        damp = varargin{3};
        stiff = varargin{4};
        if mass <= 0
            error('waypoint_models:param', 'mass must be positive.');
        end
        f1 = @(s) -(stiff * s(1:n) + damp * s(n + 1:2 * n)) / mass;
        model = struct('type', 'second_order', 'n', n, 'm', n, ...
            'output', @(s) s(1:n), 'position', @(s) s(1:n), ...
            'xIdx', 1:n, 'zIdx', n + 1:2 * n, ...
            'f0', @(s) zeros(n, 1), 'G0', @(s) eye(n), ...
            'f1', f1, 'G1', @(s) eye(n) / mass, ...
            'dyn', @(s, u) [s(n + 1:2 * n); f1(s) + u(:) / mass], 'rho', 0);

    case 'unicycle'
        l = varargin{1};
        if l <= 0
            error('waypoint_models:param', 'Offset l must be positive.');
        end
        model = struct('type', 'first_order', 'n', 2, 'm', 2, ...
            'output', @(s) s(1:2) + l * [cos(s(3)); sin(s(3))], ...
            'position', @(s) s(1:2), ...
            'f0', @(s) zeros(2, 1), ...
            'G0', @(s) [cos(s(3)) -sin(s(3)); sin(s(3)) cos(s(3))] * diag([1, l]), ...
            'dyn', @(s, u) [u(1) * cos(s(3)); u(1) * sin(s(3)); u(2)], ...
            'rho', l);

    otherwise
        error('waypoint_models:name', 'Unknown model "%s".', name);
end
end
