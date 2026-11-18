run('C:\Users\Rmssp\Documents\Projects\convex-planning-clf-cbf\startup_local.m')
% addpath(genpath('C:\Program Files\Mosek\11.2\toolbox\R2019b'));
addpath(genpath('C:\Users\Rmssp\Documents\Programs Files\MATLAB\Optimization Modeling Frameworks\YALMIP-master'));
addpath(genpath(pwd));

%savepath to add permanently

% fprintf('Checking Mosek program and license\n\n')
% mosekdiag

% fprintf('Checking Yalmip\n')
% x = sdpvar(1);
% optimize(x>= 0, x,sdpsettings('debug',0), 'verbose', (double(verbose)))

% fprintf('Checking Yalmip calling Mosek\n\n')
% optimize(x>= 0, x,sdpsettings('debug',0,'solver','mosek'))