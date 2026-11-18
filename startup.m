% STARTUP  Path setup for this repository.
% Runs ../startup_local.m (local, not in the repo: adds third_party/, nda/,
% code/), then adds YALMIP and MOSEK. Edit the folders for another machine.
% MOSEK is appended with -end: its toolbox ships its own linprog/quadprog,
% which must not shadow the Optimization Toolbox versions used by the code;
% YALMIP still finds mosekopt.

localStartup = fullfile(fileparts(mfilename('fullpath')), '..', 'startup_local.m');
if isfile(localStartup)
    run(localStartup);
else
    warning('startup:noLocal', 'startup_local.m not found; adding code/ only.');
    addpath(genpath(fileparts(mfilename('fullpath'))));
end

yalmipDir = 'C:\Users\Rmssp\Documents\Programs Files\MATLAB\Optimization Modeling Frameworks\YALMIP-master';
if isfolder(yalmipDir)
    addpath(genpath(yalmipDir));
else
    warning('startup:noYalmip', 'YALMIP folder not found: %s', yalmipDir);
end

mosekDir = 'C:\Program Files\Mosek\11.2\toolbox\R2019b';
if isfolder(mosekDir)
    addpath(mosekDir, '-end');
else
    warning('startup:noMosek', 'MOSEK toolbox folder not found: %s', mosekDir);
end
