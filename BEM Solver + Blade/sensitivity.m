%% sensitivity.m
% Compare BEM results across several blade designs at a FIXED operating
% point, so geometry changes are isolated from operating-point changes.
%
% This addresses the common pitfall in main_BEM.m: when you vary lambda
% you change BOTH the blade chord (via eq. 29) AND the rotational speed
% (Omega = lambda*U/R). Those two effects partially cancel, making it
% look like the BEM solver is insensitive to geometry. It isn't —
% you just need to isolate the variables.
%
% Usage: edit the `cases` cell array below, then run.

clear; clc; close all;

%% Fixed operating point and aerofoil
U_inf  = 10;       % m/s
Omega  = 120;      % rad/s  (equivalent to lambda=3 at R=0.25, U=10)
polar  = loadPolar('sg6043_re100k.csv');

%% Base design parameters (will be overridden per case)
base.R          = 0.25;
base.r0         = 0.05;
base.B          = 3;
base.lambda     = 3;        % only used to set chord magnitude via eq. 29
base.alpha_des  = 6;
base.taper      = 2.5;
base.N          = 20;
base.linearTwist = false;

%% Cases to compare
% Each entry is a struct of overrides on the base design.
cases = { ...
    struct('name','Base',         'B',3, 'lambda',3,   'taper',2.5), ...
    struct('name','2 blades',     'B',2, 'lambda',3,   'taper',2.5), ...
    struct('name','4 blades',     'B',4, 'lambda',3,   'taper',2.5), ...
    struct('name','Low lambda',   'B',3, 'lambda',2,   'taper',2.5), ...
    struct('name','High lambda',  'B',3, 'lambda',4.5, 'taper',2.5), ...
    struct('name','Heavy taper',  'B',3, 'lambda',3,   'taper',4.0), ...
    struct('name','Light taper',  'B',3, 'lambda',3,   'taper',1.5), ...
};

%% Run each case
fprintf('%-15s %8s %8s %8s %8s %9s %8s\n', ...
        'Case','c_R(mm)','c_T(mm)','A_B(cm2)','Q(Nm)','P(W)','CP');
fprintf('%s\n', repmat('-', 1, 75));

for k = 1:numel(cases)
    p = base;
    f = fieldnames(cases{k});
    for j = 1:numel(f)
        if ~strcmp(f{j}, 'name')
            p.(f{j}) = cases{k}.(f{j});
        end
    end

    g = buildBladeGeom(p);
    r = bemSolve(g, polar, U_inf, Omega);

    fprintf('%-15s %8.1f %8.1f %8.2f %8.4f %9.3f %8.4f\n', ...
            cases{k}.name, g.c_R*1e3, g.c_T*1e3, g.A_B*1e4, ...
            r.Q, r.P, r.CP);
end

fprintf(['\nAll runs use IDENTICAL operating point ' ...
         '(U=%.1f m/s, Omega=%.1f rad/s).\n'], U_inf, Omega);
fprintf('Any change in Q/P/CP is therefore due to geometry alone.\n');
