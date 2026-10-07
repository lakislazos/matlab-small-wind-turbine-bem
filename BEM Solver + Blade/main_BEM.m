%% main_BEM.m
% Driver script for the wind turbine blade BEM analysis.
%
% Workflow:
%   1. Choose design parameters (R, lambda, B, alpha_des).
%   2. Load aerofoil polar (CSV from AirfoilTools or similar).
%   3. Build blade geometry (chord and twist distributions).
%   4. Evaluate at the design point.
%   5. Sweep lambda to find the actual best operating point.
%   6. Check starting torque at low rotational speed.

clear; clc; close all;

%% --- 1. Design parameters --------------------------------------------
params.R          = 0.25;     % rotor tip radius (m)
params.r0         = 0.05;     % root radius (m) - measure dovetail
params.B          = 2;        % number of blades
params.lambda     = 3;        % design tip speed ratio
params.alpha_des  = 6;        % design angle of attack (deg) - from polar
params.taper      = 2.5;      % root chord / tip chord
params.N          = 20;       % number of blade elements
params.linearTwist = false;   % true => linearise twist for manufacture

U_inf_design = 10;            % design wind speed (m/s)
rho          = 1.225;         % air density (kg/m^3)

%% --- 2. Aerofoil polar -----------------------------------------------
% Download a CSV from airfoiltools.com (e.g. SG6043 at Re=100k) and
% point this at it. The reader handles the standard AirfoilTools format.
polarFile = 'sg6043_re100k.csv';
polar = loadPolar(polarFile);

fprintf('Loaded polar from %s, alpha range %.1f to %.1f deg\n', ...
        polar.file, polar.alphaMin, polar.alphaMax);

%% --- 3. Blade geometry -----------------------------------------------
geom = buildBladeGeom(params);

fprintf('\nGeometry:\n');
fprintf('  Planform area A_B = %.2f cm^2\n', geom.A_B*1e4);
fprintf('  Root chord c_R    = %.2f mm\n',   geom.c_R*1e3);
fprintf('  Tip chord  c_T    = %.2f mm\n',   geom.c_T*1e3);
fprintf('  Root twist        = %.2f deg\n',  geom.beta(1));
fprintf('  Tip twist         = %.2f deg\n',  geom.beta(end));

% Plot the geometry
figure('Name','Blade geometry','Position',[100 100 900 350]);
subplot(1,2,1);
plot(geom.r, geom.c*1e3, '-o','LineWidth',1.5); grid on;
xlabel('r (m)'); ylabel('chord c (mm)'); title('Chord distribution');
subplot(1,2,2);
plot(geom.r, geom.beta, '-o','LineWidth',1.5); grid on;
xlabel('r (m)'); ylabel('\beta (deg)'); title('Twist distribution');

%% --- 4. Evaluate at the design point ---------------------------------
Omega_design = params.lambda * U_inf_design / params.R;
opts = struct('verbose', false);

resDesign = bemSolve(geom, polar, U_inf_design, Omega_design, opts);

fprintf('\nDesign point (U=%.1f m/s, lambda=%.1f, Omega=%.1f rad/s ≈ %.0f RPM):\n', ...
        U_inf_design, params.lambda, Omega_design, Omega_design*60/(2*pi));
fprintf('  C_P   = %.4f\n', resDesign.CP);
fprintf('  C_T   = %.4f\n', resDesign.CT);
fprintf('  Power = %.3f W\n', resDesign.P);
fprintf('  Torque= %.4f Nm\n', resDesign.Q);
fprintf('  Conv. = %d / %d elements\n', sum(resDesign.converged), numel(resDesign.r));

% Plot per-element diagnostics
figure('Name','BEM design point','Position',[120 120 1000 600]);
subplot(2,2,1); plot(resDesign.r, resDesign.a, '-o', resDesign.r, resDesign.ap, '-s');
grid on; xlabel('r (m)'); ylabel('induction'); legend('a','a''');
title('Induction factors');
subplot(2,2,2); plot(resDesign.r, resDesign.alpha, '-o'); grid on;
xlabel('r (m)'); ylabel('\alpha (deg)'); title('Local angle of attack');
subplot(2,2,3); plot(resDesign.r, resDesign.CL, '-o', resDesign.r, resDesign.CD*10, '-s');
grid on; xlabel('r (m)'); ylabel('coef'); legend('C_L','10 C_D');
title('Aerofoil coefficients');
subplot(2,2,4); plot(resDesign.r, resDesign.F, '-o'); grid on;
xlabel('r (m)'); ylabel('F'); title('Prandtl loss factor');

%% --- 5. C_P-lambda sweep at design wind speed ------------------------
lambdaVec = linspace(0.5, 8, 30);
sweep = sweepLambda(geom, polar, U_inf_design, lambdaVec, opts);

[CPmax, kBest] = max(sweep.CP);
fprintf('\nC_P-lambda sweep at U=%.1f m/s:\n', U_inf_design);
fprintf('  Peak C_P = %.4f at lambda = %.2f (RPM = %.0f)\n', ...
        CPmax, sweep.lambda(kBest), sweep.Omega(kBest)*60/(2*pi));

figure('Name','C_P-lambda','Position',[140 140 700 450]);
plot(sweep.lambda, sweep.CP, '-o','LineWidth',1.5); grid on;
hold on; plot(sweep.lambda(kBest), CPmax, 'rp','MarkerSize',14,'MarkerFaceColor','r');
yline(0.593, '--', 'Betz limit');
xlabel('\lambda'); ylabel('C_P'); title('Power coefficient vs tip speed ratio');

%% --- 6. Starting torque check -----------------------------------------
% Critical: the rotor must produce more torque than dynamometer friction
% (0.03 Nm per the brief) when stationary or nearly so. We sweep at
% very low Omega (high effective alpha, mostly stalled blade).
fprintf('\nStarting torque check (U=%.1f m/s):\n', U_inf_design);
OmegaSlow = [0.1 1 5 10 20 50];   % rad/s
for k = 1:numel(OmegaSlow)
    rk = bemSolve(geom, polar, U_inf_design, OmegaSlow(k), opts);
    fprintf('  Omega=%5.1f rad/s (lambda=%.3f):  Q=%7.4f Nm\n', ...
            OmegaSlow(k), rk.lambda, rk.Q);
end
fprintf('  (must exceed 0.03 Nm at Omega -> 0 for self-start)\n');
