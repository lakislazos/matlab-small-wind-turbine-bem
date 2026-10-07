function result = bemSolve(geom, polar, U_inf, Omega, opts)
% BEMSOLVE  Full BEM with Prandtl tip+root loss and wake rotation (a').
%
%   result = bemSolve(geom, polar, U_inf, Omega)
%   result = bemSolve(geom, polar, U_inf, Omega, opts)
%
% Inputs:
%   geom    : struct with fields:
%               R      - rotor tip radius (m)
%               r0     - inboard radius where aerofoil starts (m)
%               B      - number of blades
%               r      - column vector of element radii (m), r0 < r < R
%               c      - column vector of chord at each r (m)
%               beta   - column vector of twist at each r (deg)
%   polar   : struct returned by loadPolar()
%   U_inf   : free-stream wind speed (m/s)
%   Omega   : rotational speed (rad/s)
%   opts    : (optional) struct with fields:
%               rho    - air density, default 1.225 kg/m^3
%               tol    - convergence tolerance on a, a', default 1e-5
%               maxIt  - max iterations per element, default 200
%               relax  - relaxation factor for a, a' update, default 0.3
%               useTipLoss  - default true
%               useRootLoss - default true
%               useGlauert  - Glauert high-induction correction, default true
%               verbose - print per-element diagnostics, default false
%
% Output struct fields (each is a column vector matching geom.r):
%   a       - axial induction at each element
%   ap      - tangential induction at each element
%   phi     - inflow angle (deg)
%   alpha   - local angle of attack (deg)
%   CL, CD  - looked-up coefficients
%   F       - Prandtl loss factor
%   dT      - element thrust contribution (N), per element
%   dQ      - element torque contribution (Nm), per element
%   converged - logical vector
%
% Plus scalar fields:
%   T       - total thrust (N)
%   Q       - total torque (Nm)
%   P       - power (W) = Q * Omega
%   CP      - power coefficient
%   CT      - thrust coefficient
%   lambda  - tip speed ratio Omega*R/U_inf
%   U_inf, Omega - echoed back

    % ---- option defaults ----
    if nargin < 5, opts = struct(); end
    opts = setDefault(opts, 'rho',         1.225);
    opts = setDefault(opts, 'tol',         1e-5);
    opts = setDefault(opts, 'maxIt',       200);
    opts = setDefault(opts, 'relax',       0.3);
    opts = setDefault(opts, 'useTipLoss',  true);
    opts = setDefault(opts, 'useRootLoss', true);
    opts = setDefault(opts, 'useGlauert',  true);
    opts = setDefault(opts, 'verbose',     false);

    r    = geom.r(:);
    c    = geom.c(:);
    beta = geom.beta(:);   % degrees
    B    = geom.B;
    R    = geom.R;
    r0   = geom.r0;
    rho  = opts.rho;
    N    = numel(r);

    % Pre-allocate outputs
    a      = zeros(N,1);
    ap     = zeros(N,1);
    phi    = zeros(N,1);
    alpha  = zeros(N,1);
    CL     = zeros(N,1);
    CD     = zeros(N,1);
    F      = ones(N,1);
    dT     = zeros(N,1);
    dQ     = zeros(N,1);
    conv   = false(N,1);

    % Annulus widths for integration (centred differences with endpoint
    % handling). Using trapezoidal-style widths means each element's
    % dT and dQ are already "per element including its dr", so the
    % totals are just sums.
    dr = annulusWidths(r, r0, R);

    for i = 1:N
        % Initial guess: standard a = 1/3, no swirl
        ai  = 1/3;
        api = 0;

        for it = 1:opts.maxIt
            % Local inflow angle. atan2 is robust at small numbers.
            % Note: Omega*r is the tangential blade speed, U_inf*(1-a)
            % is the axial component of the relative wind.
            phi_i = atan2(U_inf*(1 - ai), Omega*r(i)*(1 + api));

            % Local angle of attack (deg)
            alpha_i = rad2deg(phi_i) - beta(i);

            % Look up aerofoil coefficients
            CL_i = polar.CL_fn(alpha_i);
            CD_i = polar.CD_fn(alpha_i);

            % Force coefficients in the axial/tangential frame.
            % C_T (thrust-direction coeff at element) and
            % C_Q (tangential/torque-direction coeff at element).
            sinp = sin(phi_i);
            cosp = cos(phi_i);
            Cn = CL_i*cosp + CD_i*sinp;   % normal (thrust) direction
            Ct = CL_i*sinp - CD_i*cosp;   % tangential (torque) direction

            % Local solidity
            sigma_r = B*c(i)/(2*pi*r(i));

            % Prandtl tip+root loss factor
            % f_tip = (B/2) * (R - r) / (r * sin phi)
            % f_root = (B/2) * (r - r0) / (r0 * sin phi)
            % F = F_tip * F_root, each = (2/pi)*acos(exp(-f))
            % Guard against sin(phi) -> 0 at very high a or near tip
            sinp_safe = max(abs(sinp), 1e-6);

            F_tip = 1;
            if opts.useTipLoss
                f_tip = (B/2) * (R - r(i)) / (r(i) * sinp_safe);
                F_tip = (2/pi) * acos(min(1, exp(-f_tip)));
            end
            F_root = 1;
            if opts.useRootLoss
                f_root = (B/2) * (r(i) - r0) / (r0 * sinp_safe);
                F_root = (2/pi) * acos(min(1, exp(-f_root)));
            end
            F_i = max(F_tip * F_root, 1e-4);

            % Update a using axial momentum vs blade-element thrust.
            % Standard form: a/(1-a) = sigma_r*Cn / (4*F*sin^2 phi)
            % => a = 1 / (1 + 4*F*sin^2 phi / (sigma_r*Cn))
            denom_a = 4*F_i*sinp^2 / (sigma_r*Cn);
            a_new = 1 / (1 + denom_a);

            % Glauert correction for heavily-loaded rotors (a > ~0.4).
            % The simple momentum theory breaks down because the wake
            % is no longer well-defined; Glauert's empirical fix maps
            % the same C_T but allows higher a.
            if opts.useGlauert && a_new > 0.4
                % Buhl 2005 modification (smooth, monotonic):
                %   C_T = sigma_r*(1-a)^2*Cn / sin^2 phi
                % invert empirically when a>ac (ac=0.4):
                CT_local = sigma_r*(1-a_new)^2*Cn / sinp^2;
                ac = 0.4;
                % Buhl: a = (18F - 20 - 3 sqrt(CT(50-36F) + 12F(3F-4))) / (36F - 50)
                disc = CT_local*(50 - 36*F_i) + 12*F_i*(3*F_i - 4);
                if disc >= 0
                    a_new = (18*F_i - 20 - 3*sqrt(disc)) / (36*F_i - 50);
                end
                a_new = max(min(a_new, 0.95), ac);
            end

            % Update a' from angular momentum vs blade-element torque.
            % a'/(1+a') = sigma_r*Ct / (4*F*sin phi*cos phi)
            % => a' = 1 / (4*F*sin phi cos phi / (sigma_r*Ct) - 1)
            denom_ap = 4*F_i*sinp*cosp / (sigma_r*Ct);
            ap_new = 1 / (denom_ap - 1);
            % Clamp pathological values that occasionally appear
            % near stall when Ct < 0 (drag dominates lift)
            if ~isfinite(ap_new) || ap_new < -0.5
                ap_new = 0;
            end
            ap_new = min(ap_new, 1);

            % Relaxed update
            ai_old  = ai;
            api_old = api;
            ai  = (1 - opts.relax)*ai_old  + opts.relax*a_new;
            api = (1 - opts.relax)*api_old + opts.relax*ap_new;

            if abs(ai - ai_old) < opts.tol && abs(api - api_old) < opts.tol
                conv(i) = true;
                break;
            end
        end

        % Store converged state
        a(i)     = ai;
        ap(i)    = api;
        phi(i)   = rad2deg(phi_i);
        alpha(i) = alpha_i;
        CL(i)    = CL_i;
        CD(i)    = CD_i;
        F(i)     = F_i;

        % Element forces (from momentum theory, including F)
        % dT = 4*pi*r*rho*U^2*a(1-a)*F * dr
        % dQ = 4*pi*r^3*rho*U*Omega*a'(1-a)*F * dr
        dT(i) = 4*pi*r(i)*rho*U_inf^2*ai*(1-ai)*F_i * dr(i);
        dQ(i) = 4*pi*r(i)^3*rho*U_inf*Omega*api*(1-ai)*F_i * dr(i);

        if opts.verbose
            fprintf(['r=%.3f  it=%3d  a=%.4f a''=%.4f  phi=%.2f a=%.2f  ' ...
                     'CL=%.3f CD=%.4f  F=%.3f  conv=%d\n'], ...
                    r(i), it, ai, api, rad2deg(phi_i), alpha_i, ...
                    CL_i, CD_i, F_i, conv(i));
        end
    end

    % Totals
    T = sum(dT);
    Q = sum(dQ);
    P = Q * Omega;
    A = pi*R^2;
    CP = P / (0.5*rho*U_inf^3*A);
    CT = T / (0.5*rho*U_inf^2*A);

    result = struct();
    result.r         = r;
    result.a         = a;
    result.ap        = ap;
    result.phi       = phi;
    result.alpha     = alpha;
    result.CL        = CL;
    result.CD        = CD;
    result.F         = F;
    result.dT        = dT;
    result.dQ        = dQ;
    result.dr        = dr;
    result.converged = conv;
    result.T         = T;
    result.Q         = Q;
    result.P         = P;
    result.CP        = CP;
    result.CT        = CT;
    result.lambda    = Omega*R/U_inf;
    result.U_inf     = U_inf;
    result.Omega     = Omega;
end


function dr = annulusWidths(r, r0, R)
% Compute the width assigned to each element for trapezoidal-like
% integration. The total span [r0, R] is partitioned by the midpoints
% between consecutive element centres; endpoints extend to r0 and R.
    r = r(:);
    N = numel(r);
    edges = zeros(N+1, 1);
    edges(1)   = r0;
    edges(end) = R;
    for k = 2:N
        edges(k) = 0.5*(r(k-1) + r(k));
    end
    dr = diff(edges);
end


function s = setDefault(s, field, value)
    if ~isfield(s, field) || isempty(s.(field))
        s.(field) = value;
    end
end
