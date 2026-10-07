function geom = buildBladeGeom(params)
% BUILDBLADEGEOM  Generate blade geometry (chord & twist) from design params.
%
%   geom = buildBladeGeom(params)
%
% Implements the design equations from the project brief:
%   - Planform area  A_B = R^2 / (B*lambda)              (eq. 29)
%   - Linear chord   c(r) = c_R + (c_T - c_R)*(r-r0)/(R-r0)
%       with c_R/c_T = taper, and (c_R+c_T)/2 * (R-r0) = A_B
%   - Optimal twist  beta(r) = atan(2/(3*lambda*mu)) - alpha_des   (eq. 31)
%
% params struct (all required unless marked):
%   R         - tip radius (m), e.g. 0.25
%   r0        - inboard radius where aerofoil section starts (m), e.g. 0.05
%   B         - number of blades
%   lambda    - design tip speed ratio
%   alpha_des - design angle of attack (deg)
%   taper     - c_R/c_T ratio, e.g. 2.5
%   N         - number of blade elements (default 20)
%   linearTwist - logical, if true approximate optimal twist with a
%                 straight line between root and tip values (default false)

    if ~isfield(params, 'N') || isempty(params.N)
        params.N = 20;
    end
    if ~isfield(params, 'linearTwist') || isempty(params.linearTwist)
        params.linearTwist = false;
    end

    R         = params.R;
    r0        = params.r0;
    B         = params.B;
    lambda    = params.lambda;
    alpha_des = params.alpha_des;
    taper     = params.taper;
    N         = params.N;

    % Element radii: evenly spaced between r0 and R (excluding the
    % exact endpoints to avoid singularities at r0 and R).
    r = linspace(r0, R, N+2)';
    r = r(2:end-1);

    % --- Planform area target (eq. 29) ---
    A_B = R^2 / (B*lambda);

    % --- Linear chord ---
    % (c_R + c_T)/2 * (R - r0) = A_B  and  c_R = taper*c_T
    c_T = 2*A_B / ((R - r0)*(1 + taper));
    c_R = taper * c_T;
    c = c_R + (c_T - c_R)*(r - r0)/(R - r0);

    % --- Twist (eq. 31) ---
    mu = r / R;
    phi_opt = atan(2 ./ (3*lambda*mu));      % radians
    beta_opt_deg = rad2deg(phi_opt) - alpha_des;

    if params.linearTwist
        beta = beta_opt_deg(1) + ...
               (beta_opt_deg(end) - beta_opt_deg(1))*(r - r(1))/(r(end) - r(1));
    else
        beta = beta_opt_deg;
    end

    % --- Pack geometry ---
    geom.R         = R;
    geom.r0        = r0;
    geom.B         = B;
    geom.r         = r;
    geom.c         = c;
    geom.beta      = beta;
    geom.A_B       = A_B;
    geom.c_R       = c_R;
    geom.c_T       = c_T;
    geom.lambda_d  = lambda;       % design lambda (for record)
    geom.alpha_des = alpha_des;
    
end
