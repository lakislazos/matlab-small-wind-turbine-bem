function res = sweepLambda(geom, polar, U_inf, lambdaVec, opts)
% SWEEPLAMBDA  Compute CP, CT, Q across a range of tip speed ratios.
%
%   res = sweepLambda(geom, polar, U_inf, lambdaVec)
%
% At each lambda the rotational speed is set to Omega = lambda*U_inf/R,
% holding U_inf fixed. Useful for finding the operating point that
% maximises C_P for this blade shape.
%
% Output struct fields (each is a column vector of length numel(lambdaVec)):
%   lambda, Omega, CP, CT, Q, P, T

    if nargin < 5, opts = struct(); end
    R = geom.R;
    N = numel(lambdaVec);
    res.lambda = lambdaVec(:);
    res.Omega  = zeros(N,1);
    res.CP     = zeros(N,1);
    res.CT     = zeros(N,1);
    res.Q      = zeros(N,1);
    res.P      = zeros(N,1);
    res.T      = zeros(N,1);

    for k = 1:N
        Omega = lambdaVec(k)*U_inf/R;
        r = bemSolve(geom, polar, U_inf, Omega, opts);
        res.Omega(k) = Omega;
        res.CP(k)    = r.CP;
        res.CT(k)    = r.CT;
        res.Q(k)     = r.Q;
        res.P(k)     = r.P;
        res.T(k)     = r.T;
    end
end
