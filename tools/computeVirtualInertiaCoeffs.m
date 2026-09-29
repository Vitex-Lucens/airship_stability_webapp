function [kMatrix_nd, kPrimeMatrix_nd] = computeVirtualInertiaCoeffs(hullLength_m, hullRadius_m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMPUTEVIRTUALINERTIACOEFFS computes the apparent (virtual) mass and
% inertia coefficients for a prolate ellipsoid using Lamb's k-factors.
%
% The airship hull is approximated as a prolate ellipsoid with semi-major
% axis a = hullLength/2 and semi-minor axis b = hullRadius.
%
% Translational coefficients (k):
%   k_x = alpha_0 / (2 - alpha_0)       [axial]
%   k_y = beta_0  / (2 - beta_0)        [lateral]
%   k_z = k_y                           [vertical, by axial symmetry]
%
% Rotational coefficients (k'):
%   k'_x = 0                            [roll - zero for axisymmetric body]
%   k'_y = (e^2*(beta_0 - alpha_0)) / ((2-e^2)*(2*e^2 - (2-e^2)*(beta_0-alpha_0)))  [pitch]
%   k'_z = k'_y                         [yaw, by axial symmetry]
%
% where alpha_0 and beta_0 are Lamb's integrals for a prolate ellipsoid:
%   alpha_0 = (2*(1-e^2)/e^3) * (0.5*ln((1+e)/(1-e)) - e)
%   beta_0  = 1/e^2 - ((1-e^2)/(2*e^3)) * ln((1+e)/(1-e))
%   e = sqrt(1 - (b/a)^2)   [eccentricity]
%
% INPUT:
%   - hullLength_m := total hull length [m]
%   - hullRadius_m := hull max radius [m]
%
% OUTPUT:
%   - kMatrix_nd      := (3x3) diagonal translational virtual mass coefficients
%   - kPrimeMatrix_nd := (3x3) diagonal rotational virtual inertia coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = hullLength_m / 2;  % semi-major axis
b = hullRadius_m;       % semi-minor axis

% Eccentricity of the prolate ellipsoid
e = sqrt(1 - (b/a)^2);

% Lamb's integrals
alpha_0 = (2*(1 - e^2) / e^3) * (0.5 * log((1+e)/(1-e)) - e);
beta_0  = 1/e^2 - ((1 - e^2) / (2*e^3)) * log((1+e)/(1-e));

% Translational virtual mass coefficients
k_x = alpha_0 / (2 - alpha_0);
k_y = beta_0  / (2 - beta_0);
k_z = k_y;

kMatrix_nd = [k_x, 0,   0;
              0,   k_y, 0;
              0,   0,   k_z];

% Rotational virtual inertia coefficients
k_prime_x = 0;
e2 = e^2;
k_prime_y = (e2 * (beta_0 - alpha_0)) / ...
    ((2 - e2) * (2*e2 - (2 - e2)*(beta_0 - alpha_0)));
k_prime_z = k_prime_y;

kPrimeMatrix_nd = [k_prime_x, 0,         0;
                   0,         k_prime_y,  0;
                   0,         0,          k_prime_z];

end
