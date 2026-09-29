function [finVolume_m3] = finVolume(thickNessRatio, baseCord_m, ...
    LeadingEdgeTapperAngle_rad, tailingEdgeHight_m, nacaCoeffs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE THIS FUNCTION IS FULLY AI GENERATED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FINVOLUME computes the volume of a single fin by integrating the NACA
% cross-sectional area along the span (height). All fins are assumed
% identical so this only needs to be called once.
%
% The fin tapers linearly from baseCord_m at the root to zero chord at the
% tip. At each span station the cross-section is a symmetric NACA profile
% whose area can be computed analytically.
%
% For a symmetric NACA 4-digit aerofoil with thickness ratio t:
%   A(c) = 2 * integral_0^c y(x) dx
% where y(x) = (t*c/0.2)*(k1*sqrt(x/c) + k2*(x/c) + k3*(x/c)^2
%                          + k4*(x/c)^3 + k5*(x/c)^4)
% Integrating analytically (substituting u = x/c) gives:
%   A(c) = 2 * (t * c^2 / 0.2) * (k1*2/3 + k2/2 + k3/3 + k4/4 + k5/5)
%
% INPUT:
%   - thickNessRatio              := NACA thickness ratio (e.g. 0.15)
%   - baseCord_m                  := chord length at the fin root [m]
%   - LeadingEdgeTapperAngle_rad  := leading-edge taper angle [rad]
%   - tailingEdgeHight_m          := fin span (height from root to tip) [m]
%   - nacaCoeffs                  := (1x5) NACA 4-digit symmetric profile coefficients [k1..k5]
%
% OUTPUT:
%   - finVolume_m3 := volume of a single fin [m^3]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NACA 4-digit symmetric profile coefficients
k1 = nacaCoeffs(1);
k2 = nacaCoeffs(2);
k3 = nacaCoeffs(3);
k4 = nacaCoeffs(4);
k5 = nacaCoeffs(5);

% Analytical cross-section area coefficient: A(c) = areaCoeff * c^2
% A = 2 * integral_0^c y(x) dx, substituting u = x/c, dx = c*du:
% A = 2 * (t*c^2/0.2) * integral_0^1 (k1*u^0.5 + k2*u + k3*u^2 + k4*u^3 + k5*u^4) du
%   = 2 * (t*c^2/0.2) * (k1*2/3 + k2/2 + k3/3 + k4/4 + k5/5)
areaCoeff = (thickNessRatio / 0.2) * 2 * (k1*2/3 + k2/2 + k3/3 + k4/4 + k5/5);

% The fin has a full triangular planform height
fullHight_m = tan(pi/2 - LeadingEdgeTapperAngle_rad) * baseCord_m;

% Analytical integration along the span.
% A(h) = areaCoeff * tan(alpha)^2 * (H - h)^2, so:
% V = areaCoeff * tan(alpha)^2 * integral_0^span (H - h)^2 dh
%   = areaCoeff * tan(alpha)^2 * [(H-h)^3 / -3] from 0 to span
%   = areaCoeff * tan(alpha)^2 / 3 * (H^3 - (H - span)^3)
tanAlpha = tan(LeadingEdgeTapperAngle_rad);
finVolume_m3 = areaCoeff * tanAlpha^2 / 3 * ...
    (fullHight_m^3 - (fullHight_m - tailingEdgeHight_m)^3);

end

