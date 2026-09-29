function [hullPointArray] = hullOutline(designParameters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HULLOUTLINE generates the hull outline as a 2×N array [x; r] defining
% the axisymmetric profile. Three sections:
%   Nose:   elliptical cap from x=0 to x=noseLength
%   Centre: constant radius cylinder
%   Tail:   elliptical taper to a point at x=hullLength
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spaceing_m   = designParameters.hullPointSpaceing_m;
noseLength_m = designParameters.noseLength_m;
length_m     = designParameters.hullLength_m;
tailRatio    = designParameters.tailRatio;
radius_m     = designParameters.hullRadius_m;

% Section boundaries
tailSemiAxis_m   = noseLength_m + tailRatio;       % tail ellipse semi-axis along x
centreEnd_m      = length_m - tailSemiAxis_m;       % where centre section ends / tail begins

% x-coordinate arrays for each section
xNose   = 0              : spaceing_m : noseLength_m;
xCentre = noseLength_m   : spaceing_m : centreEnd_m;
xTail   = centreEnd_m    : spaceing_m : length_m;

% Radius profiles
rNose   = radius_m .* sqrt(1 - ((noseLength_m - xNose).^2) ./ noseLength_m^2);
rCentre = ones(1, length(xCentre)) .* radius_m;
rTail   = radius_m .* sqrt(1 - ((centreEnd_m - xTail).^2) ./ tailSemiAxis_m^2);

hullPointArray = [[xNose, xCentre, xTail]; [rNose, rCentre, rTail]];
end