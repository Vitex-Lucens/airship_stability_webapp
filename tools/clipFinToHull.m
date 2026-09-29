function [clippedCloud] = clipFinToHull(finCloud_m, hullPointArray_m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CLIPFINTOHULL removes fin point cloud points that lie inside the hull.
% For each fin point, the hull radius at that x-position is interpolated
% from the hull outline. Points with radial distance from the x-axis less
% than the hull radius are removed.
%
% INPUT:
%   - finCloud_m      := (3 x N) positioned fin point cloud [m]
%   - hullPointArray_m := (2 x M) hull outline [x; r] [m]
%
% OUTPUT:
%   - clippedCloud := (3 x P) fin point cloud with interior points removed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xHull = hullPointArray_m(1, :);
rHull = hullPointArray_m(2, :);

% remove duplicate x-values (occur at section boundaries)
[xHull, uniqueIdx] = unique(xHull);
rHull = rHull(uniqueIdx);

xFin = finCloud_m(1, :);
yFin = finCloud_m(2, :);
zFin = finCloud_m(3, :);

% radial distance of each fin point from the x-axis
rFin = sqrt(yFin.^2 + zFin.^2);

% linear interpolation of hull radius at each fin x-position
% points outside the hull x-range get radius = 0 (always kept)
numFin  = length(xFin);
numHull = length(xHull);
rHullAtFin = zeros(1, numFin);
for ii = 1:numFin
    xq = xFin(ii);
    if xq <= xHull(1) || xq >= xHull(numHull)
        rHullAtFin(ii) = 0;
    else
        % find the interval containing xq
        jj = 1;
        while jj < numHull && xHull(jj+1) < xq
            jj = jj + 1;
        end
        % linear interpolation
        t = (xq - xHull(jj)) / (xHull(jj+1) - xHull(jj));
        rHullAtFin(ii) = rHull(jj) + t * (rHull(jj+1) - rHull(jj));
    end
end

% keep points that are outside the hull surface
keepMask = rFin >= rHullAtFin;
clippedCloud = finCloud_m(:, keepMask);

end
