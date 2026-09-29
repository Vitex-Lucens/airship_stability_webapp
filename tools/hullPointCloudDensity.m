function [areaPointCloud] = hullPointCloudDensity(hullPointArray_m, anglesDelta_rad)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HULLPOINTCLOUDDENSITY returns the surface area represented by each point
% in the hull point cloud. The hull profile is swept around the x-axis in
% angular increments of anglesDelta_rad.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx_m = hullPointArray_m(1, 2:end) - hullPointArray_m(1, 1:end-1);
dy_m = hullPointArray_m(2, 2:end) - hullPointArray_m(2, 1:end-1);
dx_m = [dx_m, dx_m(end)];
dy_m = [dy_m, dy_m(end)];

% Arc length along the profile curve per segment
theta_rad = atan2(dy_m, dx_m);
h_m = dx_m ./ cos(theta_rad);

% Circumferential arc length per angular slice
s_m = hullPointArray_m(2, :) * anglesDelta_rad;

% Surface area per point
areaArray_m2 = h_m .* s_m;

% Replicate for each angular slice
lengthHullPoint    = length(areaArray_m2);
rollAngleArray_rad = 0:anglesDelta_rad:2*pi;
lengthRollAngle    = length(rollAngleArray_rad);
areaPointCloud     = zeros(1, lengthHullPoint * lengthRollAngle);

for ii = 1:lengthRollAngle
    startIndex = (ii - 1) * lengthHullPoint + 1;
    endIndex   = ii * lengthHullPoint;
    areaPointCloud(startIndex:endIndex) = areaArray_m2;
end

end