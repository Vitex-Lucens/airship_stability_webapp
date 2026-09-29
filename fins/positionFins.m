function [finPointCloudFinal_m] = positionFins(finPointCloud_m, xPosTailEdge_m, zPosBase_m, rotation_rad)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%POSITIONFINS positions and orients a fin.
% 1. Offset the fin radially so its base sits at the hull surface (z)
% 2. Rotate around the x-axis for dihedral placement
% 3. Translate along x to the tail position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% offset z so fin base is at hull radius
cloud_m = finPointCloud_m;
cloud_m(3,:) = cloud_m(3,:) + zPosBase_m;

% rotate around x-axis (roll) to place fin at correct angular position
DCM_nd = computeDCM(rotation_rad, 0, 0);
cloud_m = DCM_nd * cloud_m;

% translate to tail x-position
cloud_m(1,:) = cloud_m(1,:) + xPosTailEdge_m;

finPointCloudFinal_m = cloud_m;
end

