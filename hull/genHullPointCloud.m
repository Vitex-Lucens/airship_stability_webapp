function [hullPointCloud] = genHullPointCloud(hullPointArray,anglesDelta_rad)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GENHULLPOINTCLOUD generates a point cloud that can be used for CAD
%programs to generate the hull skin 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lengthHullPoint = size(hullPointArray(1,:),2);
rollAngleArray_rad = 0:anglesDelta_rad:2*pi;
lengthrollAngle =  size(rollAngleArray_rad,2);
hullPointCloud = zeros(3,lengthHullPoint*lengthrollAngle);
zeroArray = zeros(1,lengthHullPoint);

for ii =1:lengthrollAngle
    DCM_nd = computeDCM(rollAngleArray_rad(ii), 0, 0);
    hullPointArrayRot = DCM_nd * [hullPointArray;zeroArray];
    startIndex = (ii-1)*lengthHullPoint+1;
    endIndex = (ii)*lengthHullPoint;
    hullPointCloud(:,(startIndex:endIndex)) = hullPointArrayRot; 
end

end

