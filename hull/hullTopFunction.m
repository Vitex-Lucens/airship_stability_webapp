function [hullDesignStruct] = hullTopFunction(designParameters)
%HULLTOPFUNCTION runs the hull design scripts

[hullPointArray_m] = hullOutline(designParameters);
[hullVolume_m3]    = hullVolume(hullPointArray_m);
[centerOfVolume_m] = centerOfVolume(hullPointArray_m);
[hullPointCloud]   = genHullPointCloud(hullPointArray_m,designParameters.hullAngleSpaceing_rad);

hullDesignStruct.hullPointArray_m = hullPointArray_m;
hullDesignStruct.hullVolume_m3    = hullVolume_m3;
hullDesignStruct.centerOfVolume_m = centerOfVolume_m;
hullDesignStruct.hullPointCloud   = hullPointCloud;
end