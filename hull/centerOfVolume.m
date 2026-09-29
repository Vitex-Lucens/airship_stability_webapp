function [centerOfVolume_m] = centerOfVolume(hullPointArray)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CENTEROFVOLUME retuns the center of volume as a function of distance from
%the nose
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hullSliceArea_m2 = pi * hullPointArray(2,:).^2;
hullSectionLocation =  hullPointArray(1,:);
centerOfVolume_m = sum(hullSliceArea_m2 .* hullSectionLocation)/sum(hullSliceArea_m2); 
end

