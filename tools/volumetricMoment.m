function [volumetricMoment_m5] = volumetricMoment(cv_m, hullPointArray_m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VOLUMETRICMOMENT Calculate the volumetric moment of inertia about the
% centre of volume. This is the second moment of the cross-sectional area
% distribution: integral of pi*r^2*(x - cv)^2 dx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_m = hullPointArray_m(1, :);
r_m = hullPointArray_m(2, :);
dx_m = [diff(x_m), x_m(end) - x_m(end-1)];
volumetricMoment_m5 = pi * sum(r_m.^2 .* (x_m - cv_m).^2 .* dx_m);
end