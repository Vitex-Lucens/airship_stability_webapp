function [hullVolume_m3] = hullVolume(hullPointArray)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HULLVOLUME computes the hull volume by disk integration: V = sum(pi*r^2*dx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx_m = diff(hullPointArray(1, :));
r_m  = hullPointArray(2, 2:end);
hullVolume_m3 = pi * sum(r_m.^2 .* dx_m);
end