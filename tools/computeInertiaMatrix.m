function [inertiaRotCg_kgm2, centerOfGravity_m] = computeInertiaMatrix(location_m, mass_kg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMPUTEINERTIAMATRIX computes the inertia tensor about the centre of
% gravity from a set of point masses.
%
% INPUT:
%   - location_m := (3 x N) position of each point mass [m]
%   - mass_kg    := (1 x N) mass at each point [kg]
%
% OUTPUT:
%   - inertiaRotCg_kgm2 := (3x3) inertia tensor about CG [kg.m^2]
%     NOTE: off-diagonal products of inertia use the positive convention
%           (Ixy = sum(m*x*y)), not the negated convention.
%   - centerOfGravity_m  := (3x1) centre of gravity position [m]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

centerOfGravity_m = (location_m * mass_kg') / sum(mass_kg);
locationCg_m = location_m - centerOfGravity_m;

x = locationCg_m(1, :);
y = locationCg_m(2, :);
z = locationCg_m(3, :);

Ixx = dot(mass_kg, y.^2 + z.^2);
Iyy = dot(mass_kg, x.^2 + z.^2);
Izz = dot(mass_kg, x.^2 + y.^2);
Ixy = dot(mass_kg, x .* y);
Ixz = dot(mass_kg, x .* z);
Iyz = dot(mass_kg, y .* z);

inertiaRotCg_kgm2 = [Ixx, Ixy, Ixz;
                      Ixy, Iyy, Iyz;
                      Ixz, Iyz, Izz];
end

