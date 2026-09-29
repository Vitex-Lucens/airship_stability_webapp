function [netLift_N, buoyancyForce_N, gasWeight_N, gasDensity_kgpm3] = ...
    computeBuoyancy(hullVolume_m3, airDensity_kgpm3, gasDensitySL_kgpm3, ...
    temperature_K, gravity_mps2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMPUTEBUOYANCY computes the net buoyant lift of the airship.
%
% Buoyancy = weight of displaced air - weight of lifting gas inside hull.
% Gas density scales with pressure and temperature from sea level reference
% using the ideal gas law (density proportional to P/T).
%
% INPUT:
%   - hullVolume_m3       := hull volume [m^3]
%   - airDensity_kgpm3    := air density at altitude [kg/m^3]
%   - gasDensitySL_kgpm3  := lifting gas density at sea level, 288.15 K [kg/m^3]
%   - temperature_K       := ambient temperature [K]
%   - gravity_mps2        := gravitational acceleration [m/s^2]
%
% OUTPUT:
%   - netLift_N        := net upward lift force [N] (positive = upward)
%   - buoyancyForce_N  := gross buoyancy (weight of displaced air) [N]
%   - gasWeight_N      := weight of lifting gas inside hull [N]
%   - gasDensity_kgpm3 := gas density at altitude [kg/m^3]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reference conditions for gas density
T0 = 288.15;    % sea level standard temperature [K]
P0 = 101325;    % sea level standard pressure [Pa]
R_air = 287.058; % specific gas constant for dry air [J/(kg.K)]

% current pressure from air density and temperature (ideal gas)
pressure_Pa = airDensity_kgpm3 * R_air * temperature_K;

% scale gas density from sea level to current conditions
% rho_gas = rho_gas_SL * (P/P0) * (T0/T)
gasDensity_kgpm3 = gasDensitySL_kgpm3 * (pressure_Pa / P0) * (T0 / temperature_K);

% forces
buoyancyForce_N = airDensity_kgpm3 * hullVolume_m3 * gravity_mps2;
gasWeight_N     = gasDensity_kgpm3 * hullVolume_m3 * gravity_mps2;
netLift_N       = buoyancyForce_N - gasWeight_N;

end
