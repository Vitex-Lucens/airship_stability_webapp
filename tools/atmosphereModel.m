function [airDensity_kgpm3, pressure_Pa] = atmosphereModel(altitude_m, temperature_K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ATMOSPHEREMODEL computes air density and pressure using the ideal gas law
% with a standard lapse rate for the troposphere (below 11 km).
%
% Uses the barometric formula with the actual temperature rather than ISA
% temperature, allowing the user to explore temperature effects.
%
% INPUT:
%   - altitude_m    := flight altitude above sea level [m]
%   - temperature_K := ambient temperature at altitude [K]
%
% OUTPUT:
%   - airDensity_kgpm3 := air density [kg/m^3]
%   - pressure_Pa      := atmospheric pressure [Pa]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% constants
P0   = 101325;    % sea level standard pressure [Pa]
T0   = 288.15;    % sea level standard temperature [K]
L    = 0.0065;    % temperature lapse rate [K/m]
g    = 9.81;      % gravitational acceleration [m/s^2]
R    = 287.058;   % specific gas constant for dry air [J/(kg.K)]

% pressure at altitude using barometric formula (ISA lapse rate)
% P = P0 * (1 - L*h/T0)^(g/(R*L))
pressure_Pa = P0 * (1 - L * altitude_m / T0)^(g / (R * L));

% density from ideal gas law using actual temperature
airDensity_kgpm3 = pressure_Pa / (R * temperature_K);

end
