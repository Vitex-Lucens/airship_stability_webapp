function [hullAero, finAero] = runPanelAnalysis(designParameters, hullData, finData, alpha_deg)
%RUNPANELANALYSIS  Run vortex panel method on hull profile and fin airfoil.
%
%   Inputs:
%     designParameters - design parameter struct
%     hullData         - hull geometry from hullTopFunction
%     finData          - fin geometry from finTopFunction
%     alpha_deg        - angle of attack [deg]
%
%   Outputs:
%     hullAero - struct with hull panel results (Cl, Cd, Cp, xc, etc.)
%     finAero  - struct with fin panel results
%
%   The hull is treated as a 2D body in the x-r plane (axisymmetric
%   profile). The upper surface is the hull outline [x, r], the lower
%   surface is [x, -r], forming a closed shape.
%
%   The fin uses the NACA symmetric section from getFinCrossSection.

Vinf = designParameters.forwardAirSpeed_mps;

%% Hull profile as a closed 2D airfoil
hx = hullData.hullPointArray(1, :);
hr = hullData.hullPointArray(2, :);

% Build closed contour: upper surface (x, +r), then lower reversed (x, -r)
% Convention: start at trailing edge (x = hullLength, r = 0),
% go along upper surface to nose, then along lower surface back to TE.
% The hull outline starts at x=0 (nose) and ends at x=hullLength (tail).
% For VPM we need: TE -> lower -> LE -> upper -> TE

% Upper surface: nose to tail = x increasing, r positive
% Lower surface: nose to tail = x increasing, r negative
% VPM ordering: TE -> lower(tail to nose) -> LE -> upper(nose to tail) -> TE
hx_upper = hx;           % nose to tail, r > 0
hr_upper = hr;
hx_lower = hx;           % nose to tail, r < 0
hr_lower = -hr;

% Assemble: TE -> lower (reversed, tail to nose) -> upper (nose to tail) -> TE
Xb_hull = [fliplr(hx_lower), hx_upper(2:end)]';
Yb_hull = [fliplr(hr_lower), hr_upper(2:end)]';

c_hull = designParameters.hullLength_m;

[h_cl, h_cd, h_Cp, h_V, h_xc, h_yc, h_s] = vortexPanelMethod(Xb_hull, Yb_hull, Vinf, c_hull, alpha_deg);

hullAero.Cl       = h_cl;
hullAero.Cd       = h_cd;
hullAero.Cp       = h_Cp;
hullAero.V        = h_V;
hullAero.xc       = h_xc;
hullAero.yc       = h_yc;
hullAero.panelS   = h_s;
hullAero.Xb       = Xb_hull;
hullAero.Yb       = Yb_hull;
hullAero.chord    = c_hull;

%% Fin airfoil cross-section
% Generate fin cross-section (NACA symmetric, built forward from TE at origin)
finCS = getFinCrossSection(designParameters.thickNessRatio, ...
    finData.baseCord_m, designParameters.finPointSpaceing_m, ...
    designParameters.nacaCoeffs);

% finCS is [2 x N] with x, y.  x goes from -chord to 0 (upper), then 0 to -chord (lower)
% Shift so that x ranges from 0 to chord for VPM convention
Xb_fin = (finCS(1, :) + finData.baseCord_m)';
Yb_fin = finCS(2, :)';

c_fin = finData.baseCord_m;

[f_cl, f_cd, f_Cp, f_V, f_xc, f_yc, f_s] = vortexPanelMethod(Xb_fin, Yb_fin, Vinf, c_fin, alpha_deg);

finAero.Cl       = f_cl;
finAero.Cd       = f_cd;
finAero.Cp       = f_Cp;
finAero.V        = f_V;
finAero.xc       = f_xc;
finAero.yc       = f_yc;
finAero.panelS   = f_s;
finAero.Xb       = Xb_fin;
finAero.Yb       = Yb_fin;
finAero.chord    = c_fin;

end
