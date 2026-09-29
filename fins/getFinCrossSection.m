function [finCrossSectionArray] = getFinCrossSection(thickNessRatio, cordLength_m, spaceing_m, nacaCoeffs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FINCROSSSECTION defines fin cross section for a symetric NACA fin. The fin
% is built forward with the tailing edge representing the (0,0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xArray = -cordLength_m:spaceing_m:0;

k1 = nacaCoeffs(1);
k2 = nacaCoeffs(2);
k3 = nacaCoeffs(3);
k4 = nacaCoeffs(4);
k5 = nacaCoeffs(5);
var1 = ((xArray+cordLength_m)/cordLength_m);
yArray = (thickNessRatio * cordLength_m)/0.2 * (k1 * var1.^(1/2) + k2 * var1...
    + k3 * var1.^2 + k4 * var1.^3 + k5 * var1.^4);

finCrossSectionArray = [[xArray;yArray],flip([xArray;-yArray],2)];
end