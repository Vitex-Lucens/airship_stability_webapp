function [finDesignStruct] = finSizeing(designParameters,hullDesignStruct)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FINSIZEING reterns to fin geometry, this is specific for the NACA 0018 fin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
finDesignStruct = struct;
dp  = designParameters;
hd  = hullDesignStruct;
momentArm_m =   dp.xPosTailEdge_m - hd.centerOfVolume_m - dp.finSaftyFactor_m;

%required area
requiredArea_m2 = hd.hullVolume_m3 * (dp.k2 - dp.k1) / (1.5 * momentArm_m) * 0.384*2;
baseHight_m     = dp.hullRadius_m .* (1 - (dp.hullLength_m - dp.noseLength_m - ...
    dp.tailRatio - dp.xPosTailEdge_m ).^2./(dp.noseLength_m + dp.tailRatio)^2).^(1/2);

finDesignStruct.requiredArea_m2 = requiredArea_m2;
finDesignStruct.baseHight_m     = baseHight_m;

switch dp.numbFins
    case 3
        finDesignStruct.baseCord_m = requiredArea_m2/((dp.tailingEdgeHight_m + baseHight_m)*(1+cos(pi/2-dp.dihedral_rad)))...
            + dp.tailingEdgeHight_m^2 * tan(dp.LeadingEdgeTapperAngle_rad)/(2*(dp.tailingEdgeHight_m+baseHight_m));
    case 4
        finDesignStruct.baseCord_m = requiredArea_m2/((dp.tailingEdgeHight_m + baseHight_m)*(2*cos(pi/2-dp.dihedral_rad)))...
            + dp.tailingEdgeHight_m^2 * tan(dp.LeadingEdgeTapperAngle_rad)/(2*(dp.tailingEdgeHight_m+baseHight_m));
   
    otherwise
        error('numbFins must be 3 or 4, got %d', dp.numbFins)
end
end