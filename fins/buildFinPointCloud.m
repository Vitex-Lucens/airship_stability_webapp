function [finDesignStruct] = buildFinPointCloud(designParameters, finDesignStruct)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BUILDFINPOINTCLOUD build the fin as a 3D point cloud
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LeadingEdgeTapperAngle_rad = designParameters.LeadingEdgeTapperAngle_rad;
tailingEdgeHight_m         = designParameters.tailingEdgeHight_m;
thickNessRatio             = designParameters.thickNessRatio;
nacaCoeffs                 = designParameters.nacaCoeffs;
spaceing_m                 = designParameters.finPointSpaceing_m;
baseCord_m                 = finDesignStruct.baseCord_m;
numRibs                    = designParameters.numRibs;

ribSpaceing_m     = tailingEdgeHight_m/(numRibs-1);
hightArray_m      = 0:ribSpaceing_m:tailingEdgeHight_m;
fullHight_m       = tan(pi/2-LeadingEdgeTapperAngle_rad)*baseCord_m;

finPointCloud_m = [];
for ii =1:1:length(hightArray_m)

    distanceDownFin_m   = hightArray_m(ii);
    distanceRemaining_m = fullHight_m - distanceDownFin_m;
    baseCordLength_m    = tan(LeadingEdgeTapperAngle_rad) * distanceRemaining_m;
    
    finCrossSectionArray_m = getFinCrossSection(thickNessRatio, baseCordLength_m, spaceing_m, nacaCoeffs);  
    zArray_m = ones(1,size(finCrossSectionArray_m,2)) * distanceDownFin_m;
    
    finCrossSectionArray_m = [finCrossSectionArray_m;zArray_m];
    finPointCloud_m        = [finPointCloud_m,finCrossSectionArray_m];
end

finDesignStruct.finPointCloud_m = finPointCloud_m;
end