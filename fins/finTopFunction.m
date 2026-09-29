function [finDesignStruct] = finTopFunction(designParameters, hullDesignStruct)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FINTOPFUNCTION runs the fin design functions.
% Supports 3-fin (inverted-Y) and 4-fin (X-tail) configurations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

finDesignStruct = finSizeing(designParameters, hullDesignStruct);
finDesignStruct = buildFinPointCloud(designParameters, finDesignStruct);

zPosBase_m       = finDesignStruct.baseHight_m;
xPosTailEdge_m   = designParameters.xPosTailEdge_m;
finPointCloud_m  = finDesignStruct.finPointCloud_m;
dihedral_rad     = designParameters.dihedral_rad;
hullPointArray_m = hullDesignStruct.hullPointArray_m;
numPointsUnclipped = size(finPointCloud_m, 2);

% Coordinate convention: x forward, y starboard, z down
% Fin is built with span in +z. Roll around x-axis places fins:
%   roll = 0    -> fin points down (+z)
%   roll = pi   -> fin points up   (-z)
% Dihedral angle measured from the downward vertical.
switch designParameters.numbFins
    case 3
        % inverted-Y: one fin up, two lower fins splayed by dihedral
        % roll=0 -> span in +z (down), roll=pi -> span in -z (up)
        finAngles_rad = [pi, ...                         % top fin (up, -z)
                         -dihedral_rad, ...              % bottom right
                          dihedral_rad];                 % bottom left
    case 4
        % X-tail: two upper, two lower, symmetric about vertical
        finAngles_rad = [pi + dihedral_rad, ...          % upper right
                         pi - dihedral_rad, ...          % upper left
                         -dihedral_rad, ...              % lower right
                          dihedral_rad];                 % lower left
    otherwise
        error('numbFins must be 3 or 4, got %d', designParameters.numbFins);
end

% position, clip, and store each fin
finDesignStruct.finClouds = cell(1, length(finAngles_rad));
for ii = 1:length(finAngles_rad)
    rawCloud = positionFins(finPointCloud_m, xPosTailEdge_m, zPosBase_m, finAngles_rad(ii));
    finDesignStruct.finClouds{ii} = clipFinToHull(rawCloud, hullPointArray_m);
end
finDesignStruct.numbFins = length(finAngles_rad);

% fraction of fin outside the hull (all fins identical, use first)
numPointsClipped = size(finDesignStruct.finClouds{1}, 2);
finDesignStruct.externalFraction_nd = numPointsClipped / numPointsUnclipped;

end