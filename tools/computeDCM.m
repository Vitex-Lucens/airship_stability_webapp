function DCM_nd = computeDCM(roll_rad, pitch_rad, yaw_rad)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the DCM from roll --> pitch --> yaw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   - roll_rad :=(1x1 aircraft roll angle in radians)
%
% OUTPUT:=( 3x3 rotation matrix going from global to body fixed coordinate
% frame )
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matrix_roll= [1 0 0 ; 0 cos(roll_rad) sin(roll_rad);0 -sin(roll_rad) cos(roll_rad)];
matrix_pitch= [cos(pitch_rad) 0 -sin(pitch_rad) ; 0 1 0; sin(pitch_rad) 0 cos(pitch_rad)];
matrix_yaw= [cos(yaw_rad) sin(yaw_rad) 0 ; -sin(yaw_rad) cos(yaw_rad) 0;0 0 1];
%DIRECTION Z -> Y -> X 
DCM_nd = matrix_yaw * matrix_pitch * matrix_roll ;
end