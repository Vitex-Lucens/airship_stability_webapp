function [designParameters] = designParametersCaseOne()
%DESIGNPARAMETERSHL40 gives the design parameters for the airship 40 

designParameters = struct;
% hull
designParameters.hullLength_m = 12.75;
designParameters.noseLength_m = 2;
designParameters.tailRatio = 3.5;
designParameters.hullRadius_m = 1.22;
designParameters.hullMaterialMass_kgpm2 = 0.12;

% atmospheric and gas properties
designParameters.altitude_m    = 0;         % flight altitude [m]
designParameters.temperature_K = 288.15;    % ambient temperature [K] (ISA sea level = 288.15)
designParameters.liftingGas    = 'helium';  % 'helium' or 'hydrogen'
designParameters.gravity_mps2  = 9.81;

% sea level gas densities at 288.15 K [kg/m^3]
switch designParameters.liftingGas
    case 'helium'
        designParameters.liftingGasDensitySL_kgpm3 = 0.1786;
    case 'hydrogen'
        designParameters.liftingGasDensitySL_kgpm3 = 0.0899;
    otherwise
        error('liftingGas must be ''helium'' or ''hydrogen'', got ''%s''', designParameters.liftingGas);
end

%hull virtual mass and inertia coefficients (Lamb's k-factors)
[kMatrix_nd, kPrimeMatrix_nd] = computeVirtualInertiaCoeffs(designParameters.hullLength_m, designParameters.hullRadius_m);
designParameters.k1 = kMatrix_nd(1,1);       % axial virtual mass coeff
designParameters.k2 = kMatrix_nd(2,2);       % lateral virtual mass coeff
designParameters.vertialInertiaCoeff_nd = kPrimeMatrix_nd;  % rotational virtual inertia coeff matrix
designParameters.kPrime = kPrimeMatrix_nd(2,2);             % pitch/yaw rotational coeff
%designParameters.sideSlipAngle_rad = deg2rad(20);
designParameters.forwardAirSpeed_mps = 10;

% compute air density from altitude and temperature
[airDensity_kgpm3, ~] = atmosphereModel(designParameters.altitude_m, designParameters.temperature_K);
designParameters.airDensity_kgpm3 = airDensity_kgpm3;

%fins
designParameters.LeadingEdgeTapperAngle_rad = deg2rad(15);
designParameters.tailingEdgeHight_m = 1;
designParameters.thickNessRatio = 0.15;
designParameters.nacaCoeffs = [0.29690, -0.12600, -0.35160, 0.28430, -0.10150]; % NACA 4-digit symmetric profile coefficients
designParameters.xPosTailEdge_m = 11.75;
designParameters.dihedral_rad = deg2rad(50);
designParameters.numbFins = 3;
designParameters.finSaftyFactor_m = 0.25;
designParameters.finMaterialDensity_kgpm3 = 3; % fin material volumetric density [kg/m^3]
designParameters.numRibs = 100;

% ballast options: 'none', 'nose', or 'gondola'
% nose   = ballast at x = 0 (nose tip), on hull centreline
% gondola = ballast at x = noseLength, at bottom of hull (z = +radius)
designParameters.ballastOption = 'gondola';

%Motors
designParameters.motorDihdral_rad    = deg2rad(20);
designParameters.motorSafterFactor_m = 0.1;
designParameters.motorMountLength_m  = 0.75;
designParameters.motorMountWidth_m   = 0.75;

% spacing data 
designParameters.hullPointSpaceing_m   = 0.1/6;
designParameters.hullAngleSpaceing_rad = 0.01;
designParameters.finPointSpaceing_m    = 0.012;
end