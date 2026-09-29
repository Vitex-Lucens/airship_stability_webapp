clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Design Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
designParameters = designParametersCaseOne();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Generate Geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hullDesignStruct = hullTopFunction(designParameters);
finDesignStruct  = finTopFunction(designParameters, hullDesignStruct);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Mass Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% hull mass from skin surface area
[areaPointCloud] = hullPointCloudDensity(hullDesignStruct.hullPointArray_m, designParameters.hullAngleSpaceing_rad);
hullPointMass_kg = areaPointCloud * designParameters.hullMaterialMass_kgpm2;
hullTotalMass_kg = sum(hullPointMass_kg);

% fin mass from volume (scaled by fraction outside hull)
singleFinVolume_m3 = finVolume(designParameters.thickNessRatio, ...
    finDesignStruct.baseCord_m, ...
    designParameters.LeadingEdgeTapperAngle_rad, ...
    designParameters.tailingEdgeHight_m, ...
    designParameters.nacaCoeffs);
singleFinVolume_m3 = singleFinVolume_m3 * finDesignStruct.externalFraction_nd;
singleFinMass_kg   = singleFinVolume_m3 * designParameters.finMaterialDensity_kgpm3;

% treat each fin as a lumped mass at its point cloud centroid
numFins = finDesignStruct.numbFins;
finCentroids_m = zeros(3, numFins);
for ii = 1:numFins
    finCentroids_m(:, ii) = mean(finDesignStruct.finClouds{ii}, 2);
end
finMasses_kg = ones(1, numFins) * singleFinMass_kg;

% combine hull point cloud with fin lumped masses
allLocations_m = [hullDesignStruct.hullPointCloud, finCentroids_m];
allMasses_kg   = [hullPointMass_kg, finMasses_kg];
totalMass_kg   = sum(allMasses_kg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Pre-Ballast CG (needed to calculate required ballast)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, preBallastCG_m] = computeInertiaMatrix(allLocations_m, allMasses_kg);
cvX = hullDesignStruct.centerOfVolume_m;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Ballast
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate ballast mass to bring CG x-position to match CV x-position
% m_b = totalMass * (x_CG - x_CV) / (x_CV - x_b)
ballastMass_kg  = 0;
ballastPos_m    = [0; 0; 0];

switch designParameters.ballastOption
    case 'nose'
        % ballast at nose tip, on hull centreline
        ballastPos_m = [0; 0; 0];
        xBallast = ballastPos_m(1);
        if abs(cvX - xBallast) > 1e-6
            ballastMass_kg = totalMass_kg * (preBallastCG_m(1) - cvX) / (cvX - xBallast);
        end

    case 'gondola'
        % ballast at base of nose section, bottom of hull (z = +radius = down)
        ballastPos_m = [designParameters.noseLength_m; 0; designParameters.hullRadius_m];
        xBallast = ballastPos_m(1);
        if abs(cvX - xBallast) > 1e-6
            ballastMass_kg = totalMass_kg * (preBallastCG_m(1) - cvX) / (cvX - xBallast);
        end

    case 'none'
        % no ballast
        ballastMass_kg = 0;

    otherwise
        error('ballastOption must be ''none'', ''nose'', or ''gondola'', got ''%s''', ...
            designParameters.ballastOption);
end

% clamp to non-negative (if CG is already forward of CV, no ballast needed)
ballastMass_kg = max(0, ballastMass_kg);

% add ballast to the combined mass arrays
if ballastMass_kg > 0
    allLocations_m = [allLocations_m, ballastPos_m];
    allMasses_kg   = [allMasses_kg, ballastMass_kg];
    totalMass_kg   = sum(allMasses_kg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. Inertia Properties (with ballast)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rigid body inertia about CG (hull + fins + ballast)
[inertiaRotCg_kgm2, centerOfGravity_m] = computeInertiaMatrix(allLocations_m, allMasses_kg);

% virtual (apparent) mass inertia about CV
[volumetricMomentCv_m5] = volumetricMoment(hullDesignStruct.centerOfVolume_m, hullDesignStruct.hullPointArray_m);
virtualInertiaMatrixCv  = designParameters.airDensity_kgpm3 * volumetricMomentCv_m5 * designParameters.vertialInertiaCoeff_nd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 7. Buoyancy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[netLift_N, buoyancyForce_N, gasWeight_N, gasDensity_kgpm3] = ...
    computeBuoyancy(hullDesignStruct.hullVolume_m3, ...
    designParameters.airDensity_kgpm3, ...
    designParameters.liftingGasDensitySL_kgpm3, ...
    designParameters.temperature_K, ...
    designParameters.gravity_mps2);

totalWeight_N = totalMass_kg * designParameters.gravity_mps2;
excessLift_N  = netLift_N - totalWeight_N;  % positive = airship floats

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 8. Static Stability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% distance from CV to CG along x-axis (positive = CG aft of CV)
staticMargin_m = centerOfGravity_m(1) - cvX;

% static pitching moment about CV due to CG-CV offset
% buoyancy acts upward at CV, weight acts downward at CG
% moment = weight * (x_CG - x_CV), positive = nose up
staticMoment_Nm = totalWeight_N * staticMargin_m;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 9. Results Summary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
results.hullTotalMass_kg       = hullTotalMass_kg;
results.singleFinMass_kg       = singleFinMass_kg;
results.singleFinVolume_m3     = singleFinVolume_m3;
results.totalMass_kg           = totalMass_kg;
results.centerOfGravity_m      = centerOfGravity_m;
results.preBallastCG_m         = preBallastCG_m;
results.ballastMass_kg         = ballastMass_kg;
results.ballastPos_m           = ballastPos_m;
results.inertiaRotCg_kgm2     = inertiaRotCg_kgm2;
results.virtualInertiaMatrixCv = virtualInertiaMatrixCv;
results.staticMargin_m         = staticMargin_m;
results.staticMoment_Nm        = staticMoment_Nm;
results.netLift_N              = netLift_N;
results.buoyancyForce_N        = buoyancyForce_N;
results.gasWeight_N            = gasWeight_N;
results.totalWeight_N          = totalWeight_N;
results.excessLift_N           = excessLift_N;

fprintf('\n=== Airship Design Summary ===\n');
fprintf('Hull length:        %.2f m\n', designParameters.hullLength_m);
fprintf('Hull radius:        %.2f m\n', designParameters.hullRadius_m);
fprintf('Hull volume:        %.2f m^3\n', hullDesignStruct.hullVolume_m3);
fprintf('Centre of volume:   %.2f m (from nose)\n', hullDesignStruct.centerOfVolume_m);
fprintf('\n--- Atmosphere & Gas ---\n');
fprintf('Lifting gas:        %s\n', designParameters.liftingGas);
fprintf('Altitude:           %.0f m\n', designParameters.altitude_m);
fprintf('Temperature:        %.1f K (%.1f C)\n', designParameters.temperature_K, designParameters.temperature_K - 273.15);
fprintf('Air density:        %.4f kg/m^3\n', designParameters.airDensity_kgpm3);
fprintf('Gas density:        %.4f kg/m^3\n', gasDensity_kgpm3);
fprintf('\n--- Mass ---\n');
fprintf('Hull mass:          %.2f kg\n', hullTotalMass_kg);
fprintf('Fin mass (each):    %.3f kg\n', singleFinMass_kg);
fprintf('Fin mass (total):   %.3f kg  (%d fins)\n', singleFinMass_kg * numFins, numFins);
fprintf('\n--- Ballast ---\n');
fprintf('Option:             %s\n', designParameters.ballastOption);
if ballastMass_kg > 0
    fprintf('Ballast mass:       %.3f kg\n', ballastMass_kg);
    fprintf('Ballast position:   [%.2f, %.2f, %.2f] m\n', ballastPos_m);
    fprintf('Pre-ballast CG:     [%.3f, %.3f, %.3f] m\n', preBallastCG_m);
else
    fprintf('No ballast required\n');
end
fprintf('Total mass:         %.2f kg  (incl. ballast)\n', totalMass_kg);
fprintf('\n--- Buoyancy ---\n');
fprintf('Buoyancy (gross):   %.1f N\n', buoyancyForce_N);
fprintf('Gas weight:         %.1f N\n', gasWeight_N);
fprintf('Net lift:           %.1f N\n', netLift_N);
fprintf('Structure weight:   %.1f N\n', totalWeight_N);
if excessLift_N > 0
    floatStr = 'FLOATS';
else
    floatStr = 'SINKS';
end
fprintf('Excess lift:        %.1f N  (%s)\n', excessLift_N, floatStr);
fprintf('\n--- Stability ---\n');
fprintf('CG position:        [%.3f, %.3f, %.3f] m\n', centerOfGravity_m);
fprintf('Static margin:      %.3f m (CG - CV, positive = CG aft)\n', staticMargin_m);
fprintf('Static moment:      %.2f Nm (positive = nose up)\n', staticMoment_Nm);
fprintf('\n--- Inertia (about CG) ---\n');
fprintf('Ixx: %.2f  Iyy: %.2f  Izz: %.2f  [kg.m^2]\n', ...
    inertiaRotCg_kgm2(1,1), inertiaRotCg_kgm2(2,2), inertiaRotCg_kgm2(3,3));
fprintf('\n--- Fin Geometry ---\n');
fprintf('Base chord:         %.3f m\n', finDesignStruct.baseCord_m);
fprintf('Required area:      %.3f m^2\n', finDesignStruct.requiredArea_m2);
fprintf('External fraction:  %.1f%%\n', finDesignStruct.externalFraction_nd * 100);
fprintf('==============================\n\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 10. Panel Aerodynamics (Vortex Panel Method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_deg = 0;  % angle of attack for panel analysis
[hullAero, finAero] = runPanelAnalysis(designParameters, hullDesignStruct, finDesignStruct, alpha_deg);

fprintf('\n--- Panel Aerodynamics (alpha = %.1f deg) ---\n', alpha_deg);
fprintf('Hull:  Cl = %.4f   Cd = %.6f\n', hullAero.Cl, hullAero.Cd);
fprintf('Fin:   Cl = %.4f   Cd = %.6f\n', finAero.Cl, finAero.Cd);
fprintf('==============================\n\n');

results.hullAero = hullAero;
results.finAero  = finAero;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 11. Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotAirship(hullDesignStruct, finDesignStruct, designParameters, results);