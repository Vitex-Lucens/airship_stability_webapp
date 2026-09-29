function plotAirship(hullDesignStruct, finDesignStruct, designParameters, results)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOTAIRSHIP engineering visualisation of the airship design.
% Four-panel layout: 3D view, side profile, rear cross-section, data panel.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hd = hullDesignStruct;
fd = finDesignStruct;
dp = designParameters;
r  = results;

% colours
hullCol = [0.6, 0.85, 0.9];
finCol  = [0.2, 0.4, 0.8];
cvCol   = [0.0, 0.7, 0.0];
cgCol   = [0.9, 0.1, 0.1];

hullR   = hd.hullPointArray_m(2, :);
maxR    = max(hullR);

fig = figure('Name', 'Airship Design Tool', 'NumberTitle', 'off', ...
    'Position', [50, 50, 1200, 800], 'Color', 'w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PANEL 1: 3D View (top-left)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2, 2, 1);
hold on

% hull
hc = hd.hullPointCloud;
plot3(hc(1,:), hc(2,:), hc(3,:), '.', 'MarkerSize', 1, 'Color', hullCol)

% fins
for ii = 1:fd.numbFins
    fc = fd.finClouds{ii};
    plot3(fc(1,:), fc(2,:), fc(3,:), '.', 'MarkerSize', 1, 'Color', finCol)
end

% CV and CG as vertical lines spanning the hull diameter
cvX = hd.centerOfVolume_m;
cgX = r.centerOfGravity_m(1);
plot3([cvX, cvX], [-maxR, maxR], [0, 0], '-', 'Color', cvCol, 'LineWidth', 2)
plot3([cvX, cvX], [0, 0], [-maxR, maxR], '-', 'Color', cvCol, 'LineWidth', 2)
plot3([cgX, cgX], [-maxR, maxR], [0, 0], '--', 'Color', cgCol, 'LineWidth', 2)
plot3([cgX, cgX], [0, 0], [-maxR, maxR], '--', 'Color', cgCol, 'LineWidth', 2)

% ballast marker in 3D
ballastCol = [0.8, 0.5, 0.0];
if r.ballastMass_kg > 0
    plot3(r.ballastPos_m(1), r.ballastPos_m(2), r.ballastPos_m(3), 'd', ...
        'Color', ballastCol, 'MarkerSize', 10, 'MarkerFaceColor', ballastCol, 'LineWidth', 1.5)
end

xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m] (down +ve)')
title('3D View')
axis equal; grid on
set(gca, 'ZDir', 'reverse')  % z-down convention
view(30, 20)
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PANEL 2: Side Profile (top-right)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2, 2, 2);
hold on

% hull outline (upper and lower profile)
% z-down: top of hull at z = -radius, bottom at z = +radius
hullX = hd.hullPointArray_m(1, :);
fill([hullX, fliplr(hullX)], [-hullR, fliplr(hullR)], hullCol, ...
    'EdgeColor', [0.3, 0.6, 0.7], 'LineWidth', 1)

% fin side projections (x vs z)
for ii = 1:fd.numbFins
    fc = fd.finClouds{ii};
    plot(fc(1,:), fc(3,:), '.', 'MarkerSize', 1, 'Color', finCol)
end

% CV line
line([cvX, cvX], [-maxR*2, maxR*2], 'Color', cvCol, 'LineWidth', 2, 'LineStyle', '-')
% CG line
line([cgX, cgX], [-maxR*2, maxR*2], 'Color', cgCol, 'LineWidth', 2, 'LineStyle', '--')

% ballast marker
ballastCol = [0.8, 0.5, 0.0];
if r.ballastMass_kg > 0
    plot(r.ballastPos_m(1), r.ballastPos_m(3), 'd', 'Color', ballastCol, ...
        'MarkerSize', 10, 'MarkerFaceColor', ballastCol, 'LineWidth', 1.5)
end

% static margin annotation
midY = -maxR * 1.5;
plot([cvX, cgX], [midY, midY], 'k-', 'LineWidth', 1.5)
plot(cvX, midY, 'o', 'Color', cvCol, 'MarkerSize', 6, 'MarkerFaceColor', cvCol)
plot(cgX, midY, 's', 'Color', cgCol, 'MarkerSize', 6, 'MarkerFaceColor', cgCol)
text((cvX + cgX)/2, midY - maxR*0.3, ...
    sprintf('SM = %.3f m', r.staticMargin_m), ...
    'HorizontalAlignment', 'center', 'FontSize', 8)

xlabel('X [m] (nose to tail)')
ylabel('Z [m] (down +ve)')
title('Side Profile')
axis equal; grid on
set(gca, 'YDir', 'reverse')  % z down convention
legend('Hull', 'Fins', '', '', 'CV', 'CG', 'Location', 'northeast')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PANEL 3: Rear Cross-Section (bottom-left)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2, 2, 3);
hold on

% hull cross-section circle at the fin position
theta = linspace(0, 2*pi, 100);
% find hull radius at tail edge position
xTE = dp.xPosTailEdge_m;
hullXunique = unique(hullX);
hullRunique = hullR(1:length(hullXunique));
% simple interpolation for the radius at xTE
rAtTE = 0;
for jj = 1:length(hullXunique)-1
    if hullXunique(jj) <= xTE && hullXunique(jj+1) >= xTE
        t = (xTE - hullXunique(jj)) / (hullXunique(jj+1) - hullXunique(jj));
        rAtTE = hullRunique(jj) + t * (hullRunique(jj+1) - hullRunique(jj));
        break;
    end
end
plot(rAtTE * cos(theta), rAtTE * sin(theta), '-', 'Color', [0.3, 0.6, 0.7], 'LineWidth', 1.5)
fill(rAtTE * cos(theta), rAtTE * sin(theta), hullCol, 'FaceAlpha', 0.3, ...
    'EdgeColor', [0.3, 0.6, 0.7], 'LineWidth', 1.5)

% fin cross-sections (y vs z)
for ii = 1:fd.numbFins
    fc = fd.finClouds{ii};
    plot(fc(2,:), fc(3,:), '.', 'MarkerSize', 2, 'Color', finCol)
end

xlabel('Y [m]'); ylabel('Z [m] (down +ve)')
title(sprintf('Rear View at x = %.1f m', xTE))
axis equal; grid on
set(gca, 'YDir', 'reverse')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PANEL 4: Data Summary (bottom-right)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2, 2, 4);
axis off

% build text block
lines = {};
lines{end+1} = '\bf--- Geometry ---\rm';
lines{end+1} = sprintf('Hull: %.2f m x %.2f m (L x D)', dp.hullLength_m, dp.hullRadius_m*2);
lines{end+1} = sprintf('Volume: %.2f m^3', hd.hullVolume_m3);
lines{end+1} = sprintf('Fins: %d x (chord %.3f m, span %.2f m)', fd.numbFins, fd.baseCord_m, dp.tailingEdgeHight_m);
lines{end+1} = '';
lines{end+1} = '\bf--- Atmosphere ---\rm';
lines{end+1} = sprintf('Gas: %s', dp.liftingGas);
lines{end+1} = sprintf('Altitude: %.0f m', dp.altitude_m);
lines{end+1} = sprintf('Temperature: %.1f K (%.1f C)', dp.temperature_K, dp.temperature_K - 273.15);
lines{end+1} = sprintf('Air density: %.4f kg/m^3', dp.airDensity_kgpm3);
lines{end+1} = '';
lines{end+1} = '\bf--- Mass ---\rm';
lines{end+1} = sprintf('Hull: %.2f kg', r.hullTotalMass_kg);
lines{end+1} = sprintf('Fins: %.3f kg each (%.3f kg total)', r.singleFinMass_kg, r.singleFinMass_kg * fd.numbFins);
if r.ballastMass_kg > 0
    lines{end+1} = sprintf('\\color[rgb]{0.8,0.5,0}Ballast: %.3f kg (%s)\\color{black}', r.ballastMass_kg, dp.ballastOption);
end
lines{end+1} = sprintf('Total: %.2f kg', r.totalMass_kg);
lines{end+1} = '';
lines{end+1} = '\bf--- Buoyancy ---\rm';
lines{end+1} = sprintf('Net lift: %.1f N', r.netLift_N);
lines{end+1} = sprintf('Weight: %.1f N', r.totalWeight_N);
lines{end+1} = sprintf('Excess: %.1f N', r.excessLift_N);
if r.excessLift_N > 0
    lines{end+1} = '\color{green}\bfFLOATS\rm\color{black}';
else
    lines{end+1} = '\color{red}\bfSINKS\rm\color{black}';
end
lines{end+1} = '';
lines{end+1} = '\bf--- Stability ---\rm';
lines{end+1} = sprintf('CV: %.3f m  CG: %.3f m', cvX, cgX);
lines{end+1} = sprintf('Static margin: %.3f m', r.staticMargin_m);
lines{end+1} = sprintf('Static moment: %.2f Nm', r.staticMoment_Nm);
lines{end+1} = '';
lines{end+1} = '\bf--- Inertia (about CG) ---\rm';
lines{end+1} = sprintf('Ixx: %.1f  Iyy: %.1f  Izz: %.1f [kg.m^2]', ...
    r.inertiaRotCg_kgm2(1,1), r.inertiaRotCg_kgm2(2,2), r.inertiaRotCg_kgm2(3,3));

fullText = strjoin(lines, char(10));
text(0.05, 0.95, fullText, 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'FontSize', 9, ...
    'FontName', 'FixedWidth', 'Interpreter', 'tex')
title('Design Summary')

end