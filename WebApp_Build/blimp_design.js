/**
 * blimp_design.js
 *
 * Vanilla ES6 JavaScript translation of the MATLAB airship parameterised
 * design code. Runs entirely client-side with zero external dependencies.
 *
 * Coordinate convention: x forward (nose to tail), y starboard, z down.
 *
 * PUBLIC API:
 *   window.runBlimpDesign(inputs) -> results object
 *
 * -----------------------------------------------------------------------
 * ACCEPTED INPUT KEYS (with defaults, units, and expected ranges):
 *
 * @param {number} hullLength_m            - Hull length [m] (5..50)
 * @param {number} noseLength_m            - Nose ellipse semi-axis [m] (0.5..10)
 * @param {number} tailRatio               - Tail ellipse parameter [m] (1..10)
 * @param {number} hullRadius_m            - Hull max radius [m] (0.3..5)
 * @param {number} hullMaterialMass_kgpm2  - Hull skin areal density [kg/m^2] (0.01..5)
 * @param {number} altitude_m              - Flight altitude [m] (0..11000)
 * @param {number} temperature_K           - Ambient temperature [K] (200..330)
 * @param {string} liftingGas              - 'helium' or 'hydrogen'
 * @param {number} gravity_mps2            - Gravitational accel [m/s^2] (9.78..9.83)
 * @param {number} forwardAirSpeed_mps     - Forward airspeed [m/s] (0..100)
 * @param {number} LeadingEdgeTapperAngle_rad - Fin LE taper angle [rad] (0.05..0.8)
 * @param {number} tailingEdgeHight_m      - Fin span [m] (0.3..5)
 * @param {number} thickNessRatio          - NACA thickness ratio (0.05..0.30)
 * @param {number[]} nacaCoeffs            - 5-element NACA profile coefficients
 * @param {number} xPosTailEdge_m          - Fin trailing-edge x-position [m]
 * @param {number} dihedral_rad            - Fin dihedral angle [rad] (0..pi/2)
 * @param {number} numbFins                - Number of fins: 3 or 4
 * @param {number} finSaftyFactor_m        - Fin safety factor [m] (0..2)
 * @param {number} finMaterialDensity_kgpm3 - Fin material density [kg/m^3] (0.5..500)
 * @param {number} numRibs                 - Fin point-cloud rib count (10..500)
 * @param {string} ballastOption           - 'none', 'nose', or 'gondola'
 * @param {number} hullPointSpaceing_m     - Hull outline spacing [m]
 * @param {number} hullAngleSpaceing_rad   - Hull angular spacing [rad]
 * @param {number} finPointSpaceing_m      - Fin chordwise spacing [m]
 *
 * -----------------------------------------------------------------------
 * RETURNED OUTPUT KEYS (with units):
 *
 * @returns {number} hullVolume_m3         - Hull volume [m^3]
 * @returns {number} centerOfVolume_m      - Centre of volume x-position [m]
 * @returns {number} hullTotalMass_kg      - Hull skin mass [kg]
 * @returns {number} singleFinVolume_m3    - Single fin clipped volume [m^3]
 * @returns {number} singleFinMass_kg      - Single fin mass [kg]
 * @returns {number} totalMass_kg          - Total mass incl. ballast [kg]
 * @returns {number[]} centerOfGravity_m   - CG position [x,y,z] [m]
 * @returns {number[]} preBallastCG_m      - Pre-ballast CG [x,y,z] [m]
 * @returns {number} ballastMass_kg        - Required ballast [kg]
 * @returns {number[]} ballastPos_m        - Ballast position [x,y,z] [m]
 * @returns {number[][]} inertiaRotCg_kgm2 - 3x3 inertia tensor about CG [kg.m^2]
 * @returns {number[][]} virtualInertiaMatrixCv - 3x3 virtual inertia about CV
 * @returns {number} staticMargin_m        - Static margin CG-CV [m]
 * @returns {number} staticMoment_Nm       - Static pitching moment [Nm]
 * @returns {number} netLift_N             - Net buoyant lift [N]
 * @returns {number} buoyancyForce_N       - Gross buoyancy [N]
 * @returns {number} gasWeight_N           - Gas weight [N]
 * @returns {number} gasDensity_kgpm3      - Gas density at altitude [kg/m^3]
 * @returns {number} totalWeight_N         - Total weight [N]
 * @returns {number} excessLift_N          - Excess lift [N] (>0 = floats)
 * @returns {number} airDensity_kgpm3      - Air density at altitude [kg/m^3]
 * @returns {number} baseCord_m            - Fin base chord [m]
 * @returns {number} requiredArea_m2       - Required total fin area [m^2]
 * @returns {number} externalFraction_nd   - Fraction of fin outside hull
 * @returns {number} numbFins              - Number of fins
 * @returns {number} k1                    - Axial virtual mass coefficient
 * @returns {number} k2                    - Lateral virtual mass coefficient
 * @returns {number} kPrime                - Pitch/yaw rotational coefficient
 * @returns {object} hullPointArray        - {x:[], r:[]} hull outline
 * @returns {object[]} finClouds           - Array of {x:[],y:[],z:[]} per fin
 * @returns {object} hullPointCloud        - {x:[],y:[],z:[]} hull surface cloud
 */

// =========================================================================
// Helper: build an evenly-spaced array from start to end (inclusive, like MATLAB colon)
// =========================================================================
function linspace(start, stop, n) {
  const arr = new Array(n);
  if (n === 1) { arr[0] = start; return arr; }
  const step = (stop - start) / (n - 1);
  for (let i = 0; i < n; i++) arr[i] = start + i * step;
  return arr;
}

function colonArray(start, step, stop) {
  // Mimics MATLAB start:step:stop
  const arr = [];
  if (step > 0) {
    for (let v = start; v <= stop + step * 1e-10; v += step) arr.push(v);
  } else if (step < 0) {
    for (let v = start; v >= stop + step * 1e-10; v += step) arr.push(v);
  }
  return arr;
}

// =========================================================================
// Helper: 3x3 matrix operations (row-major as [row][col])
// =========================================================================
function mat3Mul(A, B) {
  // 3x3 * 3x3
  const C = [[0,0,0],[0,0,0],[0,0,0]];
  for (let i = 0; i < 3; i++)
    for (let j = 0; j < 3; j++)
      for (let k = 0; k < 3; k++)
        C[i][j] += A[i][k] * B[k][j];
  return C;
}

function mat3MulVec(M, v) {
  // 3x3 * 3x1 -> 3x1
  return [
    M[0][0]*v[0] + M[0][1]*v[1] + M[0][2]*v[2],
    M[1][0]*v[0] + M[1][1]*v[1] + M[1][2]*v[2],
    M[2][0]*v[0] + M[2][1]*v[1] + M[2][2]*v[2]
  ];
}

function scalarMulMat3(s, M) {
  return M.map(row => row.map(v => v * s));
}

// =========================================================================
// Source: tools/computeDCM.m / computeDCM
// =========================================================================
function computeDCM(roll_rad, pitch_rad, yaw_rad) {
  const cr = Math.cos(roll_rad),  sr = Math.sin(roll_rad);
  const cp = Math.cos(pitch_rad), sp = Math.sin(pitch_rad);
  const cy = Math.cos(yaw_rad),   sy = Math.sin(yaw_rad);

  const matrix_roll  = [[1, 0,  0 ],
                         [0, cr, sr],
                         [0,-sr, cr]];
  const matrix_pitch = [[cp, 0, -sp],
                         [0,  1,  0 ],
                         [sp, 0,  cp]];
  const matrix_yaw   = [[ cy, sy, 0],
                         [-sy, cy, 0],
                         [  0,  0, 1]];
  // DIRECTION Z -> Y -> X
  return mat3Mul(matrix_yaw, mat3Mul(matrix_pitch, matrix_roll));
}

// =========================================================================
// Source: tools/atmosphereModel.m / atmosphereModel
// =========================================================================
function atmosphereModel(altitude_m, temperature_K) {
  const P0 = 101325;    // sea level standard pressure [Pa]
  const T0 = 288.15;    // sea level standard temperature [K]
  const L  = 0.0065;    // temperature lapse rate [K/m]
  const g  = 9.81;      // gravitational acceleration [m/s^2]
  const R  = 287.058;   // specific gas constant for dry air [J/(kg.K)]

  const pressure_Pa = P0 * Math.pow(1 - L * altitude_m / T0, g / (R * L));
  const airDensity_kgpm3 = pressure_Pa / (R * temperature_K);

  return { airDensity_kgpm3, pressure_Pa };
}

// =========================================================================
// Source: tools/computeVirtualInertiaCoeffs.m / computeVirtualInertiaCoeffs
// =========================================================================
function computeVirtualInertiaCoeffs(hullLength_m, hullRadius_m) {
  const a = hullLength_m / 2;  // semi-major axis
  const b = hullRadius_m;       // semi-minor axis

  const e = Math.sqrt(1 - (b / a) * (b / a));
  const e2 = e * e;
  const e3 = e2 * e;

  // Lamb's integrals
  const alpha_0 = (2 * (1 - e2) / e3) * (0.5 * Math.log((1 + e) / (1 - e)) - e);
  const beta_0  = 1 / e2 - ((1 - e2) / (2 * e3)) * Math.log((1 + e) / (1 - e));

  // Translational virtual mass coefficients
  const k_x = alpha_0 / (2 - alpha_0);
  const k_y = beta_0  / (2 - beta_0);
  const k_z = k_y;

  const kMatrix_nd = [[k_x, 0,   0  ],
                       [0,   k_y, 0  ],
                       [0,   0,   k_z]];

  // Rotational virtual inertia coefficients
  const k_prime_x = 0;
  const k_prime_y = (e2 * (beta_0 - alpha_0)) /
    ((2 - e2) * (2 * e2 - (2 - e2) * (beta_0 - alpha_0)));
  const k_prime_z = k_prime_y;

  const kPrimeMatrix_nd = [[k_prime_x, 0,         0        ],
                            [0,         k_prime_y, 0        ],
                            [0,         0,         k_prime_z]];

  return { kMatrix_nd, kPrimeMatrix_nd };
}

// =========================================================================
// Source: tools/computeBuoyancy.m / computeBuoyancy
// =========================================================================
function computeBuoyancy(hullVolume_m3, airDensity_kgpm3, gasDensitySL_kgpm3,
                          temperature_K, gravity_mps2) {
  const T0 = 288.15;
  const P0 = 101325;
  const R_air = 287.058;

  const pressure_Pa = airDensity_kgpm3 * R_air * temperature_K;
  const gasDensity_kgpm3 = gasDensitySL_kgpm3 * (pressure_Pa / P0) * (T0 / temperature_K);

  const buoyancyForce_N = airDensity_kgpm3 * hullVolume_m3 * gravity_mps2;
  const gasWeight_N     = gasDensity_kgpm3 * hullVolume_m3 * gravity_mps2;
  const netLift_N       = buoyancyForce_N - gasWeight_N;

  return { netLift_N, buoyancyForce_N, gasWeight_N, gasDensity_kgpm3 };
}

// =========================================================================
// Source: tools/finVolume.m / finVolume
// =========================================================================
function finVolume(thickNessRatio, baseCord_m, LeadingEdgeTapperAngle_rad,
                   tailingEdgeHight_m, nacaCoeffs) {
  const k1 = nacaCoeffs[0]; // 0-based indexing (MATLAB: nacaCoeffs(1))
  const k2 = nacaCoeffs[1];
  const k3 = nacaCoeffs[2];
  const k4 = nacaCoeffs[3];
  const k5 = nacaCoeffs[4];

  // A(c) = areaCoeff * c^2
  // A = 2 * (t*c^2/0.2) * (k1*2/3 + k2/2 + k3/3 + k4/4 + k5/5)
  const areaCoeff = (thickNessRatio / 0.2) * 2 *
    (k1 * 2 / 3 + k2 / 2 + k3 / 3 + k4 / 4 + k5 / 5);

  const fullHight_m = Math.tan(Math.PI / 2 - LeadingEdgeTapperAngle_rad) * baseCord_m;

  const tanAlpha = Math.tan(LeadingEdgeTapperAngle_rad);
  const finVolume_m3 = areaCoeff * tanAlpha * tanAlpha / 3 *
    (Math.pow(fullHight_m, 3) - Math.pow(fullHight_m - tailingEdgeHight_m, 3));

  return finVolume_m3;
}

// =========================================================================
// Source: tools/volumetricMoment.m / volumetricMoment
// =========================================================================
function volumetricMoment(cv_m, hullX, hullR) {
  const n = hullX.length;
  let sum = 0;
  for (let i = 0; i < n; i++) {
    const dx = (i < n - 1) ? (hullX[i + 1] - hullX[i]) : (hullX[n - 1] - hullX[n - 2]);
    sum += hullR[i] * hullR[i] * (hullX[i] - cv_m) * (hullX[i] - cv_m) * dx;
  }
  return Math.PI * sum;
}

// =========================================================================
// Source: tools/computeInertiaMatrix.m / computeInertiaMatrix
// Point masses: locX[i], locY[i], locZ[i], mass[i]
// =========================================================================
function computeInertiaMatrix(locX, locY, locZ, mass) {
  const N = mass.length;
  const totalMass = mass.reduce((a, b) => a + b, 0);

  // CG
  let cgx = 0, cgy = 0, cgz = 0;
  for (let i = 0; i < N; i++) {
    cgx += locX[i] * mass[i];
    cgy += locY[i] * mass[i];
    cgz += locZ[i] * mass[i];
  }
  cgx /= totalMass;
  cgy /= totalMass;
  cgz /= totalMass;

  // Inertia about CG
  let Ixx = 0, Iyy = 0, Izz = 0, Ixy = 0, Ixz = 0, Iyz = 0;
  for (let i = 0; i < N; i++) {
    const x = locX[i] - cgx;
    const y = locY[i] - cgy;
    const z = locZ[i] - cgz;
    const m = mass[i];
    Ixx += m * (y * y + z * z);
    Iyy += m * (x * x + z * z);
    Izz += m * (x * x + y * y);
    Ixy += m * x * y;
    Ixz += m * x * z;
    Iyz += m * y * z;
  }

  // NOTE: off-diagonal products use the positive convention (Ixy = sum(m*x*y))
  const inertia = [[Ixx, Ixy, Ixz],
                    [Ixy, Iyy, Iyz],
                    [Ixz, Iyz, Izz]];

  return { inertia, cg: [cgx, cgy, cgz] };
}

// =========================================================================
// Source: hull/hullOutline.m / hullOutline
// =========================================================================
function hullOutline(dp) {
  const spaceing_m   = dp.hullPointSpaceing_m;
  const noseLength_m = dp.noseLength_m;
  const length_m     = dp.hullLength_m;
  const tailRatio    = dp.tailRatio;
  const radius_m     = dp.hullRadius_m;

  const tailSemiAxis_m = noseLength_m + tailRatio;
  const centreEnd_m    = length_m - tailSemiAxis_m;

  const xNose   = colonArray(0,              spaceing_m, noseLength_m);
  const xCentre = colonArray(noseLength_m,   spaceing_m, centreEnd_m);
  const xTail   = colonArray(centreEnd_m,    spaceing_m, length_m);

  const rNose = xNose.map(x => {
    const d = noseLength_m - x;
    return radius_m * Math.sqrt(Math.max(0, 1 - (d * d) / (noseLength_m * noseLength_m)));
  });
  const rCentre = xCentre.map(() => radius_m);
  const rTail = xTail.map(x => {
    const d = centreEnd_m - x;
    return radius_m * Math.sqrt(Math.max(0, 1 - (d * d) / (tailSemiAxis_m * tailSemiAxis_m)));
  });

  return {
    x: xNose.concat(xCentre, xTail),
    r: rNose.concat(rCentre, rTail)
  };
}

// =========================================================================
// Source: hull/hullVolume.m / hullVolume
// =========================================================================
function hullVolume(hullX, hullR) {
  let vol = 0;
  for (let i = 1; i < hullX.length; i++) { // 0-based: starts at 1 (MATLAB: 2:end)
    const dx = hullX[i] - hullX[i - 1];
    vol += hullR[i] * hullR[i] * dx;
  }
  return Math.PI * vol;
}

// =========================================================================
// Source: hull/centerOfVolume.m / centerOfVolume
// =========================================================================
function centerOfVolume(hullX, hullR) {
  let sumAx = 0, sumA = 0;
  for (let i = 0; i < hullX.length; i++) {
    const a = Math.PI * hullR[i] * hullR[i];
    sumAx += a * hullX[i];
    sumA  += a;
  }
  return sumAx / sumA;
}

// =========================================================================
// Source: hull/genHullPointCloud.m / genHullPointCloud
// =========================================================================
function genHullPointCloud(hullX, hullR, anglesDelta_rad) {
  const n = hullX.length;
  const rollAngles = colonArray(0, anglesDelta_rad, 2 * Math.PI);
  const totalPts = n * rollAngles.length;
  const cx = new Float64Array(totalPts);
  const cy = new Float64Array(totalPts);
  const cz = new Float64Array(totalPts);

  for (let ii = 0; ii < rollAngles.length; ii++) {
    const DCM = computeDCM(rollAngles[ii], 0, 0);
    const base = ii * n;
    for (let j = 0; j < n; j++) {
      // Input vector: [hullX[j], hullR[j], 0]
      const v = mat3MulVec(DCM, [hullX[j], hullR[j], 0]);
      cx[base + j] = v[0];
      cy[base + j] = v[1];
      cz[base + j] = v[2];
    }
  }
  return { x: cx, y: cy, z: cz, length: totalPts };
}

// =========================================================================
// Source: hull/hullTopFunction.m / hullTopFunction
// =========================================================================
function hullTopFunction(dp) {
  const outline = hullOutline(dp);
  const vol     = hullVolume(outline.x, outline.r);
  const cv      = centerOfVolume(outline.x, outline.r);
  const cloud   = genHullPointCloud(outline.x, outline.r, dp.hullAngleSpaceing_rad);

  return {
    hullPointArray: outline,    // {x:[], r:[]}
    hullVolume_m3: vol,
    centerOfVolume_m: cv,
    hullPointCloud: cloud       // {x:F64, y:F64, z:F64, length}
  };
}

// =========================================================================
// Source: tools/hullPointCloudDensity.m / hullPointCloudDensity
// =========================================================================
function hullPointCloudDensity(hullX, hullR, anglesDelta_rad) {
  const n = hullX.length;

  // dx and dr along the profile
  const dx = new Array(n);
  const dr = new Array(n);
  for (let i = 0; i < n - 1; i++) {
    dx[i] = hullX[i + 1] - hullX[i];
    dr[i] = hullR[i + 1] - hullR[i];
  }
  dx[n - 1] = dx[n - 2]; // replicate last (MATLAB: [dx_m, dx_m(end)])
  dr[n - 1] = dr[n - 2];

  // Arc length along profile curve per segment
  const h = new Array(n);
  for (let i = 0; i < n; i++) {
    const theta = Math.atan2(dr[i], dx[i]);
    const cosT = Math.cos(theta);
    h[i] = (Math.abs(cosT) > 1e-15) ? dx[i] / cosT : 0;
  }

  // Circumferential arc length per angular slice
  const s = hullR.map(r => r * anglesDelta_rad);

  // Area per point
  const areaArray = new Array(n);
  for (let i = 0; i < n; i++) areaArray[i] = h[i] * s[i];

  // Replicate for each angular slice
  const rollAngles = colonArray(0, anglesDelta_rad, 2 * Math.PI);
  const totalPts = n * rollAngles.length;
  const areaPc = new Float64Array(totalPts);
  for (let ii = 0; ii < rollAngles.length; ii++) {
    const base = ii * n;
    for (let j = 0; j < n; j++) {
      areaPc[base + j] = areaArray[j];
    }
  }
  return areaPc;
}

// =========================================================================
// Source: fins/finSizeing.m / finSizeing
// =========================================================================
function finSizeing(dp, hd) {
  const momentArm_m = dp.xPosTailEdge_m - hd.centerOfVolume_m - dp.finSaftyFactor_m;

  // required area
  const requiredArea_m2 = hd.hullVolume_m3 * (dp.k2 - dp.k1) / (1.5 * momentArm_m) * 0.384 * 2;

  const inner = dp.hullLength_m - dp.noseLength_m - dp.tailRatio - dp.xPosTailEdge_m;
  const baseHight_m = dp.hullRadius_m *
    Math.sqrt(Math.max(0, 1 - (inner * inner) / ((dp.noseLength_m + dp.tailRatio) * (dp.noseLength_m + dp.tailRatio))));

  let baseCord_m;
  const teh = dp.tailingEdgeHight_m;
  const tanLE = Math.tan(dp.LeadingEdgeTapperAngle_rad);
  const sumH = teh + baseHight_m;

  if (dp.numbFins === 3) {
    baseCord_m = requiredArea_m2 / (sumH * (1 + Math.cos(Math.PI / 2 - dp.dihedral_rad)))
      + teh * teh * tanLE / (2 * sumH);
  } else if (dp.numbFins === 4) {
    baseCord_m = requiredArea_m2 / (sumH * (2 * Math.cos(Math.PI / 2 - dp.dihedral_rad)))
      + teh * teh * tanLE / (2 * sumH);
  } else {
    throw new Error('numbFins must be 3 or 4, got ' + dp.numbFins);
  }

  return { requiredArea_m2, baseHight_m, baseCord_m };
}

// =========================================================================
// Source: fins/getFinCrossSection.m / getFinCrossSection
// =========================================================================
function getFinCrossSection(thickNessRatio, cordLength_m, spaceing_m, nacaCoeffs) {
  const k1 = nacaCoeffs[0]; // 0-based (MATLAB: nacaCoeffs(1))
  const k2 = nacaCoeffs[1];
  const k3 = nacaCoeffs[2];
  const k4 = nacaCoeffs[3];
  const k5 = nacaCoeffs[4];

  // x from -cordLength to 0
  const xArr = colonArray(-cordLength_m, spaceing_m, 0);
  const nPts = xArr.length;

  // upper surface
  const xOut = [];
  const yOut = [];
  const yUpper = new Array(nPts);
  for (let i = 0; i < nPts; i++) {
    const var1 = (xArr[i] + cordLength_m) / cordLength_m;
    yUpper[i] = (thickNessRatio * cordLength_m) / 0.2 *
      (k1 * Math.pow(Math.max(0, var1), 0.5) + k2 * var1 +
       k3 * var1 * var1 + k4 * var1 * var1 * var1 +
       k5 * var1 * var1 * var1 * var1);
    xOut.push(xArr[i]);
    yOut.push(yUpper[i]);
  }

  // lower surface (reversed, like MATLAB flip([xArray;-yArray],2))
  for (let i = nPts - 1; i >= 0; i--) {
    xOut.push(xArr[i]);
    yOut.push(-yUpper[i]);
  }

  return { x: xOut, y: yOut };
}

// =========================================================================
// Source: fins/buildFinPointCloud.m / buildFinPointCloud
// =========================================================================
function buildFinPointCloud(dp, finSizeResult) {
  const LeadingEdgeTapperAngle_rad = dp.LeadingEdgeTapperAngle_rad;
  const tailingEdgeHight_m = dp.tailingEdgeHight_m;
  const thickNessRatio = dp.thickNessRatio;
  const nacaCoeffs = dp.nacaCoeffs;
  const spaceing_m = dp.finPointSpaceing_m;
  const baseCord_m = finSizeResult.baseCord_m;
  const numRibs = dp.numRibs;

  const ribSpaceing_m = tailingEdgeHight_m / (numRibs - 1);
  const hightArray_m = colonArray(0, ribSpaceing_m, tailingEdgeHight_m);
  const fullHight_m = Math.tan(Math.PI / 2 - LeadingEdgeTapperAngle_rad) * baseCord_m;

  const cloudX = [];
  const cloudY = [];
  const cloudZ = [];

  for (let ii = 0; ii < hightArray_m.length; ii++) { // 0-based (MATLAB: 1:length)
    const distanceDownFin_m = hightArray_m[ii];
    const distanceRemaining_m = fullHight_m - distanceDownFin_m;
    const baseCordLength_m = Math.tan(LeadingEdgeTapperAngle_rad) * distanceRemaining_m;

    if (baseCordLength_m <= spaceing_m) continue; // skip degenerate ribs

    const cs = getFinCrossSection(thickNessRatio, baseCordLength_m, spaceing_m, nacaCoeffs);
    for (let j = 0; j < cs.x.length; j++) {
      cloudX.push(cs.x[j]);
      cloudY.push(cs.y[j]);
      cloudZ.push(distanceDownFin_m);
    }
  }

  return { x: cloudX, y: cloudY, z: cloudZ };
}

// =========================================================================
// Source: fins/positionFins.m / positionFins
// =========================================================================
function positionFins(finCloud, xPosTailEdge_m, zPosBase_m, rotation_rad) {
  const n = finCloud.x.length;
  const outX = new Array(n);
  const outY = new Array(n);
  const outZ = new Array(n);
  const DCM = computeDCM(rotation_rad, 0, 0);

  for (let i = 0; i < n; i++) {
    // 1. offset z so fin base is at hull radius
    const px = finCloud.x[i];
    const py = finCloud.y[i];
    const pz = finCloud.z[i] + zPosBase_m;

    // 2. rotate around x-axis
    const rot = mat3MulVec(DCM, [px, py, pz]);

    // 3. translate x to tail position
    outX[i] = rot[0] + xPosTailEdge_m;
    outY[i] = rot[1];
    outZ[i] = rot[2];
  }

  return { x: outX, y: outY, z: outZ };
}

// =========================================================================
// Source: tools/clipFinToHull.m / clipFinToHull
// =========================================================================
function clipFinToHull(finCloud, hullX, hullR) {
  // Remove duplicate x-values (section boundaries)
  const ux = [], ur = [];
  for (let i = 0; i < hullX.length; i++) {
    if (i === 0 || hullX[i] !== ux[ux.length - 1]) {
      ux.push(hullX[i]);
      ur.push(hullR[i]);
    }
  }

  const n = finCloud.x.length;
  const numHull = ux.length;
  const outX = [], outY = [], outZ = [];

  for (let i = 0; i < n; i++) {
    const xq = finCloud.x[i];
    const rFin = Math.sqrt(finCloud.y[i] * finCloud.y[i] + finCloud.z[i] * finCloud.z[i]);

    let rHullAt = 0;
    if (xq > ux[0] && xq < ux[numHull - 1]) {
      // find interval (0-based, MATLAB jj=1 -> j=0)
      let j = 0;
      while (j < numHull - 1 && ux[j + 1] < xq) j++;
      const t = (xq - ux[j]) / (ux[j + 1] - ux[j]);
      rHullAt = ur[j] + t * (ur[j + 1] - ur[j]);
    }

    if (rFin >= rHullAt) {
      outX.push(finCloud.x[i]);
      outY.push(finCloud.y[i]);
      outZ.push(finCloud.z[i]);
    }
  }

  return { x: outX, y: outY, z: outZ };
}

// =========================================================================
// Source: fins/finTopFunction.m / finTopFunction
// =========================================================================
function finTopFunction(dp, hd) {
  const finSize = finSizeing(dp, hd);
  const finCloud = buildFinPointCloud(dp, finSize);

  const zPosBase_m = finSize.baseHight_m;
  const xPosTailEdge_m = dp.xPosTailEdge_m;
  const dihedral_rad = dp.dihedral_rad;
  const numPointsUnclipped = finCloud.x.length;

  // Coordinate convention: x forward, y starboard, z down
  // Fin is built with span in +z. Roll around x-axis places fins:
  //   roll = 0    -> fin points down (+z)
  //   roll = pi   -> fin points up   (-z)
  let finAngles_rad;
  if (dp.numbFins === 3) {
    // inverted-Y: one fin up, two lower fins splayed by dihedral
    finAngles_rad = [
      Math.PI,                  // top fin (up, -z)
      -dihedral_rad,            // bottom right
       dihedral_rad             // bottom left
    ];
  } else if (dp.numbFins === 4) {
    // X-tail: two upper, two lower
    finAngles_rad = [
      Math.PI + dihedral_rad,   // upper right
      Math.PI - dihedral_rad,   // upper left
      -dihedral_rad,            // lower right
       dihedral_rad             // lower left
    ];
  } else {
    throw new Error('numbFins must be 3 or 4, got ' + dp.numbFins);
  }

  const finClouds = [];
  for (let ii = 0; ii < finAngles_rad.length; ii++) {
    const raw = positionFins(finCloud, xPosTailEdge_m, zPosBase_m, finAngles_rad[ii]);
    const clipped = clipFinToHull(raw, hd.hullPointArray.x, hd.hullPointArray.r);
    finClouds.push(clipped);
  }

  const numPointsClipped = finClouds[0].x.length;
  const externalFraction_nd = numPointsClipped / numPointsUnclipped;

  return {
    baseCord_m: finSize.baseCord_m,
    requiredArea_m2: finSize.requiredArea_m2,
    baseHight_m: finSize.baseHight_m,
    finClouds,
    numbFins: finAngles_rad.length,
    externalFraction_nd
  };
}

// =========================================================================
// Source: tools/vortexPanelMethod.m  (Vortex Panel Method)
// Adapted from Vortex-Panel-Method-for-Airfoil-Aerodynamic-Coefficients
// =========================================================================
/**
 * Vortex panel method for 2-D airfoil aerodynamic coefficients.
 * @param {number[]} Xb - boundary point x-coordinates (Np)
 * @param {number[]} Yb - boundary point y-coordinates (Np)
 * @param {number} Vinf - freestream velocity [m/s]
 * @param {number} c - chord length [m]
 * @param {number} alpha_deg - angle of attack [deg]
 * @returns {object} { Cl, Cd, Cp[], V[], xc[], yc[], panelS[] }
 */
function vortexPanelMethod(Xb, Yb, Vinf, c, alpha_deg) {
  var Mp = Xb.length;
  var M  = Mp - 1;
  var alpha = alpha_deg * Math.PI / 180;

  // Panel geometry
  var x = new Array(M), y = new Array(M), s = new Array(M);
  var theta = new Array(M), sine = new Array(M), cosine = new Array(M);
  var RHS = new Array(Mp);
  for (var i = 0; i < M; i++) {
    var ip = i + 1;
    x[i] = 0.5 * (Xb[i] + Xb[ip]);
    y[i] = 0.5 * (Yb[i] + Yb[ip]);
    s[i] = Math.sqrt((Xb[ip] - Xb[i]) * (Xb[ip] - Xb[i]) + (Yb[ip] - Yb[i]) * (Yb[ip] - Yb[i]));
    theta[i]  = Math.atan2(Yb[ip] - Yb[i], Xb[ip] - Xb[i]);
    sine[i]   = Math.sin(theta[i]);
    cosine[i] = Math.cos(theta[i]);
    RHS[i]    = Math.sin(theta[i] - alpha);
  }

  // Influence coefficient sub-matrices
  var CN1 = [], CN2 = [], CT1 = [], CT2 = [];
  for (var i = 0; i < M; i++) {
    CN1[i] = new Array(M); CN2[i] = new Array(M);
    CT1[i] = new Array(M); CT2[i] = new Array(M);
    for (var j = 0; j < M; j++) {
      if (i === j) {
        CN1[i][j] = -1.0; CN2[i][j] = 1.0;
        CT1[i][j] = 0.5 * Math.PI; CT2[i][j] = 0.5 * Math.PI;
      } else {
        var A = -(x[i] - Xb[j]) * cosine[j] - (y[i] - Yb[j]) * sine[j];
        var B = (x[i] - Xb[j]) * (x[i] - Xb[j]) + (y[i] - Yb[j]) * (y[i] - Yb[j]);
        var C = Math.sin(theta[i] - theta[j]);
        var D = Math.cos(theta[i] - theta[j]);
        var E = (x[i] - Xb[j]) * sine[j] - (y[i] - Yb[j]) * cosine[j];
        var F = Math.log(1.0 + s[j] * (s[j] + 2 * A) / B);
        var G = Math.atan2(E * s[j], B + A * s[j]);
        var P = (x[i] - Xb[j]) * Math.sin(theta[i] - 2 * theta[j])
              + (y[i] - Yb[j]) * Math.cos(theta[i] - 2 * theta[j]);
        var Q = (x[i] - Xb[j]) * Math.cos(theta[i] - 2 * theta[j])
              - (y[i] - Yb[j]) * Math.sin(theta[i] - 2 * theta[j]);
        CN2[i][j] = D + 0.5 * Q * F / s[j] - (A * C + D * E) * G / s[j];
        CN1[i][j] = 0.5 * D * F + C * G - CN2[i][j];
        CT2[i][j] = C + 0.5 * P * F / s[j] + (A * D - C * E) * G / s[j];
        CT1[i][j] = 0.5 * C * F - D * G - CT2[i][j];
      }
    }
  }

  // Assemble AN (Mp x Mp) and AT (M x Mp)
  var AN = [], AT = [];
  for (var i = 0; i < Mp; i++) { AN[i] = new Array(Mp).fill(0); }
  for (var i = 0; i < M; i++) { AT[i] = new Array(Mp).fill(0); }

  for (var i = 0; i < M; i++) {
    AN[i][0]    = CN1[i][0];
    AN[i][Mp-1] = CN2[i][M-1];
    AT[i][0]    = CT1[i][0];
    AT[i][Mp-1] = CT2[i][M-1];
    for (var j = 1; j < M; j++) { // j from 1 to M-1 (0-based), MATLAB j from 2 to M
      AN[i][j] = CN1[i][j] + CN2[i][j - 1];
      AT[i][j] = CT1[i][j] + CT2[i][j - 1];
    }
  }
  // Kutta condition
  AN[Mp-1][0]    = 1.0;
  AN[Mp-1][Mp-1] = 1.0;
  for (var j = 1; j < M; j++) AN[Mp-1][j] = 0.0;
  RHS[Mp-1] = 0.0;

  // Solve AN * Gama = RHS using Gaussian elimination with partial pivoting
  var n = Mp;
  var aug = [];
  for (var i = 0; i < n; i++) {
    aug[i] = new Array(n + 1);
    for (var j = 0; j < n; j++) aug[i][j] = AN[i][j];
    aug[i][n] = RHS[i];
  }
  for (var col = 0; col < n; col++) {
    // Partial pivot
    var maxRow = col, maxVal = Math.abs(aug[col][col]);
    for (var row = col + 1; row < n; row++) {
      if (Math.abs(aug[row][col]) > maxVal) { maxVal = Math.abs(aug[row][col]); maxRow = row; }
    }
    if (maxRow !== col) { var tmp = aug[col]; aug[col] = aug[maxRow]; aug[maxRow] = tmp; }
    var piv = aug[col][col];
    if (Math.abs(piv) < 1e-30) continue;
    for (var row = col + 1; row < n; row++) {
      var factor = aug[row][col] / piv;
      for (var k = col; k <= n; k++) aug[row][k] -= factor * aug[col][k];
    }
  }
  var Gama = new Array(n);
  for (var i = n - 1; i >= 0; i--) {
    var sum = aug[i][n];
    for (var j = i + 1; j < n; j++) sum -= aug[i][j] * Gama[j];
    Gama[i] = (Math.abs(aug[i][i]) > 1e-30) ? sum / aug[i][i] : 0;
  }

  // Tangential velocity and pressure coefficient
  var V  = new Array(M);
  var Cp = new Array(M);
  for (var i = 0; i < M; i++) {
    V[i] = Math.cos(theta[i] - alpha);
    for (var j = 0; j < Mp; j++) V[i] += AT[i][j] * Gama[j];
    Cp[i] = 1.0 - V[i] * V[i];
  }

  // Force coefficients from pressure integration
  var cy = 0, cx = 0;
  for (var i = 0; i < M; i++) {
    cy -= Cp[i] * s[i] * cosine[i];
    cx -= Cp[i] * s[i] * (-sine[i]);
  }
  cy /= c;
  cx /= c;
  var Cd =  cx * Math.cos(alpha) + cy * Math.sin(alpha);
  var Cl = -cx * Math.sin(alpha) + cy * Math.cos(alpha);

  return { Cl: Cl, Cd: Cd, Cp: Cp, V: V, xc: x, yc: y, panelS: s, Gama: Gama };
}

// =========================================================================
// Source: tools/runPanelAnalysis.m  (run VPM on hull profile & fin airfoil)
// =========================================================================
function runPanelAnalysis(dp, hd, fd, alpha_deg) {
  var Vinf = dp.forwardAirSpeed_mps;

  // --- Hull profile as 2D closed body (x-r plane) ---
  var hx = hd.hullPointArray.x;
  var hr = hd.hullPointArray.r;
  var nH = hx.length;

  // Assemble: TE -> lower(tail to nose, -r) -> upper(nose to tail, +r) -> TE
  var Xb_raw = [], Yb_raw = [];
  for (var i = nH - 1; i >= 0; i--) { Xb_raw.push(hx[i]); Yb_raw.push(-hr[i]); }
  for (var i = 1; i < nH; i++)      { Xb_raw.push(hx[i]); Yb_raw.push(hr[i]); }

  // Remove duplicate consecutive points (zero-length panels cause NaN)
  var Xb_hull = [Xb_raw[0]], Yb_hull = [Yb_raw[0]];
  for (var i = 1; i < Xb_raw.length; i++) {
    var dx = Xb_raw[i] - Xb_raw[i-1], dy = Yb_raw[i] - Yb_raw[i-1];
    if (dx * dx + dy * dy > 1e-20) { Xb_hull.push(Xb_raw[i]); Yb_hull.push(Yb_raw[i]); }
  }

  var c_hull = dp.hullLength_m;
  var hullRes = vortexPanelMethod(Xb_hull, Yb_hull, Vinf, c_hull, alpha_deg);

  // --- Fin airfoil cross-section ---
  var finCS = getFinCrossSection(dp.thickNessRatio, fd.baseCord_m, dp.finPointSpaceing_m, dp.nacaCoeffs);

  // Shift x so range is 0 to chord (finCS.x goes from -chord to 0, then 0 to -chord)
  var Xb_fin = [], Yb_fin = [];
  for (var i = 0; i < finCS.x.length; i++) {
    Xb_fin.push(finCS.x[i] + fd.baseCord_m);
    Yb_fin.push(finCS.y[i]);
  }

  var c_fin = fd.baseCord_m;
  var finRes = vortexPanelMethod(Xb_fin, Yb_fin, Vinf, c_fin, alpha_deg);

  return {
    hull: {
      Cl: hullRes.Cl, Cd: hullRes.Cd, Cp: hullRes.Cp, V: hullRes.V,
      xc: hullRes.xc, yc: hullRes.yc, panelS: hullRes.panelS,
      Gama: hullRes.Gama, Xb: Xb_hull, Yb: Yb_hull, chord: c_hull, Vinf: Vinf
    },
    fin: {
      Cl: finRes.Cl, Cd: finRes.Cd, Cp: finRes.Cp, V: finRes.V,
      xc: finRes.xc, yc: finRes.yc, panelS: finRes.panelS,
      Gama: finRes.Gama, Xb: Xb_fin, Yb: Yb_fin, chord: c_fin, Vinf: Vinf
    },
    alpha_deg: alpha_deg
  };
}

// =========================================================================
// Source: runDesign.m / runBlimpDesign (entry point)
// =========================================================================
window.runBlimpDesign = function(userInputs) {
  // -------------------------------------------------------------------
  // 1. Design Parameters (merge user inputs with defaults)
  // Source: designPamameters/designParametersCaseOne.m
  // -------------------------------------------------------------------
  const dp = {
    // hull
    hullLength_m:            12.75,
    noseLength_m:            2,
    tailRatio:               3.5,
    hullRadius_m:            1.22,
    hullMaterialMass_kgpm2:  0.12,
    // atmosphere
    altitude_m:              0,
    temperature_K:           288.15,
    liftingGas:              'helium',
    gravity_mps2:            9.81,
    forwardAirSpeed_mps:     10,
    // fins
    LeadingEdgeTapperAngle_rad: 15 * Math.PI / 180,  // deg2rad(15)
    tailingEdgeHight_m:      1,
    thickNessRatio:          0.15,
    nacaCoeffs:              [0.29690, -0.12600, -0.35160, 0.28430, -0.10150],
    xPosTailEdge_m:          11.75,
    dihedral_rad:            50 * Math.PI / 180,       // deg2rad(50)
    numbFins:                3,
    finSaftyFactor_m:        0.25,
    finMaterialDensity_kgpm3: 3,
    numRibs:                 100,
    // ballast
    ballastOption:           'gondola',
    // spacing
    hullPointSpaceing_m:     0.1 / 6,
    hullAngleSpaceing_rad:   0.01,
    finPointSpaceing_m:      0.012
  };

  // Override defaults with user-supplied values
  if (userInputs) {
    for (const key in userInputs) {
      if (userInputs.hasOwnProperty(key)) {
        dp[key] = userInputs[key];
      }
    }
  }

  // Gas density at sea level
  if (dp.liftingGas === 'helium') {
    dp.liftingGasDensitySL_kgpm3 = 0.1786;
  } else if (dp.liftingGas === 'hydrogen') {
    dp.liftingGasDensitySL_kgpm3 = 0.0899;
  } else {
    throw new Error("liftingGas must be 'helium' or 'hydrogen', got '" + dp.liftingGas + "'");
  }

  // Virtual inertia coefficients (Lamb's k-factors)
  const { kMatrix_nd, kPrimeMatrix_nd } = computeVirtualInertiaCoeffs(dp.hullLength_m, dp.hullRadius_m);
  dp.k1 = kMatrix_nd[0][0]; // 0-based (MATLAB: kMatrix_nd(1,1))
  dp.k2 = kMatrix_nd[1][1]; // 0-based (MATLAB: kMatrix_nd(2,2))
  dp.vertialInertiaCoeff_nd = kPrimeMatrix_nd;
  dp.kPrime = kPrimeMatrix_nd[1][1]; // 0-based (MATLAB: kPrimeMatrix_nd(2,2))

  // Air density from atmosphere model
  const atm = atmosphereModel(dp.altitude_m, dp.temperature_K);
  dp.airDensity_kgpm3 = atm.airDensity_kgpm3;

  // -------------------------------------------------------------------
  // 2. Generate Geometry
  // -------------------------------------------------------------------
  const hd = hullTopFunction(dp);
  const fd = finTopFunction(dp, hd);

  // -------------------------------------------------------------------
  // 3. Mass Properties
  // -------------------------------------------------------------------
  // Hull mass from skin surface area
  const areaPc = hullPointCloudDensity(hd.hullPointArray.x, hd.hullPointArray.r, dp.hullAngleSpaceing_rad);
  const hullCloudLen = hd.hullPointCloud.length;
  const hullPointMass = new Float64Array(hullCloudLen);
  let hullTotalMass_kg = 0;
  for (let i = 0; i < hullCloudLen; i++) {
    hullPointMass[i] = areaPc[i] * dp.hullMaterialMass_kgpm2;
    hullTotalMass_kg += hullPointMass[i];
  }

  // Fin mass from volume (scaled by external fraction)
  let singleFinVolume_m3 = finVolume(dp.thickNessRatio, fd.baseCord_m,
    dp.LeadingEdgeTapperAngle_rad, dp.tailingEdgeHight_m, dp.nacaCoeffs);
  singleFinVolume_m3 *= fd.externalFraction_nd;
  const singleFinMass_kg = singleFinVolume_m3 * dp.finMaterialDensity_kgpm3;

  // Fin centroids (lumped mass approach)
  const numFins = fd.numbFins;
  const finCentX = new Array(numFins);
  const finCentY = new Array(numFins);
  const finCentZ = new Array(numFins);
  for (let ii = 0; ii < numFins; ii++) {
    const fc = fd.finClouds[ii];
    const nPts = fc.x.length;
    let sx = 0, sy = 0, sz = 0;
    for (let j = 0; j < nPts; j++) {
      sx += fc.x[j]; sy += fc.y[j]; sz += fc.z[j];
    }
    finCentX[ii] = sx / nPts;
    finCentY[ii] = sy / nPts;
    finCentZ[ii] = sz / nPts;
  }

  // Combine hull point cloud with fin lumped masses
  const totalPts = hullCloudLen + numFins;
  const allX = new Float64Array(totalPts);
  const allY = new Float64Array(totalPts);
  const allZ = new Float64Array(totalPts);
  const allM = new Float64Array(totalPts);

  for (let i = 0; i < hullCloudLen; i++) {
    allX[i] = hd.hullPointCloud.x[i];
    allY[i] = hd.hullPointCloud.y[i];
    allZ[i] = hd.hullPointCloud.z[i];
    allM[i] = hullPointMass[i];
  }
  for (let i = 0; i < numFins; i++) {
    allX[hullCloudLen + i] = finCentX[i];
    allY[hullCloudLen + i] = finCentY[i];
    allZ[hullCloudLen + i] = finCentZ[i];
    allM[hullCloudLen + i] = singleFinMass_kg;
  }

  let totalMass_kg = 0;
  for (let i = 0; i < totalPts; i++) totalMass_kg += allM[i];

  // -------------------------------------------------------------------
  // 4. Pre-Ballast CG
  // -------------------------------------------------------------------
  const preBallastResult = computeInertiaMatrix(allX, allY, allZ, allM);
  const preBallastCG_m = preBallastResult.cg;
  const cvX = hd.centerOfVolume_m;

  // -------------------------------------------------------------------
  // 5. Ballast
  // -------------------------------------------------------------------
  let ballastMass_kg = 0;
  let ballastPos_m = [0, 0, 0];

  if (dp.ballastOption === 'nose') {
    ballastPos_m = [0, 0, 0];
    const xBallast = ballastPos_m[0];
    if (Math.abs(cvX - xBallast) > 1e-6) {
      ballastMass_kg = totalMass_kg * (preBallastCG_m[0] - cvX) / (cvX - xBallast);
    }
  } else if (dp.ballastOption === 'gondola') {
    ballastPos_m = [dp.noseLength_m, 0, dp.hullRadius_m];
    const xBallast = ballastPos_m[0];
    if (Math.abs(cvX - xBallast) > 1e-6) {
      ballastMass_kg = totalMass_kg * (preBallastCG_m[0] - cvX) / (cvX - xBallast);
    }
  } else if (dp.ballastOption !== 'none') {
    throw new Error("ballastOption must be 'none', 'nose', or 'gondola'");
  }

  ballastMass_kg = Math.max(0, ballastMass_kg);

  // Build final mass arrays (with ballast if applicable)
  let finalX, finalY, finalZ, finalM;
  if (ballastMass_kg > 0) {
    const n2 = totalPts + 1;
    finalX = new Float64Array(n2);
    finalY = new Float64Array(n2);
    finalZ = new Float64Array(n2);
    finalM = new Float64Array(n2);
    for (let i = 0; i < totalPts; i++) {
      finalX[i] = allX[i]; finalY[i] = allY[i]; finalZ[i] = allZ[i]; finalM[i] = allM[i];
    }
    finalX[totalPts] = ballastPos_m[0];
    finalY[totalPts] = ballastPos_m[1];
    finalZ[totalPts] = ballastPos_m[2];
    finalM[totalPts] = ballastMass_kg;
    totalMass_kg += ballastMass_kg;
  } else {
    finalX = allX; finalY = allY; finalZ = allZ; finalM = allM;
  }

  // -------------------------------------------------------------------
  // 6. Inertia Properties (with ballast)
  // -------------------------------------------------------------------
  const inertiaResult = computeInertiaMatrix(finalX, finalY, finalZ, finalM);
  const inertiaRotCg_kgm2 = inertiaResult.inertia;
  const centerOfGravity_m = inertiaResult.cg;

  // Virtual (apparent) mass inertia about CV
  const volumetricMomentCv_m5 = volumetricMoment(cvX, hd.hullPointArray.x, hd.hullPointArray.r);
  const virtualInertiaMatrixCv = scalarMulMat3(
    dp.airDensity_kgpm3 * volumetricMomentCv_m5,
    dp.vertialInertiaCoeff_nd
  );

  // -------------------------------------------------------------------
  // 7. Buoyancy
  // -------------------------------------------------------------------
  const buoy = computeBuoyancy(hd.hullVolume_m3, dp.airDensity_kgpm3,
    dp.liftingGasDensitySL_kgpm3, dp.temperature_K, dp.gravity_mps2);

  const totalWeight_N = totalMass_kg * dp.gravity_mps2;
  const excessLift_N = buoy.netLift_N - totalWeight_N;

  // -------------------------------------------------------------------
  // 8. Static Stability
  // -------------------------------------------------------------------
  const staticMargin_m = centerOfGravity_m[0] - cvX;
  const staticMoment_Nm = totalWeight_N * staticMargin_m;

  // -------------------------------------------------------------------
  // 9. Return results
  // -------------------------------------------------------------------
  return {
    // geometry
    hullVolume_m3:        hd.hullVolume_m3,
    centerOfVolume_m:     hd.centerOfVolume_m,
    hullPointArray:       hd.hullPointArray,
    hullPointCloud:       hd.hullPointCloud,
    finClouds:            fd.finClouds,
    baseCord_m:           fd.baseCord_m,
    requiredArea_m2:      fd.requiredArea_m2,
    externalFraction_nd:  fd.externalFraction_nd,
    numbFins:             fd.numbFins,
    // mass
    hullTotalMass_kg,
    singleFinVolume_m3,
    singleFinMass_kg,
    totalMass_kg,
    centerOfGravity_m,
    preBallastCG_m,
    ballastMass_kg,
    ballastPos_m,
    // inertia
    inertiaRotCg_kgm2,
    virtualInertiaMatrixCv,
    k1: dp.k1,
    k2: dp.k2,
    kPrime: dp.kPrime,
    // buoyancy
    airDensity_kgpm3:     dp.airDensity_kgpm3,
    gasDensity_kgpm3:     buoy.gasDensity_kgpm3,
    netLift_N:            buoy.netLift_N,
    buoyancyForce_N:      buoy.buoyancyForce_N,
    gasWeight_N:          buoy.gasWeight_N,
    totalWeight_N,
    excessLift_N,
    // stability
    staticMargin_m,
    staticMoment_Nm,
    // input echo
    liftingGas:           dp.liftingGas,
    ballastOption:        dp.ballastOption,
    altitude_m:           dp.altitude_m,
    temperature_K:        dp.temperature_K,
    hullLength_m:         dp.hullLength_m,
    hullRadius_m:         dp.hullRadius_m,
    tailingEdgeHight_m:   dp.tailingEdgeHight_m,
    dihedral_rad:         dp.dihedral_rad
  };
};

// Expose panel analysis separately (expensive, called on demand)
window.runPanelAnalysis = runPanelAnalysis;

// =========================================================================
// TEST HARNESS (commented out - uncomment to run verification)
// =========================================================================
/*
// Instructions: Run these 5 test cases through both the MATLAB code and
// this JavaScript translation. Compare all output values. Agreement should
// be within 1e-8 relative error on all scalar outputs.
//
// To run: uncomment this block and open the HTML page with a browser console.

const testCases = [
  {
    name: "Case 1: Default (helium, sea level, 3 fins, gondola ballast)",
    inputs: {}
  },
  {
    name: "Case 2: Hydrogen, sea level, 4 fins, nose ballast",
    inputs: {
      liftingGas: 'hydrogen',
      numbFins: 4,
      ballastOption: 'nose'
    }
  },
  {
    name: "Case 3: Helium, 2000m altitude, 300K, no ballast",
    inputs: {
      altitude_m: 2000,
      temperature_K: 300,
      ballastOption: 'none'
    }
  },
  {
    name: "Case 4: Large airship, helium, sea level",
    inputs: {
      hullLength_m: 25,
      hullRadius_m: 2.5,
      noseLength_m: 4,
      tailRatio: 5,
      xPosTailEdge_m: 23,
      tailingEdgeHight_m: 2,
      numbFins: 4,
      ballastOption: 'gondola'
    }
  },
  {
    name: "Case 5: Small blimp, hydrogen, high altitude",
    inputs: {
      hullLength_m: 8,
      hullRadius_m: 0.8,
      noseLength_m: 1.5,
      tailRatio: 2.5,
      xPosTailEdge_m: 7,
      tailingEdgeHight_m: 0.6,
      liftingGas: 'hydrogen',
      altitude_m: 5000,
      temperature_K: 255,
      numbFins: 3,
      ballastOption: 'nose',
      finMaterialDensity_kgpm3: 5
    }
  }
];

testCases.forEach(tc => {
  console.log('=== ' + tc.name + ' ===');
  const r = window.runBlimpDesign(tc.inputs);
  console.log('Hull volume:', r.hullVolume_m3, 'm^3');
  console.log('Centre of volume:', r.centerOfVolume_m, 'm');
  console.log('Hull mass:', r.hullTotalMass_kg, 'kg');
  console.log('Fin mass (each):', r.singleFinMass_kg, 'kg');
  console.log('Total mass:', r.totalMass_kg, 'kg');
  console.log('Ballast:', r.ballastMass_kg, 'kg');
  console.log('CG:', r.centerOfGravity_m);
  console.log('Net lift:', r.netLift_N, 'N');
  console.log('Excess lift:', r.excessLift_N, 'N');
  console.log('Static margin:', r.staticMargin_m, 'm');
  console.log('Static moment:', r.staticMoment_Nm, 'Nm');
  console.log('Ixx:', r.inertiaRotCg_kgm2[0][0],
              'Iyy:', r.inertiaRotCg_kgm2[1][1],
              'Izz:', r.inertiaRotCg_kgm2[2][2]);
  console.log('');
});
*/
