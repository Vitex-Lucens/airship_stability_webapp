function [c_l, c_d, Cp, Vtan, xc, yc, panelS] = vortexPanelMethod(Xb, Yb, Vinf, c, alpha_deg)
%VORTEXPANELMETHOD  Vortex panel method for 2-D airfoil aerodynamic coefficients.
%   Adapted from Vortex-Panel-Method-for-Airfoil-Aerodynamic-Coefficients.
%
%   Inputs:
%     Xb, Yb     - boundary point coordinates (Np x 1), ordered nose -> lower
%                   surface -> trailing edge -> upper surface -> nose (closed loop)
%     Vinf       - freestream velocity [m/s]
%     c          - chord length [m]
%     alpha_deg  - angle of attack [deg]
%
%   Outputs:
%     c_l    - sectional lift coefficient
%     c_d    - sectional drag coefficient
%     Cp     - pressure coefficient at each panel control point (M x 1)
%     Vtan   - tangential velocity ratio V/Vinf at each panel (M x 1)
%     xc, yc - control point coordinates (M x 1)
%     panelS - panel lengths (M x 1)

Mp = length(Xb);           % number of boundary points
M  = Mp - 1;               % number of panels
alpha = alpha_deg * pi / 180;

% Panel geometry: control points, lengths, angles
x     = zeros(1, M);
y     = zeros(1, M);
s     = zeros(1, M);
theta = zeros(1, M);
sine  = zeros(1, M);
cosine= zeros(1, M);
RHS   = zeros(1, Mp);

for i = 1:M
    ip = i + 1;
    x(i) = 0.5 * (Xb(i) + Xb(ip));
    y(i) = 0.5 * (Yb(i) + Yb(ip));
    s(i) = sqrt((Xb(ip) - Xb(i))^2 + (Yb(ip) - Yb(i))^2);
    theta(i)  = atan2(Yb(ip) - Yb(i), Xb(ip) - Xb(i));
    sine(i)   = sin(theta(i));
    cosine(i) = cos(theta(i));
    RHS(i)    = sin(theta(i) - alpha);
end

% Normal and tangential influence coefficient sub-matrices
CN1 = zeros(M, M);
CN2 = zeros(M, M);
CT1 = zeros(M, M);
CT2 = zeros(M, M);

for i = 1:M
    for j = 1:M
        if i == j
            CN1(i, j) = -1.0;
            CN2(i, j) =  1.0;
            CT1(i, j) =  0.5 * pi;
            CT2(i, j) =  0.5 * pi;
        else
            A = -(x(i) - Xb(j)) * cosine(j) - (y(i) - Yb(j)) * sine(j);
            B = (x(i) - Xb(j))^2 + (y(i) - Yb(j))^2;
            C = sin(theta(i) - theta(j));
            D = cos(theta(i) - theta(j));
            E = (x(i) - Xb(j)) * sine(j) - (y(i) - Yb(j)) * cosine(j);
            F = log(1.0 + s(j) * (s(j) + 2 * A) / B);
            G = atan2(E * s(j), B + A * s(j));
            P = (x(i) - Xb(j)) * sin(theta(i) - 2 * theta(j)) ...
              + (y(i) - Yb(j)) * cos(theta(i) - 2 * theta(j));
            Q = (x(i) - Xb(j)) * cos(theta(i) - 2 * theta(j)) ...
              - (y(i) - Yb(j)) * sin(theta(i) - 2 * theta(j));
            CN2(i, j) = D + 0.5 * Q * F / s(j) - (A * C + D * E) * G / s(j);
            CN1(i, j) = 0.5 * D * F + C * G - CN2(i, j);
            CT2(i, j) = C + 0.5 * P * F / s(j) + (A * D - C * E) * G / s(j);
            CT1(i, j) = 0.5 * C * F - D * G - CT2(i, j);
        end
    end
end

% Assemble influence coefficient matrix AN and tangential matrix AT
AN = zeros(Mp, Mp);
AT = zeros(M, Mp);

for i = 1:M
    AN(i, 1)  = CN1(i, 1);
    AN(i, Mp) = CN2(i, M);
    AT(i, 1)  = CT1(i, 1);
    AT(i, Mp) = CT2(i, M);
    for j = 2:M
        AN(i, j) = CN1(i, j) + CN2(i, j - 1);
        AT(i, j) = CT1(i, j) + CT2(i, j - 1);
    end
end

% Kutta condition: gamma(1) + gamma(Mp) = 0
AN(Mp, 1)  = 1.0;
AN(Mp, Mp) = 1.0;
for j = 2:M
    AN(Mp, j) = 0.0;
end
RHS(Mp) = 0.0;

% Solve for vortex strengths
Gama = AN \ RHS';

% Tangential velocity and pressure coefficient at each panel
V  = zeros(1, M);
Cp_arr = zeros(1, M);
for i = 1:M
    V(i) = cos(theta(i) - alpha);
    for j = 1:Mp
        V(i) = V(i) + AT(i, j) * Gama(j);
    end
    Cp_arr(i) = 1.0 - V(i)^2;
end

% Force coefficients from pressure integration
cy = 0.0;
cx = 0.0;
for i = 1:M
    cy = cy - Cp_arr(i) * s(i) * cosine(i);
    cx = cx - Cp_arr(i) * s(i) * (-sine(i));
end
cy = cy / c;
cx = cx / c;

c_d =  cx * cos(alpha) + cy * sin(alpha);
c_l = -cx * sin(alpha) + cy * cos(alpha);

% Outputs
Cp     = Cp_arr(:);
Vtan   = V(:);
xc     = x(:);
yc     = y(:);
panelS = s(:);
end
