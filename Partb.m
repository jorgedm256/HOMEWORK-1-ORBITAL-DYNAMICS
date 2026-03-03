clear
clc

%% Part a)    This part is correct. Exercise IOD1 was very similar, so we used it
muE=3.986E5;
t=13.5;
rS=0.8;
phi=50;
rP_2_SEZ=[5369.09; 2332.14; 656.480];
vP_2_SEZ=[-3.37027; 4.90364; 0.106703];
r1=rP_2_SEZ;
v1=vP_2_SEZ;

[r0,v0] = SEZ2ECI(r1,v1,rS,phi,t);

% --- SHOW RESULTS ---
fprintf('\n--- Part (a) Results ---\n');

fprintf('r0:%.6f \n', r0);
fprintf('v0:%.6f\n', v0);



%% Part b) Verify
%% COE

r_vec = r0; 
v_vec = v0;

[a, e, i_rad, Omega_rad, omega_rad, nu_rad] = rv2coe(r_vec, v_vec, muE);
COE_part_b=[a, e, i_rad*180/pi, Omega_rad*180/pi, omega_rad*180/pi, nu_rad*180/pi];

if COE_part_b(3) > 180
    COE_part_b(3)=COE_part_b(3)-360;
end
if COE_part_b(4) > 180
    COE_part_b(4)=COE_part_b(4)-360;
end
if COE_part_b(5) > 180
    COE_part_b(5)=COE_part_b(5)-360;
end
if COE_part_b(6) > 180
    COE_part_b(6)=COE_part_b(6)-360;
end



% --- SHOW RESULTS ---
fprintf('\n--- Part (b) Results ---\n');
fprintf('\n--- Classical orbital elements (COEs) ---\n');
fprintf('Semi-major axis (a):            %.4f km\n', COE_part_b(1));
fprintf('Eccentricity (e):            %.6f\n', COE_part_b(2));
fprintf('Inclination (i):              %.4f deg\n', COE_part_b(3));
fprintf('RAAN (Omega):                 %.4f deg\n', COE_part_b(4));
fprintf('Argument of periapsis (w):    %.4f deg\n', COE_part_b(5));
fprintf('True anomaly (nu):      %.4f deg\n', COE_part_b(6));


% --- PLOTTING THE ORBIT ---

% 1. Create an array of true anomalies to draw the full ellipse
theta_plot = linspace(0, 2*pi, 200);

% 2. Calculate the distance (r) for each true anomaly
p = a * (1 - e^2); % Semi-latus rectum
r_plot = p ./ (1 + e * cos(theta_plot));

% 3. Coordinates in the Perifocal Frame (P, Q, W)
X_peri = r_plot .* cos(theta_plot);
Y_peri = r_plot .* sin(theta_plot);
Z_peri = zeros(1, length(theta_plot)); % Orbit lies flat in the perifocal frame

% 4. Rotation Matrix from Perifocal to ECI (R3(-Omega) * R1(-i) * R3(-omega))
cw = cos(omega_rad); sw = sin(omega_rad);
cO = cos(Omega_rad); sO = sin(Omega_rad);
ci = cos(i_rad);     si = sin(i_rad);

Q_pX = [cO*cw - sO*ci*sw,  -cO*sw - sO*ci*cw,  sO*si;
        sO*cw + cO*ci*sw,  -sO*sw + cO*ci*cw, -cO*si;
        si*sw,              si*cw,             ci];

% 5. Rotate the perifocal coordinates into the ECI frame
r_ECI_orbit = Q_pX * [X_peri; Y_peri; Z_peri];

% 6. Plotting
figure;
hold on; grid on; axis equal;

% Plot the Earth as a sphere for reference
[X_E, Y_E, Z_E] = sphere(50);
surf(X_E*6371, Y_E*6371, Z_E*6371, 'EdgeColor', 'none', 'FaceColor', 'c', 'FaceAlpha', 0.2);

% Plot the orbit trajectory
plot3(r_ECI_orbit(1,:), r_ECI_orbit(2,:), r_ECI_orbit(3,:), 'b-', 'LineWidth', 1.5);

% Plot the spacecraft's current position at the observation time
plot3(r0(1), r0(2), r0(3), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');

% Formatting
xlabel('X_{ECI} (km)'); ylabel('Y_{ECI} (km)'); zlabel('Z_{ECI} (km)');
title('Spacecraft Orbit in ECI Frame');
legend('Earth', 'Orbit Trajectory', 'Current S/C Position');
view(3);
hold off;




%% --- PART (C) START --- Verify
fprintf('\n--- Part (c) Results ---\n');

% -------------------------------------------------------------------------
% REQUIREMENT 1: Local time when the spacecraft reaches theta = theta_0 + pi/2
% -------------------------------------------------------------------------

% 1.1 Calculate current Eccentric Anomaly (E_0) and new true anomaly (nu_1)
nu_0 = nu_rad; % From part b
nu_1 = nu_0 + pi/2;

% Using your code for 2BP13
ne=sqrt(muE/a^3);
E0=theta2E(nu_0,e);
M_0=E0-e*sin(E0);
tfromP=M_0/ne;  % this is the time that passes from passing from the pericenter to theta=nu_0

theta = nu_1;
E = theta2E(theta, e);
M = E - e * sin(E);
tfromPdef = M/ne-tfromP;  % time from theta0 to current theta
delta_t1_sec=tfromPdef;
delta_t1_hr = delta_t1_sec / 3600;

% % This alternative method should give you the same result
% % Robust way to calculate E_0 using atan2 to avoid quadrant issues
% sin_E0 = (sqrt(1-e^2)*sin(nu_0)) / (1 + e*cos(nu_0));
% cos_E0 = (e + cos(nu_0)) / (1 + e*cos(nu_0));
% E_0 = atan2(sin_E0, cos_E0);
% 
% % Robust way to calculate E_1
% sin_E1 = (sqrt(1-e^2)*sin(nu_1)) / (1 + e*cos(nu_1));
% cos_E1 = (e + cos(nu_1)) / (1 + e*cos(nu_1));
% E_1 = atan2(sin_E1, cos_E1);
% 
% % Ensure E_1 is greater than E_0 for forward time propagation
% if E_1 < E_0
%     E_1 = E_1 + 2*pi;
% end
% 
% % 1.2 Calculate Mean Anomalies (M_0 and M_1) using Kepler's Equation
% M_0 = E_0 - e*sin(E_0);
% M_1 = E_1 - e*sin(E_1);
% 
% % 1.3 Calculate Mean Motion (n) and Time of Flight (delta_t)
% n = sqrt(muE / a^3); % Mean motion in rad/s
% delta_t1_sec = (M_1 - M_0) / n; % Time of flight in seconds
% delta_t1_hr = delta_t1_sec / 3600;

% 1.4 Calculate new local time
local_time_initial = 22.5; % Initial observation at 22:30
local_time_new_hr = mod(local_time_initial + delta_t1_hr, 24);

% Convert decimal hours to HH:MM:SS format for displaying
hrs = floor(local_time_new_hr);
mins = floor((local_time_new_hr - hrs) * 60);
secs = round((local_time_new_hr - hrs - mins/60) * 3600);

fprintf('1. Local time at theta = theta_0 + pi/2 is: %02d:%02d:%02d (next day if past midnight)\n', hrs, mins, secs);

% -------------------------------------------------------------------------
% REQUIREMENT 2: Position and velocity in ECI 6 hours later
% -------------------------------------------------------------------------

% This part was done using 2BP12
t0=M_0/ne;  % t0 is the time that passes from passing from the pericenter to theta0

t1=6*3600+t0;
%%% This was done because we are not starting at theta =0
M1=ne*t1;
Eguess=M1;
fE= @(E) E - e*sin(E) - ne*t1;
dfE= @(E) 1-e*cos(E);
tol=10^-8;
maxiter=100;

Eattime=newton(fE,dfE,Eguess,tol,maxiter);
nu_2=E2theta(Eattime,e);

% This alternative method should give you the same result
% delta_t2_sec = 6 * 3600; % 6 hours in seconds
% 
% % 2.1 Calculate new Mean Anomaly 6 hours later
% M_2 = M_0 + ne * delta_t2_sec;
% M_2 = mod(M_2, 2*pi); % Keep within 0 to 2*pi
% 
% % 2.2 Solve Kepler's Equation (M = E - e*sin(E)) for E_2 using Newton-Raphson
% E_guess = M_2;
% tol = 1e-8;
% err = 1;
% while err > tol
%     f = E_guess - e*sin(E_guess) - M_2;
%     fp = 1 - e*cos(E_guess);
%     E_new = E_guess - f/fp;
%     err = abs(E_new - E_guess);
%     E_guess = E_new;
% end
% E_2 = E_guess;
% 
% % 2.3 Calculate new true anomaly (nu_2)
% sin_nu2 = (sqrt(1-e^2)*sin(E_2)) / (1 - e*cos(E_2));
% cos_nu2 = (cos(E_2) - e) / (1 - e*cos(E_2));
% nu_2 = atan2(sin_nu2, cos_nu2);

% 2.4 Use Orbital Elements to get ECI state (using helper function at bottom)
p_semi = a * (1 - e^2); % Semi-latus rectum
[r_ECI_6h, v_ECI_6h] = coe2rv(muE, p_semi, e, i_rad, Omega_rad, omega_rad, nu_2);

fprintf('2. State in ECI 6 hours later:\n');
fprintf('   r_ECI = [%.3f, %.3f, %.3f] km\n', r_ECI_6h(1), r_ECI_6h(2), r_ECI_6h(3));
fprintf('   v_ECI = [%.5f, %.5f, %.5f] km/s\n', v_ECI_6h(1), v_ECI_6h(2), v_ECI_6h(3));

% -------------------------------------------------------------------------
% REQUIREMENT 4: Compute state in SEZ after 6 hours
% (Calls the function from Requirement 3)
% -------------------------------------------------------------------------

% The time 't' since Oxz crossing is initial t (13.5h) + 6 hours
t_6h = t + 6; % In hours

[r_SEZ_6h, v_SEZ_6h] = ECI2SEZ(r_ECI_6h, v_ECI_6h, rS, phi, t_6h);

fprintf('3. State in SEZ 6 hours later:\n');
fprintf('   r_SEZ = [%.3f, %.3f, %.3f] km\n', r_SEZ_6h(1), r_SEZ_6h(2), r_SEZ_6h(3));
fprintf('   v_SEZ = [%.5f, %.5f, %.5f] km/s\n', v_SEZ_6h(1), v_SEZ_6h(2), v_SEZ_6h(3));





%% --- PART (D) START --- Verify
fprintf('\n--- Part (d) Results ---\n');

% 1. Define masses (in kg)
m_total = 2000;      % Initial spacecraft mass (from problem description)
m_payload = 500;     % Payload mass
m_sc_new = 1500;     % Remaining spacecraft mass (2000 - 500)

% 2. Define payload velocity in ECI (ensure it's a column vector)
v_payload = [-6.0; 3.0; -1.6]; % km/s

% 3. Apply Conservation of Linear Momentum to find new S/C velocity
% Ensure v_ECI_6h is a column vector to prevent matrix dimension errors
if size(v_ECI_6h, 2) > 1
    v_ECI_6h = v_ECI_6h';
end

v_sc_new = (m_total * v_ECI_6h - m_payload * v_payload) / m_sc_new;

% Position remains identical instantaneously
r_sc_new = r_ECI_6h;
if size(r_sc_new, 2) > 1
    r_sc_new = r_sc_new';
end

% 4. Compute new Classical Orbital Elements using rv2coe
[a_new, e_new, i_new_rad, Omega_new_rad, omega_new_rad, nu_new_rad] = rv2coe(r_sc_new, v_sc_new, muE);

% Convert angles to degrees for display
i_new_deg = i_new_rad * 180/pi;
Omega_new_deg = Omega_new_rad * 180/pi;
omega_new_deg = omega_new_rad * 180/pi;
nu_new_deg = nu_new_rad * 180/pi;

if i_new_deg >180
    i_new_deg=i_new_deg-360;
end

if Omega_new_deg >180
    Omega_new_deg=Omega_new_deg-360;
end
if omega_new_deg >180
    omega_new_deg=omega_new_deg-360;
end
if nu_new_deg >180
    nu_new_deg=nu_new_deg-360;
end

% --- SHOW RESULTS ---
fprintf('New Spacecraft Velocity after release:\n');
fprintf('   v_sc_new = [%.5f, %.5f, %.5f] km/s\n', v_sc_new(1), v_sc_new(2), v_sc_new(3));

fprintf('\n--- New Classical Orbital Elements (Post-Release) ---\n');
fprintf('Semi-major axis (a):          %.4f km\n', a_new);
fprintf('Eccentricity (e):             %.6f\n', e_new);
fprintf('Inclination (i):              %.4f deg\n', i_new_deg);
fprintf('RAAN (Omega):                 %.4f deg\n', Omega_new_deg);
fprintf('Argument of periapsis (w):    %.4f deg\n', omega_new_deg);
fprintf('True anomaly (nu):            %.4f deg\n', nu_new_deg);


%% --- PART (E) START --- Verify
fprintf('\n--- Part (e) Results ---\n');

% 1. Get current radius magnitude (position doesn't change during maneuver)
r_mag = norm(r_sc_new);

% 2. Calculate the required magnitude of the circular velocity
v_circ_mag = sqrt(muE / r_mag);

% 3. Determine the direction of the new velocity
% The orbital plane stays the same, so the specific angular momentum 
% direction (h_hat) remains unchanged.
h_vec_new = cross(r_sc_new, v_sc_new);
h_hat = h_vec_new / norm(h_vec_new);

r_hat = r_sc_new / r_mag;

% The velocity in a circular orbit is perfectly perpendicular to 'r' 
% and lies in the orbital plane (perpendicular to 'h')
v_circ_vec = v_circ_mag * cross(h_hat, r_hat);

% 4. Calculate the Delta-V vector required for the maneuver
delta_v_vec = v_circ_vec - v_sc_new;

% Calculate the total Delta-V magnitude (cost of the maneuver)
delta_v_mag = norm(delta_v_vec);

% --- MOSTRAR RESULTADOS EN CONSOLA ---
fprintf('Target circular velocity vector:\n');
fprintf('   v_circ = [%.5f, %.5f, %.5f] km/s\n', v_circ_vec(1), v_circ_vec(2), v_circ_vec(3));
fprintf('Required Delta-V vector:\n');
fprintf('   Delta_V = [%.5f, %.5f, %.5f] km/s\n', delta_v_vec(1), delta_v_vec(2), delta_v_vec(3));
fprintf('Total Delta-V magnitude:      %.4f km/s\n', delta_v_mag);