clear

clc
%%
%[text] # **Problem (\#1)**
% *Student 1:* 100591752

%

% *Student 2:* 100495586
%%
%[text] ## Data of the problem:
%[text] 
%[text] \-Spacecraft mass = 2000 kg
%[text] \-Ground station location: 50 degrees North and 800 m above sea level
%[text] \-At 22:30, the state of the aircraft with respect to the SEZ reference frame is:
%[text] r1=\[5369.09; 2332.14; 656.480\];
%[text] v2=\[-3.37027; 4.90364; 0.106703\];
%[text] \-The station crossed the Oxz plane of the ECI reference frame at 9:00
%%
%[text] ## a)

muE=3.986E5;


% The time passed from 9:00 to 22:30 in hours is 13:30 hours:


t=13.5;

rS=0.8;


% The latitude is passed to radians:

phi=50*pi/180;

% The position r1 and velocity v1:

rP_2_SEZ=[5369.09; 2332.14; 656.480];

vP_2_SEZ=[-3.37027; 4.90364; 0.106703];

r1=rP_2_SEZ;

v1=vP_2_SEZ;


% SEZ2ECI function is used to get the position and velocity at the ECI reference

% frame:


[r0,v0] = SEZ2ECI(r1,v1,rS,phi,t);

fprintf('\n--- Part (a) Results ---\n'); %[output:8ef14807]

fprintf(' r0 = [%.5f, %.5f, %.5f] km\n', r0(1), r0(2), r0(3)); %[output:1d16c762]

fprintf(' v0 = [%.5f, %.5f, %.5f] km/s\n', v0(1), v0(2), v0(3)); %[output:6bf0c1da]
%%
%[text] ## b)
r_vec = r0;

v_vec = v0;


% Having the position and velocity, rv2coe function is used to get the classical

% orbital elements:

[a, e, i_rad, Omega_rad, omega_rad, nu_rad] = rv2coe(r_vec, v_vec, muE);

COE_part_b=[a, e, i_rad*180/pi, Omega_rad*180/pi, omega_rad*180/pi, nu_rad*180/pi];


% In order to show all the angles between -180 and 180 degrees:

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

fprintf('\n--- Part (b) Results ---\n'); %[output:5a595f0c]

fprintf('\n--- Classical orbital elements (COEs) ---\n'); %[output:0e562bd3]

fprintf('Semi-major axis (a): %.4f km\n', COE_part_b(1)); %[output:8ed01f17]

fprintf('Eccentricity (e): %.6f\n', COE_part_b(2)); %[output:3ffad306]

fprintf('Inclination (i): %.4f deg\n', COE_part_b(3)); %[output:4dec4604]

fprintf('RAAN (Omega): %.4f deg\n', COE_part_b(4)); %[output:7f597f7a]

fprintf('Argument of periapsis (w): %.4f deg\n', COE_part_b(5)); %[output:0beaa84d]

fprintf('True anomaly (nu): %.4f deg\n', COE_part_b(6)); %[output:0d7ad85d]

% In order to plot the orbit, first an array of true anomalies is made:


theta_plot = linspace(0, 2*pi, 200);

% The distance (r) is calculated for each true anomaly:

p = a * (1 - e^2); % Semi-latus rectum

r_plot = p ./ (1 + e * cos(theta_plot));

% The coordinates are drawn into the perifocal frame

X_peri = r_plot .* cos(theta_plot);

Y_peri = r_plot .* sin(theta_plot);

Z_peri = zeros(1, length(theta_plot)); % Orbit lies flat in the perifocal frame

% The rotation matrix from perifocal to ECI

% The rotaion is (R3(-Omega) * R1(-i) *R3(-omega)):

cw = cos(omega_rad); sw = sin(omega_rad);

cO = cos(Omega_rad); sO = sin(Omega_rad);

ci = cos(i_rad); si = sin(i_rad);

Q_pX = [cO*cw - sO*ci*sw, -cO*sw - sO*ci*cw, sO*si;

sO*cw + cO*ci*sw, -sO*sw + cO*ci*cw, -cO*si;

si*sw, si*cw, ci];

% The perifocal coordinates are rotated into the ECI reference frame:

r_ECI_orbit = Q_pX * [X_peri; Y_peri; Z_peri];

% Plotting

figure(); %[output:1f16fe28]

hold on; grid on; axis equal; %[output:1f16fe28]

% The Earth is plotted as a sphere for reference:

[X_E, Y_E, Z_E] = sphere(50);

surf(X_E*6371, Y_E*6371, Z_E*6371, 'EdgeColor', 'none', 'FaceColor', 'c', 'FaceAlpha', 0.2); %[output:1f16fe28]


% The following curve shows the trajectory:

plot3(r_ECI_orbit(1,:), r_ECI_orbit(2,:), r_ECI_orbit(3,:), 'b-', 'LineWidth', 1.5); %[output:1f16fe28]

% The following shows the spacecraft's position at the observation time:

plot3(r0(1), r0(2), r0(3), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r'); %[output:1f16fe28]

xlabel('X_{ECI} (km)'); ylabel('Y_{ECI} (km)'); zlabel('Z_{ECI} (km)'); %[output:1f16fe28]

title('Spacecraft Orbit in ECI Frame'); %[output:1f16fe28]

legend('Earth', 'Orbit Trajectory', 'Current S/C Position'); %[output:1f16fe28]

view(3); %[output:1f16fe28]

hold off; %[output:1f16fe28]
%%
%[text] ## **c)**
% c.1) Local time when the spacecraft reaches theta = theta_0 + pi/2

% The current eccentric anomaly (E_0) and new true anomaly (nu_1) are first

% calculated:

nu_0 = nu_rad; % From part b

nu_1 = nu_0 + pi/2;

% Because nu_0 is not equal to 0, the time that has passed from the orbit pericenter

% to nu_0 needs to be calculated, named as tfromP:

ne=sqrt(muE/a^3);

E0=theta2E(nu_0,e);

M_0=E0-e*sin(E0);

tfromP=M_0/ne;

% Now the time that has passed from nu_0 to nu_0 + pi/2 is calculated:

theta = nu_1;

E = theta2E(theta, e);

M = E - e * sin(E);

tfromPdef = M/ne-tfromP; % time from theta0 to current theta (time from theta0 to theta minus time from pericenter to theta0)

delta_t1_sec=tfromPdef;

delta_t1_hr = delta_t1_sec / 3600;

% The new local time is calculated:

local_time_initial = 22.5; % Initial observation at 22:30

local_time_new_hr = mod(local_time_initial + delta_t1_hr, 24);

% Time is converted to decimal hours to HH:MM:SS format for displaying

hrs = floor(local_time_new_hr);

mins = floor((local_time_new_hr - hrs) * 60);

secs = round((local_time_new_hr - hrs - mins/60) * 3600);
fprintf('\n--- Part (c) Results ---\n'); %[output:4fdbba66]

fprintf('1. Local time at theta = theta_0 + pi/2 is: %02d:%02d:%02d (next day if past midnight)\n', hrs, mins, secs); %[output:91609af9]

% c.2) Position and velocity in ECI 6 hours later

% First the time that has passed from crossing the pericenter to nu_0 (initial

% observation) is calculated:

t0=M_0/ne;

% The 6 hours are added to the initial time:

t1=6*3600+t0;

% Newton method is used to iterate the E:

M1=ne*t1;

Eguess=M1;

fE= @(E) E - e*sin(E) - ne*t1;

dfE= @(E) 1-e*cos(E);

tol=10^-8;

maxiter=100;

Eattime=newton(fE,dfE,Eguess,tol,maxiter); %[output:60d354e5]

nu_2=E2theta(Eattime,e);

% Now that the true anomaly has changed, the position and velocity need to be

% recalculated using the new orbital elements (true anomaly):

p_semi = a * (1 - e^2); % Semi-latus rectum

[r_ECI_6h, v_ECI_6h] = coe2rv(muE, p_semi, e, i_rad, Omega_rad, omega_rad, nu_2);

fprintf('\n 2. State in ECI 6 hours later:\n'); %[output:3e42949e]

fprintf(' r0 = [%.3f, %.3f, %.3f] km\n', r_ECI_6h(1), r_ECI_6h(2), r_ECI_6h(3)); %[output:3b1afd5c]

fprintf(' v0 = [%.5f, %.5f, %.5f] km/s\n', v_ECI_6h(1), v_ECI_6h(2), v_ECI_6h(3)); %[output:936e1791]

% c.3) Compute state in SEZ after 6 hours

% The time that has passed since Oxz crossing is initial t (13.5h) + 6 hours:

t_6h = t + 6; % In hours

t=t_6h;

r0=r_ECI_6h;

v0=v_ECI_6h;

% Applying ECI2SEZ function:

[r1, v1] = ECI2SEZ(r0, v0, rS, phi, t);

r_SEZ_6h=r1;

v_SEZ_6h=v1;

fprintf('3. State in SEZ 6 hours later:\n'); %[output:6b30419d]

fprintf(' r1 = [%.3f, %.3f, %.3f] km\n', r_SEZ_6h(1), r_SEZ_6h(2), r_SEZ_6h(3)); %[output:9e4e0fd0]

fprintf(' v1 = [%.5f, %.5f, %.5f] km/s\n', v_SEZ_6h(1), v_SEZ_6h(2), v_SEZ_6h(3)); %[output:78f51728]
%%
%[text] ## d)
% The initial mass is 2000 kg:

m = 2000; % Initial spacecraft mass

% The payload mass is 500 kg:

m_payload = 500; % Payload mass

% Then the remaining aircraft mass is:

m_sc_new = 1500; % Remaining spacecraft mass (2000 - 500)

% The payload velocity in ECI is:

v_payload = [-6.0; 3.0; -1.6]; % km/s

% Ensure v_ECI_6h is a column vector to prevent matrix dimension errors

if size(v_ECI_6h, 2) > 1

v_ECI_6h = v_ECI_6h';

end

% Conservation of linear momentum is applied to find the new spacecraft velocity:

v_sc_new = (m * v_ECI_6h - m_payload * v_payload) / m_sc_new;

% Position remains identical instantaneously because the separation takes place

% in a very small amount of time.

r_sc_new = r_ECI_6h;

if size(r_sc_new, 2) > 1

r_sc_new = r_sc_new';

end

% Because the velocity and position have changed, the new classical orbital

% elements are calculated using rv2coe:

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
fprintf('\n--- Part (d) Results ---\n'); %[output:6e1601cf]

fprintf('New Spacecraft Velocity after release:\n'); %[output:9d111e2f]

fprintf(' v_sc_new = [%.5f, %.5f, %.5f] km/s\n', v_sc_new(1), v_sc_new(2), v_sc_new(3)); %[output:23b00857]

fprintf('\n--- New Classical Orbital Elements (Post-Release) ---\n'); %[output:65f0f99f]

fprintf('Semi-major axis (a): %.4f km\n', a_new); %[output:2405e50f]

fprintf('Eccentricity (e): %.6f\n', e_new); %[output:9771e70c]

fprintf('Inclination (i): %.4f deg\n', i_new_deg); %[output:4334c999]

fprintf('RAAN (Omega): %.4f deg\n', Omega_new_deg); %[output:67339525]

fprintf('Argument of periapsis (w): %.4f deg\n', omega_new_deg); %[output:25b18ce5]

fprintf('True anomaly (nu): %.4f deg\n', nu_new_deg); %[output:2b2c80a5]
%%
%[text] ## e)
% Using the semi-major axis (a_new) and eccentricity (e_new) from part (d),

% the new radius of pericenter is obtained:

r_p = a_new * (1 - e_new);

% Using the Vis-Viva equation, the velocity at the pericenter is found:

v_p_elliptical = sqrt(muE * (2/r_p - 1/a_new));

% The velocity required for a circular orbit with that same radius is calculated:

% For a circular orbit, v = sqrt(mu/r)

v_p_circular = sqrt(muE / r_p);

% Finally, it is computed the delta-V required to circularize the orbit. Since

% both velocities are strictly tangential at pericenter, we just subtract their

% magnitudes.

delta_V = abs(v_p_circular - v_p_elliptical);

% --- SHOW RESULTS ---
fprintf('\n--- Part (e) Results ---\n'); %[output:4c056766]

fprintf('Radius of new pericenter (r_p): %.4f km\n', r_p); %[output:044c00ed]

fprintf('Velocity at pericenter (current): %.5f km/s\n', v_p_elliptical); %[output:17d4a4ba]

fprintf('Velocity for circular orbit: %.5f km/s\n', v_p_circular); %[output:32a6778c]

fprintf('Required Delta-V to circularize: %.5f km/s\n', delta_V); %[output:9d4f7384]
%%
%[text] ## Discussion
%[text] This exercise successfully modeled the trajectory, time propagation, and impulsive maneuvers of a 2000 kg spacecraft in Earth orbit.
%[text] **Initial Orbit (Parts a & b):** By transforming the topocentric SEZ observation into the inertial ECI frame, we determined the spacecraft's initial orbit. The resulting Classical Orbital Elements reveal a relatively low-eccentricity orbit ($e \\approx 0.1001$) with a semi-major axis of $a \\approx 8995.6$ km and an inclination of $25^\\circ$. As verified by the 3D plot, the orbit is elliptical but visually quite close to circular.
%[text] **Time Propagation (Part c):** Using Kepler's equation, we propagated the spacecraft's state 6 hours forward. Because the eccentricity is small, the Newton-Raphson numerical method was highly efficient, converging on the new eccentric anomaly in just 2 iterations. This gave us the exact new state vectors and local time (23:02:03) at the required future true anomaly.
%[text] **Payload Separation (Part d):** The instantaneous release of a 500 kg payload was calculated using the conservation of linear momentum. Because the payload was ejected with a specific velocity vector, the remaining 1500 kg spacecraft experienced a reactive change in velocity. This maneuver fundamentally altered its orbit, decreasing both its semi-major axis (to $8778.7$ km) and its eccentricity (dropping to $0.0613$).
%[text] **Circularization Maneuver (Part e):** Finally, we calculated the cost to circularize the new orbit at its pericenter ($r\_p \\approx 8240.8$ km). Because a maneuver at an apsidal point only requires a change in velocity magnitude (not direction), we used the scalar difference from the Vis-Viva equation. The current velocity at pericenter was $7.164$ km/s, while a circular orbit required $6.955$ km/s. The required $\\Delta V$ is therefore $\\approx 0.2099$ km/s. This relatively small required "burn" makes perfect physical sense, as the post-release orbit was already highly circular ($e \\approx 0.061$).
%[text] ## AI Declaration
%[text] We declare that no AI has been used for the development of the code nor discussion of the results of this homework.

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":45.1}
%---
%[output:8ef14807]
%   data: {"dataType":"text","outputData":{"text":"\n--- Part (a) Results ---\n","truncated":false}}
%---
%[output:1d16c762]
%   data: {"dataType":"text","outputData":{"text":" r0 = [-7081.21219, -5457.42427, 1932.79031] km\n","truncated":false}}
%---
%[output:6bf0c1da]
%   data: {"dataType":"text","outputData":{"text":" v0 = [4.59530, -4.08358, 2.24811] km\/s\n","truncated":false}}
%---
%[output:5a595f0c]
%   data: {"dataType":"text","outputData":{"text":"\n--- Part (b) Results ---\n","truncated":false}}
%---
%[output:0e562bd3]
%   data: {"dataType":"text","outputData":{"text":"\n--- Classical orbital elements (COEs) ---\n","truncated":false}}
%---
%[output:8ed01f17]
%   data: {"dataType":"text","outputData":{"text":"Semi-major axis (a): 8995.6258 km\n","truncated":false}}
%---
%[output:3ffad306]
%   data: {"dataType":"text","outputData":{"text":"Eccentricity (e): 0.100106\n","truncated":false}}
%---
%[output:4dec4604]
%   data: {"dataType":"text","outputData":{"text":"Inclination (i): 25.0051 deg\n","truncated":false}}
%---
%[output:7f597f7a]
%   data: {"dataType":"text","outputData":{"text":"RAAN (Omega): -169.9931 deg\n","truncated":false}}
%---
%[output:0beaa84d]
%   data: {"dataType":"text","outputData":{"text":"Argument of periapsis (w): 135.2704 deg\n","truncated":false}}
%---
%[output:0d7ad85d]
%   data: {"dataType":"text","outputData":{"text":"True anomaly (nu): -105.2767 deg\n","truncated":false}}
%---
%[output:1f16fe28]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAArwAAAGlCAYAAAAYiyWNAAAAAXNSR0IArs4c6QAAIABJREFUeF7tvQvwH1V5\/3\/CPcQEDARMiBjFBMXWqPGSwVp+1kunl9CqtYmZVhsjYkeDiIEkWlFASYCgVawaMaZazYRe7N\/EvzMt0v4YnAz5a9QIrS3UEDEkCgQTAoSU23\/ei+eT892c3T3X3XN23zuTgeSze\/Y5r+fs7nuffc5zxj311FNPCW4kQAIkQAIkQAIkQAIk0FMC4yh4e+pZdosESIAESIAESIAESKAgQMHLgUACJEACJEACJEACJNBrAhS8vXYvO0cCJEACJEACJEACJEDByzFAAiRAAiRAAiRAAiTQawIUvL12LztHAiRAAiRAAiRAAiRAwcsxQAIkQAIkQAIkQAIk0GsCFLy9di87RwIkQAIkQAIkQAIkQMHLMUACJEACJEACJEACJNBrAhS8vXYvO0cCJEACJEACJEACJEDByzFAAiRAAiRAAiRAAiTQawIUvL12LztHAiRAAiRAAiRAAiRAwcsxQAKZEDhw4IBYvny52LRp0xiLV65cKebPn59JL7o1U8dw\/fr1Yu7cueKGG24QmzdvFqtWrRLjx49vNLTKH+eff75YtmxZ4\/HqDldddZVYs2aNmDdvXuP577zzTrFo0SKxa9cuIW2XbZn0QbXbxdaqjt16661i4cKFtf2eNm2aWLdunZg5c+aY\/XTH6lg88MADYvHixWLbtm2izvYq36gnnT17tli7dq2YPHmyla+4MwmQQJ4EKHjz9ButHhgBVeTouh5SuPQZLQThihUrRl2UAuwb3\/iGseDEwSbirixG67iGELymbaQkeJuEaVkgU\/D2+epk30ggLgEK3rh82ToJBCEgxUxdYzYCK4hRGTYiBW85umcqFtHlppcPiaUqmqnDZnP+Kuwh2vBxqclLQJmJybhWfUXB6+MhHksCwyZAwTts\/7P3GRBQo2Dlz7yq+FKjvKr4+aM\/+iPxrne9a9RTXQqEKiTkjlWffHWCTye2y\/vVfUIuC59yP9X+vOQlLxFXXHFFYabsS5P9VZHEN77xjWL37t3itttuGzMS6iLmqq1lO8vnqfJJuQ87duwYRZjf+c53ive+971FyoLaR2lgOaVh1qxZo8\/8aieq+lAV4VUZl20w+YIgBa+p0Ff7URfJRZ\/k+HIRvE22q21+5CMfET\/60Y+KtCHVJt2YL\/tefZlavXq1WLp0aZF6ofpQ\/cJgmt5hkuaSwW2MJpJA5wQoeDt3AQ0ggXoCZRFlEsltipypIkAnFqtEb10UTxXSVfuVH\/J1n7RVgazrj2zrpJNO0go+9EG2gZxcXf6zreBVWVUJEbVPpn2QKRVVI0H1V2zBa2KDbh9bwSvFX5Xwk\/2EeESONbbYglftl\/Tvzp07RznT5X6rY6CcLlPe9zWveY245ZZbxvxz+SWwqg3mG\/MpQQL+BCh4\/RmyBRKITqBKQFZFr1SBqIsyqiJDPmTVh7d6Pl10Te6LjkshKR\/K+LfyxKKqqKL6gJfnUc8tRbTan3KE2tR+2OWb0lA3YUwdBDrxV9eHJn+hbcmnygbTlIamCK8akVQFZpPoMklpqIp4m04UdBG8uotTHf9qmzrxLbmqtuvEujqW5b7lyLA6QRK55FVRZN2XC05OjX6b5Ql6ToCCt+cOZvf6Q6AuglSVAlAWKerDXfcAratiUBXBK1cGqNpP\/rsu6qqLOKPSgYzsSdHRJLrq7E9F8Or6UNU\/nb9iC96qqGMT+z4IXpv0B4wn3Yuj+m9V6Ug6H1a9jOle6PpzV2NPSKA9AhS87bHmmUggCAGTaG9VtE8X3WuaKV+OvDblaDZ92pXHT58+fRQdbhIaddFLU\/tTEby6VAgbf8UWvHW5qXVlvPogeE3z2+WFrBO86otBVTRd58OmNKSmF44gNxc2QgI9JkDB22Pnsmv9J9CUL1o3qUqKTN0kLEy2kTVVUxe8pvaHELwhcnj7LnibXojkVZliDm9Z8FZNQqxLaaDg7f99lz3MkwAFb55+o9UDIqDLaVW7LwWfKjRMP5Gfe+65oyhr1aQzXW6tOnEO9mFGulxQwHTyUtXnXvnv6I9cwMEkAtpkfwjBizZUgV2OuplWaSjnrJr6CwuMpB7hNRW8NlUadDmtpgtPNH09qEvzUX9Tx3wMwVuV0jCgWx27SgJRCVDwRsXLxknAn0BdFQW1dTVy2DQJSpdWoJuIhvbrJq2h+kFZcKtVE5omKekmrelKrZkI3ib7Qwle1zq8dWkZpv7CCmWpC96mEa++mDR9xkdbvnV4Qwle2U7VJDedYLVJaVD9KhnFWiSkyUf8nQT6SICCt49eZZ96R8AkP1KNQDUJCd2M8ypo5Whu1fKxpm2q7dXl3+oi1rp0gKacYV1krqkcVJNIsvWHGhmuS2mo8kFdWTI5sa\/MwacOrxqBNo08mjBB\/1TB25R\/rb5w4f9dqjQ0+bJpImfdtRQyh7f89UAdC6ZR897d+NghEghIgII3IEw2RQIxCVSJg6ac0KaFDMoPWojBz372s+Lqq68uCvCXBYMuwqmb7FPer+6hXRYVtiuhldMMquyvEm\/lKLppsX+dGGoqFefrr6oIr2kfmsqSxZy0Vha88nrRiWUdpy4Eb\/m6g12XXHKJeN\/73lcsLCHHvm+Et4oFJ6vFvKuy7SERoOAdkrfZ18EQMK3JOhgg7CgJkAAJkMCgCVDwDtr97HxfCVDw9tWz7BcJkAAJkIALAQpeF2o8hgQSJ0DBm7iDaB4JkAAJkECrBCh4W8XNk5FAOwQoeNvhzLOQAAmQAAnkQYCCNw8\/0UoSIAESIAESIAESIAFHAhS8juB4GAmQAAmQAAmQAAmQQB4EKHjz8BOtJAESIAESIAESIAEScCRAwesIjoeRAAmQAAmQAAmQAAnkQYCCNw8\/0UoSIAESIAESIAESIAFHAhS8juB4GAmQAAmQAAmQAAmQQB4EKHjz8BOtJAESIAESIAESIAEScCRAwesIjoeRAAmQAAmQAAmQAAnkQYCCNw8\/0UoSIAESIAESIAESIAFHAhS8juB4GAmQAAmQAAmQAAmQQB4EKHjz8BOtJAESIIHkCTz++ONi\/\/79YuLEieKoo45K3l4aSAIkMBwCFLzD8TV7SgIkQAJRCTz66KNi9+7dYurUqeK4446Lei42TgIkQAI2BCh4bWhxXxIgARIggUoCFLwcHCRAAqkSoOBN1TO0iwRIgAQyI0DBm5nDaC4JDIgABe+AnM2ukgAJkEBMAhS8MemybRIgAR8CFLw+9HgsCZAACZDAiICP4N25c6fAH24kQAKHE5g+fbrAH27uBCh43dnxSBIgARIgAYWAq+CF0L344ovFli1byJMESEBD4FWvepW45pprKHo9RgcFrwc8HkoCJEACJHCIgKvgvfXWW8XChQuLB\/ppp51GpCRAAgoBvAh++tOfFuvXrxdz584lG0cCFLyO4HgYCZAACZDAWAK+gpcPdI4oEjicgHwh5PXhNzooeP348WgSIAESIIFfE6Dg5VAggfAEKHjDMKXgDcORrZAACZDA4AmEFrx7ekL0pJ70g93ohgAFbxjuFLxhOLIVEiABEhg8gdCCd2tPiP6mEOIYpS8HDhwQy5cvF5s2bdL2cPbs2WLt2rVi8uTJzgS2bdsmZs2aJcaPHy9uuOEGsXnzZrFq1ari79zyIkDBG8ZfFLxhOLIVEiABEhg8AQpe\/RCoErxnn322mD9\/fvBxUxa4FLzBEbfaIAVvGNwUvGE4shUSIAESGDwBCl4K3sFfBBEAUPCGgUrBG4YjWyEBEiCBwROg4A0veO+8806xaNEisWvXrlHj559\/vli2bFnx96uuukpMmDBB3HTTTQJpDO9+97vFF7\/4xeK3adOmiXXr1okf\/OAHRUrDxIkTi9JW2ObNm8cUh0yuWAreMI6i4A3Dka2QAAmQwOAJUPCGFbwPPPCAWLx4cSFuZf3VsviB4EUuMITtzJkzCwN0KQ0rVqwQK1euLFIoZA4xBLEUzoMfvAkDoOAN4xwK3jAc2QoJkAAJDJ4ABa+d4K2atCaFqa61sgiG4EX0V52QphO8GzZsGDMRjnm9+VyuFLxhfEXBG4YjWyEBEiCBwROg4LUTvKaT1nRVHeQiBBC82NRIrcmkNQrefC5XCt4wvqLgDcORrZAACZDA4AlQ8IYVvDKai9xcmberi\/BS8Pb70qPgDeNfCt4wHNkKCZAACXgTePzxx8X+\/fuLyUVHHXWUd3ttN0DBG1bwQugggqvW5JWT2FavXl3k9TLC2\/Yob\/98FLxhmFPwhuHIVkiABEjAm4CrYPQ+caAGXO2veqD3feGJppQGcFm6dOloQpoa8a1LaSgLZV36AlMaAg36Fpqh4A0DmYI3DEe2QgIkQALeBFwFo\/eJAzXgav9QBW\/VpDVZTgxVFyBMUWFBbl\/60pfEN7\/5zaLkGPJ2dRHesjD+2c9+dthKaxS8gQZ9C81Q8IaBTMEbhiNbIQESIAFvAq6C0fvEgRpwtX9ogjcQbjYzEAIUvGEcTcEbhiNbIQESIAFvAq6C0fvEgRpwtb\/qgb4\/kF1dNzOxawN4\/qwJUPCGcR8FbxiObIUESIAEvAm4CkbvEwdqwNV+PtADOYDN9JIAr48wbqXgDcORrZAACZCANwFXweh94kANuNrPB3ogB7CZXhLg9RHGrRS8YTiyFRIgARLwJuAqGL1PHKgBV\/v5QA\/kADbTSwK8PsK4lYI3DEe2QgIkQALeBFwFo\/eJAzXgYj9qD3\/3u98V73znO4UstRXIHDZDAr0gQMEbxo0UvGE4shUSIAES8CbgIhi9TxqwAVP7IXKxLxbZwH+xkthFF100SMFbLjumliSrcw0WoLjkkkvE1VdfLVC+rLyVV2TD7+A8a9YsMX78+DG765YuLrcnV3qzHS6hy5+hX7\/61a\/EGWecYWtKtvtT8IZxHQVvGI5shQRIgAS8CZgKRu8TRWqgyX4IXYiVhx56qLDguOOOK1aVu\/3228XChQsHJXilyLz77ru1K6ktWbJEzJ8\/v9JTTYK3fKCN8NSt8BZpyFg1qxPxVg1kujMFbxjHUfCG4chWSIAESMCbQJNg9D5B5AZ09svlkvEb\/mDJ5Gc84xljlk8e4gO9ToCWlw\/WuY2Cd27k0ZxO80O8PmLQp+CNQZVtkgAJkIADgb4I3ilTphS9RxQT0VyIXBnNxX\/L29Ae6CaRSnUFNYjbK6+8UsyZM0d88pOfFPPmzRPvete7xKWXXirOPfdcccUVVxRI1XQI9RxYaU2u1maSMlEV4YVNEyZMEDfddFORHoGc69mzZ4vly5cLddU4\/NvatWvF5MmTi5XiNm\/eLFatWjVKpUA7a9asOcxmOS7KaR4rV64Ub3jDG8TixYuL82JTUyzq2ivb\/JGPfERs3LjxsKh6XXqIw6Uc9JChXR9B4SmNUfDGIst2SYAESMCSQF8Er+y2jOY+85nPrCUxtAe6SXRWFZ179uwRixYtKoQulhPGJqPAEMFSTEIobtiwoRBz2CAQsf\/cuXO1wrPKKXWCF8J23bp1o7zh8tLGUmjjnDi3KnhxPohjuSwy\/o5zLV26dNQm9r\/uuutGf1ej3cg\/Vvsk00LQjspAPR72qTbrXjZs0j0sL+kguw\/t+ggCTdMIBW8ssmyXBEiABCwJ5Ch41Qlo+\/btK3p80kknFSkLumiuDknVA\/21r7UEmOju69YJMWPGIeNkxPbaa68toqBVTCDWIF6l4F29enUhXqXgLUclpQA8++yzRxHR0IJ3165dY6K1OttVEayKyZ07dxaR6nK\/5f4XXHBBIYhhvy5\/uSxWdS8OKgO0gbbLNqv2lfdPcQhR8IbxCgVvGI5shQRIgAS8CeQkeGVu7t69e4t+y7QFpDBMnTrVWOzi2KoH+rhx3kiTaOCuuw4XvE2f0MsR3vL+VaJZirnzzjsvSoQXQGWUWYVbTkOQKQeq4EU6AiYn6jbsX7a5vF9Z8FZFonHOHTt2FHaWI9ByvKkvE02+6HoQUfCG8QAFbxiObIUESIAEvAmkLnjL5cTKE9Bc7R9ahNclhzdVwSvzZ9W83aoIL3JnZcqFLrLdxMVU8Krn1wleGdV929veJpDfXM4x9r6QAzdAwRsGKAVvGI5shQRIgAS8CbgKRu8TNzQga+aqE9BQyxXVFtTN1f4hPtBtqjToPt03fc6Xk7xCpzSoEV6dQJViUubpliO8ar5uedg1pRe4pjTootKw67bbbitqQVelUMS+rkzbH+L1YcrGZj8KXhta3JcESCBpAvIzO\/JHEX3MbXMVjDH6KVlC5OL\/TSagudo\/xAe6TR3eKsGLiWymk9ZsauvWTVrTCd4FCxaMcm5lxFeX0oBjkaOLTU4ykyJWtqFOvEMUWBXBqEih5viaTlrTCV45GQ6\/qZPwYlxLvm0O8frwZaY7noI3BlW2SQIk0AkBV8HVibGak3Ztvy5loa6cWLkLrvYP+YEu+y5Z6sqG1UV41bJkalpBORoq\/y7LicnJb7qxbyp4cawUjpgYhg1Cd8aMGaPUhRtvvHFMyoBuVTeUHVMnqenKksnf5W+oWCFFc1NZMp3g1YnlVO4DZTuGfH2E9AkFb0iabIsESKBTAq6Cq1OjlZN3Zb9uAhpKiZVTFpo4udrPB3oT2Xx\/T7XkV1P6RErEeX2E8QYFbxiObIUESCABAq6CKwHTCxPatL9qAhpyc03LiTHCm8rISdcO3aSxFKw1KQ2Xgp2wgYI3jCcoeMNwZCskQAIJEGhTMMbobhv2Q+j+6le\/GrMCmm4Cmkv\/XO3nA92FdtrHpJwjK1MgsFJcXWpHKoR5fYTxBAVvGI5shQRIIAECroIrAdOjRnirJqCFntznyp8P9FRGIO1IkQCvjzBeoeANw5GtkAAJJEDAVXAlYHpwwes7Ac2FiSt\/PtBdaPOYoRDg9RHG0xS8YTiyFRIggQQIuAquBEwPJnh1E9Aw+QyT0GJvrvz5QI\/tGbafMwFeH2G8R8EbhiNbIQESSICAq+BKwHQvwdtFNFfHzJU\/H+ipjEDakSIBXh9hvELBG4YjWyEBEkiAgKvgSsB0J8FbnoCGxSGQl2tbTixU\/135D\/mBXq45q9bSDeWXkO2gni8mPZ5xxhmVzao1f+VOciGK8kGYQHbOOeeMJo+V6\/pi\/6bJZWWG6jnKNX59WJQrTvz0pz8tvpxggQxdrWSfc6nHDvn6CMUQ7VDwhqTJtkiABDol4Cq4OjVaObmJ\/TJlAfvij1wBLfQENBcmJvbr2m39gb5jhxBf+YoQ\/\/f\/CoH\/x\/YXfyHERz\/q0m2nY6QoPP3000cLKKAhiLfrrrsuydW\/dEsJlzsvBevq1atHIrZqkQe0d9lll4mPfvSjhWiEoNy0adOYvsv2lixZMmZxCvW8VbV+dbY4OUtzkM3Kdb7nbP368DU40eMpeBN1DM0iARKwJ+AquOzPFOeIKvtlygKEA5b6hci1WQEtjrWHt+rKv9UHOgTua197SOiq3ZgxA+vMCvF\/\/k90ZBB3WJ1MrhamnrDut+iG1ZzARPDWic9LLrlEXH311WLmzJnFWeD3m2++WSxbtqz4\/6VLl2qFft1v8iVh8+bNlSyxD84RcqPgDUmznbYoeNvhzLOQAAm0QMBVcLVgmtEpyvZXlRNrYwKakcGlnVz5tyZ468Su7AtE7113uXTf+BgT4Sgbq1tWWIpHCOQJEyaIm266ScilgyEky\/+GmrPlpYzVT\/5SrL7kJS8RV1xxRWGCTLHA\/y9evLhoH1tVigLa2LBhg1i7dm0Rta3bsC82LBvss0BF3Wpu5Xbr+i\/F84oVK0Zmq3xkW0jBWLhw4WgfpFycdNJJoizo65Y8rmNd5tba9WE8gvPckYI3T7\/RahIgAQ0BV8GVCkxpP3JwZVQ31Wiujpkr\/9Ye6H\/7t0IsWtTsbkR5keIQabPJ9zQVvOVUAF16QFmMSuG9YMGCQnTKXFgp8mQqwrRp04oIqYlQL+fvVuXQou1PfOIT4h3veIeYPn26WL58uTj77LMr0xaahLMuwlseV+V0kXK6RDlqW2aviue6fWV\/YLOM4JfP3cRa7W9r10ek8Z5KsxS8qXiCdpAACXgTcBVc3icO0EC5nJhMWehqAppLl1z5t\/ZARyoD8nabNqQ0\/Pu\/N+3l\/LvNsramgrecHlFOi5Di9W1ve9uY1cVU4XbjjTceFp1Vo6doA1FeiN+mFcp0E8nUyWfo11e+8hXx4Q9\/uOAIwVu2zRRw3aQ1eU7Z\/7KoVvu3cePG2ui0qeCF3eVob\/n8ukh4VaS6tevDFHim+1HwZuo4mk0CJHA4AVfB1RXL8gQ0\/P2JJ54QM2bMKHJ0c9tc+bf2QH\/uc\/W5u2XQkdMaYkR40QU1T7X8KV9XOUF2W6YtQPCWI6WugldFKs993333jXJ01fzdKjFqOv7rUhpkG1XRafXlA\/uqaRvl6LSp4N2zZ0+RolFO64CdO3bsKPyks5mC19TjbvtR8Lpx41EkQAIJEnAVXG12BaJWzc1VUxZgx+7du8XUqVMpeGM4xVTwRo7wNqUGqCIM4qkcLaz71C6xlQWvScWCJhFmE+Etu6+cPlEuR1aXw1sVnZbn8BG8VZPP1KixFL6+glc9von1+PHjRwhbeyGMcc0l1CYFb0LOoCkkQAJ+BFIWvCYT0FK238Qzrva39kA3zeH92MeilyhrmmglUxR27tx5mOAtVy3QicWqCK\/M19X5s0mENQnesqhVz6FGcd\/whjeMKUeG\/eqi3j5VGqQNJikNqsjUienPfOYzxT\/LqhJqFFe1H\/uYpDTURdMpeE3uOHb7UPDa8eLeJNBrAlKUpVDT1QW0q+ByOZfJMbYroKVmv0kf1X1c7W9N8MLYpihv5HQGycu0Dq\/cDzmzEFpSuG3dunWUHmAieHFeXY1fNdcXOax1IgxtNE0uk5HR8mIRaiQVUWuZv6sKu9B1eMvjt2nSWjmvtiySVc5lga7+3XTSGgWv7R3Gb38KXj9+PJoEekXAVbCkAiEV+8sT0JC2gFJiTRPQUrHf1Z+u9rcqeDFpDZPXdBvELiar4b8tbWrpKpxSt9KauvoYKiZ86EMfEtdff\/2opq2p4JWiVy27NW\/evDGVBJpEmBS06nFlVLrV0srnwTGoDFHeymXD8LvJSmtVdXib2tfl6a5Zs2Z0mFp+TeUsxTCqY6CNl73sZU5lydQ6zMzhjXvRUfDG5cvWSSArAq6CJZVOdmm\/LpoLgYsIlukEtC7tD+FDV\/tbFbzoaAIrrYXgzTaGQaD166OnWCl4e+pYdosEXAi4ChaXc8U4pgv7cc79+\/ePWQENIrcpmqvrfxf2h\/SDq\/18oIf0AtvqGwFeH2E8SsEbhiNbIYFeEHAVLKl0vi37qyag+eY+t2V\/LH+52s8HeiyPsN0+EOD1EcaLFLxhOLIVEugFAVfBkkrnY9pvOwHNhUlM+13ssT3G1X4+0G1Jc\/8hEeD1EcbbFLxhOLIVEugFAVfBkkrnY9jvOgHNhUkM+13scD3G1X75QH\/\/+98vXvWqV7menseRQC8J3HPPPeLiiy9unLzXy84H7BQFb0CYbIoEcifgKlhS6Xco+0NMQHNhEsp+l3OHOMbVftSaxQN9y5YtIcxgGyTQOwJ4EbzmmmsESp5xcyNAwevGjUeRQC8JuAqWVGD42g+h+6tf\/SrIBDQXJr72u5wz5DE+9kP04k9XG2zfu3evOPnkkwXKyLW5yXPjnCeeeOKoqkeXNun6n5o9bfqo63NB6FLs+nmBgtePH48mgV4R8BEsKYBwsV+mLOBY\/IHYQYUF3wloLjxc7Hc5T6xjcra\/a9sxDu+7775iDKKMHZaX7tqm8jhJzZ5Y45jt9pMABW8\/\/cpekYATgdwfaKb2y5QFFI9\/6KGHCpELkQGRa1oz1wlww0Gm9sc4d4g2c7Y\/FdulHXJMYnxC\/HY5LuXYSIVRiLHKNoZHgIJ3eD5nj0mgkkDuD7Qm+6vKiWEVtBS2JvtTsLHOhpztT8l2daLkww8\/LGbMmOFU1zn0eEmJUei+sb3+E6Dg7b+P2UMSMCaQ+wNNZ38b5cSMATPCGwpV8HZSHPvIKb7tttvElClTxKmnnlosT93lliKjLnnw3HkRoODNy1+0lgSiEsj9gabaj0\/CcgIaoMmUBZcV0KJCVxrvE\/8UPsHb+C1F9tImrNyH9BuMaaQ3tD2pTnJMkZGNj7nvsAlQ8A7b\/+w9CYwhkPsDTc1\/RGS3ywloLkOrL\/xTyTm18UGK7FWb0BdMasO4RiWHLqK9KTKy8TH3HTYBCt5h+5+9J4HsBa86AW3Pnj3i4MGDYtq0aZ1PQHMZWrkLipztT9F2nU34aoFUB7zMIdWhzUh6ioxcrjMeM0wCFLzD9Dt7TQJaAjk90HQT0CAC0IccI4xwSE78dQMoZ\/tTtL3KJrWEWZvR3hQZ8VZOAqYEKHhNSXE\/EhgAgdQfaE0T0FK3v2kI0f4mQvF+T5F9k01tR3ub7InnHbZMAv4EKHj9GbIFEugNgVQfaGqZJsBGJBc5jOUJaKnabzpAaL8pqfD7pcjexKY2o70m9oT3DFskgTAEKHjDcGQrJFAQkMKsi1W6QrggpQeaLpoLgYsZ61V5iynZ7+IP2u9CLcwxKbK3samNaK+NPWG8wlZIIBwBCt5wLNkSCTAHM8AYwEN1\/\/79Y1ZAg8g1KSeW+wOZ9gcYQI5NpMje1iY12ovrBZPaQm629oQ8N9siAV8CFLy+BHk8CSgEcn8gdGV\/1QpotpHyruwPdRHQ\/lAk7dtJkb2rTViOGCXMQldycLXH3hs8ggTCE6DgDc+ULQ6YQO4PhDbtb5qA5jKM2rTfxb6mY2h\/E6F4v6fI3scmNdqLFCBULvHdfOzxPTePJwFfAhS8vgR5PAkwwms1BkwnoFk1+uudc38g0\/6nHalG\/MeNGyeOPPLIIIKtbkylyD6ETbINuQiLz4IVIexxua55DAmEIEDBG4Ii2yABCq7aMeAyAc1lUOX+QB6y\/VUR\/6OPPrpbpT2sAAAgAElEQVRYIjqEYBui4FVfILBgBaK9yO11WZ449\/Hpck\/hMf0hQMHbH1+yJwkQyP2BENp+iBiIFeQU4gGLh63pBDQXd4a238UGn2OGaL8u4o8JV2okUt3HR7ANVfDKfsvxhb+7LFiR+\/j0uTZ5bP4EKHjz9yF7kBCB3B8IIewPNQHNxa0h7Hc5b6hjhmK\/a8TfV7ANXfDK\/qslzJDbaxrtzX18hrpO2U6eBCh48\/QbrU6UQO4PBFf7pYA5cODAmGguqixU1cyN4UJX+2PY4tJmn+3HGFFfhnwi\/q6CjYL3EAGMNVRygE9Mo725j0+Xa5LH9IcABW9\/fMmeJEAg9weCrf0mn6PbdIut\/W3aZnKuPtofK+LvItgoeA8nYLNgRe7j0+Qa5D79JUDB21\/fsmcdEMj9gWBif4xyYqFcZWJ\/qHPFaKcv9ssFD7CACPoko7kxIv42go2CV0\/AdHni3MdnjGuWbeZDgII3H1\/R0gwI5P5AqLNfnYAGVyBVAQLGZAW0tlzXZ\/5tMfQ5j5pji3YgdDH5LPYYMRVsFLz13m16ecj9+vIZ2zw2fwIUvPn7kD1IiEDuD4Sy\/fJzNP5dRuogXmxXQGvLRX3j3xY3n\/OoY2Tfvn1FTuhpp53WyRhpEmwUvM2ernt5yP36au499+gzAQrePnuXfWudQO4PBGm\/jMip5cRifI4O7aC+8MfM+brJflsVcD+pgfjC0m9zAgGvmqSIiC5qvTbZH8gMbTOu0d4Ux06XNuleHrq0J+aYYdvDIEDBOww\/D6KXz3ve88T27ds77WvODwQIhfvvv1\/88pe\/FJMmTYpe6D+Go3LmDx7S\/t3Tpoljjz22QFQnaEMwlKLYRAw3TVJMib9ttDcl26Vfu7ZJfXmQX3Z2797d6QtNiDHPNoZJgIJ3mH7vZa8peO3dWp6Ahr8\/+OCD4owzzihKFeW2dS0QXHkhYgth+9hjjxURUrDHCmNdbaoItpmkmBr\/smCTk+l0XFOzXX0B6jJiDjvwpQclzBDBB9Ou7enquuB58yZAwZu3\/2i9QoCC13w46CJ1mFyEB1rOEZwURYvOKzIloRy9TUXwwuYnnnhCPPLII0XU+cwnnxQvffLJxkmKqfJXBRtEry5dJEXbU7JJvjyAJaK9EL3cSCAnAhS8OXmLttYSoOCtHyAmE9BSesC6DPeU7a8SuWo\/uxa8ELlgePDgwULwHnnkkUVqBQQi\/h8bor9V6Q8p81ejvehPWbClaHtqNsGeHTt2FONgwoQJhfBVl4B2uWZ5DAm0RYCCty3SPE90AhS8hyOuWt2qagJaag9Y20GTov1fs+hEF4IXwhbnBTv8F8IW6RQQhU1pFWXxmyL\/Mn5pI75mqIItRdtTs0nagyi5HKsYJ\/i76fLEFpcDdyWBoAQoeIPiZGNdEqDgPUS\/anWrpmhMag9Y2\/GUiv0yJ9fW\/jYFr4zmIm0BG4Tu8ccf77wU9J8pk+5Sz\/FUU3qkYMO\/pZbOk8p4luO4bI\/8O343XZ7Y9prg\/iQQigAFbyiSbKdzAkMXvDaTi6qcldoD1nZQdW2\/q9CV\/YwteHXRXKQsHHPMMY3RXBNfwP6T77tPvGHyZGfhbHKeUPuogg3RXuSnpiTWux7PVdHxMiO1IgZ+Y7Q31AhlOyEJUPCGpMm2OiUwVMFbNQHNZXWr1B6wtgOqK\/tt0hbq+hRD8EqRK9MWZMqCTFuwZWxq\/4uPProy1zfkOUO0BcGGknyIdqNCicu1E8IOU4EZ41wmbdZdX\/gNlRxwP2K014Qm92mbAAVv28R5vmgEhiR4ddFcPKTHjx\/vFVnrSjCGGhRt2+8b0S33O6TgrZqAhslGsTad\/Uh1yGFDObj\/\/M\/\/FJMnTxannnpqEpOx2h7PTX4ysce2\/nHTOfk7CYQiQMEbiiTb6ZzAEAQvHjj79+8vPr3isyHyDyFyQ0WkTB5onTu6xoC27A8tdGWXfAWvzwS0EH6tsz914SvHDq6nAwcOFNdXVQmzEKxM2mhrPJvYgn1M7XFd7c7UDu5HAi4EKHhdqPGYJAn0VfBWTUBDpYXQuXKmD7QkB4DFA9nH\/lDpCzobXAWvWjMX7cp0hbrliX0YVB3bZH9dSbMY9ti0qY59XFf4PI9\/6\/LzfGrXo609jPbajEDuG5sABW9swmy\/NQJ9ErwhJqC5gLd9oLmcI+YxMe2PFdVVeTQJRnVfmbKAY2Q5sXLN3JisfQR7isJXN3a6Fmwxx7PL2HCxh9FeF9I8JgYBCt4YVNlmJwT6IHh1E9DaLO7u8kDrxNkVJ41lf8yoro3grZqAZlIztw0\/2Qj21ERv1djpUrDFGs+uY8HHnq5fHlz7zOP6Q4CCtz++HHxPchW8sSaguQwInweay\/lCHxPa\/jaiuiaCV1czF9HcmBPQXHxjI3hl+6nk9jaNnS4EW5NNLj7yOcbXHvXlAS\/yyJHmRgJtEaDgbYs0zxOdQG6CFzd\/PETlBDTkDSIvN9QENBfgvg80l3OGPCak\/W1FdasE7xFHHOG8AlpIpjZtuQjeVISvydhpW7CZ2GTjH999Q9mDex5ypFOYGOjLhMfnQ4CCNx9f0dIGAjkIXpmygAcH\/sjlTWNMQHMZMKEeaC7nDnFMCPvbjurqBC8mnanL\/MaomRuCd7kNH8GLtrqM9tqMnbYEm41NMfxZbjOkPfLlAek4TStAttE3nqP\/BCh4++\/jwfQwVcErUxZQ6kgtJwaR2\/Ys+qbBEPKB1nSuGL\/72t9FVBccZMoCxgf+H+Oi6wloLv7xFbxdil7bsaNGe+EvrDAWerO1KfT5Ywre2LayfRIoE6Dg5ZjoDYHUBC+it7JmLh6OMpqbcjQjtQes7eD0sb9tsaurmYv+og\/IbURUN7cthOBFn7uY0OY6duRxMa5vV5tijZvU7InVT7bbTwIUvP306yB7lYLgVT91SpGL6E+K0VzdIMn9geZqf5tiVzcB7fjjjy+iuqEEY1c3gJD2ty16XccOWKvVVeBHvLCEqJHtY1OMMZCaPTH6yDb7S4CCt7++derZAw88IBYvXiy2bds25vjZs2eLtWvXFstuYsPn+eXLl4tNmzYVfy\/\/Lg++6qqrxJo1a4q\/Tps2Taxbt07MnDlzTNu33nqrWLhw4ejfVq5cKebPn29tf5eCV05A27NnjwDD008\/vShY3+UENGuALS3c4GKX6TEuD+Q2xG5VzdxjjjlmTCQ3pGA0ZRZyv9D2tyl6XcZOmZ1sA\/8eYsGKEDaF9G9q9oTsG9vqPwEK3v772KqHd955p1i0aJFYvXq1mDt3rvZYKXYhYJctW1bsA2EL4aqKYvzbrl27xKpVq4rlb2+44QZx3XXXjRG9OGbp0qWjf5PnX7JkibXobVvw6iagIaqDKO+zn\/3s5PJzTQZC7g80W\/tjit2qmrl1E9BCC0YTn4fcJ5b9bUxmsx07ddzUEmbI7XWN9oa0KYSfU7MnRJ\/YxnAIUPAOx9dGPYUAhVBVhWv5QJ1wlZHhBQsWFEJVRm3Xr18\/Es5loSyPgbCWwhnnQvsbNmyotUHXmTYEb9MEtNwfCEOyP5bYldHcgwcPFhPQjjzyyGICmknN3FiC0ejiD7BTTPtji97QYx\/tofQW7hmu0d7QNvm6ODV7fPvD44dFgIJ3WP5u7C3E5ubNm0dRWd0B5cit3Ef9940bN2pFqypm8flfF002iTK3LXhlNBfR27oJaLk\/EIZif2ixq5uAJiO5NpPPYgrGxos\/wA6x7Y8pemONfZ8FK2LZ5Orq1Oxx7QePGyYBCt5h+l3baxmBvfvuu4vIBNIRsM2bN28kgHXpDKrglWkN119\/\/Zh0BrmPGh2G4FXTGeQ+rmkNoSO8uhXQmiag5f5AGIL9IcUuhO4jjzxSVFbAJkWua7m52IIx9u0utv0xc3pjjn3X5Ylj2uQyFlKzx6UPPGa4BCh4h+v7w3ouUwww4Urm3aoiGGkOyMXFZDU1f1cVszIVoUrwqjm7VYK3nB5h6qJQgledcY1zI\/8OpcRMJqDl\/kDou\/0hxG7VBDSIXKQv+GyxBaOPbSbHtmF\/LNHbxti3jfa2YZOJX+U+qdljYzv3JQEKXo6BRgJqxPXcc8\/tpeDVRXMhcCHwbaJ1uT8Q+my\/j9gNlbLQdLG1IRibbPD5vS37Y4jetsa+TbS3LZtMfZ6aPaZ2cz8SAAEKXo6DRgLq5LILLrigUvCqlRpySWnADVwuDoFILsQtRK5JNFcHLvcHQl\/tdxW7upq5phPQGi8szQ5tCUYX20yOadP+0KK37bFvEu1t26YmH6dmT5O9\/J0EVAIUvBwPjQTK1RRyn7RWNQENi0O4lg\/qyye\/3B9oOvu3CiF+0jjKD+2gi+ZC5JZr5lo0abxrm4LR2CiLHdu2P+Qkti7GvhrtxUs2FqxQty5sqnN3avZYDE3uSgKM8HIMHCJQNVmsXDVBVzZMV5asPCGtqiyZLGUmLWkqS1aVq1uXw+syAc1lbOT+QOib\/TZiF2IN\/ccf5OJiAlpdzVyX8dF0TNuCscke29+7sD+U6O1y7KsrNEL0yjSqLm3S+T41e2zHJ\/cfNgFGeIft\/8N6j+gtVk+TK6LpJrLp\/q288ER5shtWaNPV78W\/rVixQsh6vaYVGiBusW3fvn3UB53g9ZmA5jI0cn8g9M3+plSGqpq5ISaguYyfLgSji51Vx3RlfwjR2\/XYV6O9GH9YsKJrm8p+Ts2ekGOXbfWfAAVv\/31s3UMpQuWB559\/\/piFIfDvqSwtrIpc+f+hJqBZgxvo0rwunGIdoz6Q\/\/G447SnaWsCmksfuxKMLrbqjunK\/hD5vKmIOWmHnFOA6C\/Er83k2VD+pOCNRZLtdkGAgrcL6jxnUAIy2otGt2zZUiztG2ICmouRqTw0XWzHMX2xf\/e0aWL7sceOwaCbgHb88ccnISSkoV0JRtfxUj6uS\/t9o7wpjX35Zer+++8XDz74oDjzzDOdJ9KG8m0f7g8hWbCt\/AhQ8ObnM1pcQQDC9+abby4eDCEmoLmATumhOVT7b3zgAXH\/lClF\/m1Vzdw2JqC58O9SMLrYm5LghS0+ojfFa3fv3r3itttuKyaznXrqqUU98C63FBl1yYPnzosABW9e\/qK1NQRCLTzhAzn3B0If7P\/8\/v1jJv3ICWj4JGyzzK\/POHA9loLXldzTx\/mkNqQ49qVNKJWINDJ8uUJ6g281GVfKKTJy7QuPGx4BCt7h+by3Pabg9Xdtzg80fAb+t337xPcffriopQyhG7Nmrj\/tw1ug4PWn6hrlTXHsqzaBDJZ8xzg\/8cQTO4n2psjIf8SwhaEQoOAdiqcH0E8KXn8n5\/ZAK09Q\/MbxxxeRsFNOOUUgPze3jYI3jMdcRG+KY19nk8mCFWEoHt5Kioxi9ZXt9o8ABW\/\/fDrYHlHw+rs+lwcahC4e\/HKCIj7x\/r\/PfGYR1UXeIyJgqacv6LxFwes\/hl1TG1Ic+1U22SxPHIbo062kyChk\/9hWvwlQ8Pbbv4PqHQWvv7tTfqDVrZC37aijitXUKBj9x4BPCynxt43ypjj2m2xqO9rbZI\/P2OGxJBCbAAVvbMJsvzUCFLz+qFN7oJmukCcXmEhJcLl4g\/bbUfvzP6\/ZH29AWGovwPbUUwEacWjC5HpsM9prYo9DN3kICbRCgIK3Fcw8SRsEKHj9KafyQNOtkIdyc7qyTOrywRSM\/mPAp4XQ\/GsFrYmhTUvtmbRR2qdN8WtzPbYR7bWxxwEtDyGBqAQoeKPiZeNtEqDg9afd5QPNdYU8VdOEFlz+RO1aGKr93sJWwfx3f3foLzZlytSxP368fpU+1ZttCF\/b61GN9uIFEfV7Q2629oQ8N9siAV8CFLy+BHl8MgQoeP1d0fYDDQ9oNTfXdoU8NbqL3g9VMPp7PkwLJvx9xK0qZk0tNs3lbRr748bpzxhT+DbZVMUAkzlRwgzXE0RvqGWJXe0x9RX3I4GYBCh4Y9Jl260SoOD1x93WA61uAppNUf3yF2sTweVPKV4LfbPfRdxC1Kor5D355JPiiCOOKCpvuG4motdm7KviN0XBC05qtBeCFwtW+G42jHzPxeNJIDQBCt7QRNleZwQoeP3Rx3ygmU5AM+1FObrLCK8puTj7uYpbaQ1ELgQ\/\/mAcyhXy8N+HH37YayGR0IIXNkvRm6rglVzlNY0Xyao8eNMREfP+YGoD9yMBVwIUvK7keFxyBCh4\/V0S44Gm1syFhXjwYvIZHr4+m24+Ut8ipD58Yh9rI3DrUhFkNPeRRx4pTNatkKfug\/rKEydOLPaz2ZpEr8vYjx3pdbFJx0SdBIpoL9IcbL6klMUzosWh0iRsfMh9ScCHAAWvDz0emxQBCl5\/d4R+wKI9\/JHRJQgVlwdtuWe66C72oeD1HwNVLdgI3C9\/+bHahT9kNBdjAz6T0VyIqLoFQ6R\/YSNW0pswYYJxh2MK3lhR3lDXY1mw4u8uyxOHtsfYedyRBAIQoOANAJFNpEGAgtffDz4PNJmygKV95QpoEDAQuaGjQVXVpih4\/ceAbMFU4KrR2yb+ELqI5GKcYYO4xdiwHR9IcUA7EMonnHCCcbS3TvS6jP1cIrzlUaGWMEO01vQl1IVRuBHJlkjAjwAFrx+\/KEdDLNx+++3i+9\/\/vrjrrrvE9u3bBcTcc5\/7XPHyl79czJ49W4wfPz7KuXNulILX33suDzSbmrn+Fj69lgDWFNBtTYIrxPljttGl\/S4Ct8xCZ786AU1Gc4899thC5NqmJajnQ1v79+8vJriZRnspeA8RxLWOSg64fk2jvS73h5jXC9smARsCFLw2tCLve\/fdd4vPfe5z4p\/\/+Z+Lz3zYEB07\/fTTBX7DzV1GRRYsWCDOO+88MX369MhW5dM8Ba+\/r0wfaKEnoNlYXreWQJeC0aYPVfu2aX8IgVsleHHfwlaegNaUsuDCUI324rx1KRF1dXlNx75qY64RXrUPNgtWuDBy8SmPIYEYBCh4Y1C1bBM37HXr1ok1a9aIV77yleKtb32rmDNnjjjppJOKcjxyQ3meX\/ziF2Lr1q1i\/fr1RRQYwnfJkiVi0qRJlmft3+4UvP4+bXqgqRPQ8BkUfyAyfCeg2VhOwWtD69C+MQRuleCV\/66bgOZmff1RiPIiIACB3RTtrYryNo39sgWxxS7OZ2uTK1vT5Ynbsse1HzyOBOoIUPB2PD5wk\/7ABz5QRHHf\/e53i2c961nGFu3cuVNcf\/314p577hGf+tSnCuEx5I2C19\/7ugdaqJq5\/tbVpzOg\/TYjpCH6UyUY8Ym5LlJpem4TkeuymIN6fnUCmqy0IPO2Q\/TBtK\/YzyTaS8FbTbQp2kvBazMauW9qBCh4O\/bIwYMHi8jEySef7GzJ\/fffX4hd5MUNeaPg9fe+fKDJJUkxNmWVhVgT0GysrovuUvA+TbINkYvz6MqJQeBivIQS7DZjQ+7bFO2tSmuwEXNtRHfRHxubXFjpjqmL9nZhT6h+sR0SoODlGOgNAQpef1diwuTu3btHs7ZD1cz1t+zpFih4DyfZlsCVIlfm5eomoKUUYa+L9uqivKZiri2x25XglSNMF+01ZRTqemc7JBCSAAVvSJqB2kKu7r59+8RTFcUdx40bV5TiUfN7A50662YoeN3cp05Aw7h74IEHxKxZs4q8XNtyUW4WmB1VV51BtpCS4DLr1di9TO1vU+TCQilyIXhkzVxZUkztgan9LmxcjlGjvfILBdpxFbyq2EU7servyr52LTDVaC\/uB\/iSiBdiLjzhMhp5TNcEKHi79kDp\/HfccYf44Ac\/KP7jP\/6j0jKUJVu7dq2YPHlyYtZ3aw4Frzn\/qpq5aAFR3hQfaE3RXSnM9u7d2+kndXMvHL5nnWBsErm+ubhla2TKAtKu8P8mE9BSE7yqcER6DvoA0fbio48Wc0odbhKXbYtdmNdkk89YszkW9wSUMMMXH9w7Urw\/2PSH+w6TAAVvQn5Hwf7ly5cXVRguuugiccopp2itQ67ui1\/84sHn7JbhUPA2D+aqCWhY6jelB6yuJ0MUvF2I3HLKgozkmkxAS1XwYjyp0d6zxo0Tv1uaN1EnLrsQu6ldjzLae++99xZfGGfMmNF8w+EeJJAQAQrehJyBT8mLFy8uSo3Nnz8\/IcvyMIWCV+8nm5q5qUSUhix43\/nOo2svuNCRXCkG4XtZZQGRUJT3sk1pSVnwSqjSxj959NEibafpZa8rsZua4JX27Nixo0CJZZ1Vfnk8JWjlkAlQ8CbkfXw2QmT3937v98Sb3vSmhCzLwxQK3rF+Umvm4heZw1hXMzdVwWuSv4s+5iC4dFdT25FcKXJ1E9COOeYY55JoufBHtPfcBx8USH\/BdYGqJLheyvmpXYrdVAUvGIGX9LXkZ7o8cR5PE1rZRwIUvIl59dvf\/rb42te+Jq655hpx2mmnJWZd2uZQ8IrioS1LiclyYnKyickDiYK3vTHehciVQheR3KYJaC4kchG86BvKk73o0UcLkYsN14mav9612E1Z8MocXnm\/gK2myxO7jCseQwIhCFDwhqAYsA1EHj7\/+c8XC0ogR0pXiQHLCV9++eWjT3EBT591U0MVvFUT0GTxfxunpip4TfJ3c4jwNoncL3\/5sSLqGLqObdUENETnkL4QastJ8KLPsloDSnChnjleBM444wwxceIzxiCJXY2hin9q12OVPWoJM4hhk5frUGOO7ZCAKQEKXlNSLeyHMmRf+MIXiuguJog8\/\/nPF\/i8WN4oePXOGJrgldFcCCRseMj45tSl9oCVns5Z8DaJXDUnN6RgVFdAkzVzbSagudzyQtrvcn7bY9TyZLiOnvnME5MQutKI1K7HOnvwGyo54L7EaK\/tSOT+bRCg4G2DsuE5HnzwQfG+971PvOhFLxIXXnghqzAYcpO7DUHw2kxAs8RX7J7aAxY2mebvYt9UBJeNyFX9FMJ+3QpoqOyCSUaxtxD2x7ZRbV8VvOUUhrvv\/nmRq2o7cS+k\/aldjyb2NC1PHJIP2yIBGwIUvDa0Iu8rqzRg4tprXvOayGfrX\/N9Frx40CA3FzmGiOTiITx+\/PgiohtyM3mghTyfSVs5Cd46oWtSXcFVMHYRzdX5ztV+k3EQYx8peMti97HHHi+ilXKZZFnJIYYNdW2mdj2a2lO3PHHbDHk+EpAEKHgTGgu4SVx11VVi0qRJRaQXK6pxMyfQN8FbVTMXubmxcuRMH2jmXvHfM3XB6ytyVUK2ghFCV52Ahpx\/vAx1FZW0td9\/dPi18OelW+z27XeNWVSh62hlatejrT1d8\/MbHTy6bwQoeBPzKIp6I8L70pe+VLzlLW8pxG9549LCeqf1QfDGTlloGu62D7Sm9kL8bpq\/i3O1JbhCilxbwStTFtBXmZuLlIXQE9BcfNcWfxfb1GPKQhe\/HTjwdMWG8ipiXUYrU7seXezpkp\/vOOHx\/SJAwZuQP2VKw7Zt22qt4tLC\/RO8uglo+IwaOmWhabi7PNCa2vT9PSXBG0voSkZVglGmLMi6uaisEHsCmovfchC8OrGLKgxNY7+LaGWTTS4+8jnGx54u+Pn0lcf2jwAFb0I+xZr1P\/7xjwX+W7dxaeF+CN6qmrnIze3qk7TPAy3WpdS14I0tcusivFXlxNqYgObiz5QFr07o\/t1Th0qTmYx9NVqJl1FMaou5mdgU8\/zltn3taZtfm2x4rvQJUPCm76PDLMRNA7l6uhq9GXYnmMk5pDRU1cyNMQHNBazvA83lnE3HdCV4q4SuyeSzpj5V\/a6uXiWjuqlGc3V9SFXwlsUuhK7c5MQ1m7GPyaOY1IZc+piVHGxsch1zNseFsqctfjZ94779J0DBm5CP8ckHEd7f\/u3f1k5Ye\/LJJ8XNN98svv71r4urr75aTJ48OSHruzclZcFbNQGtq9nfVd4K9UALNRpsJqzhnL6Cq81obpkRBC6EAHyAl1mZstBVtN\/Fh778Xc5Zd0xVVFc9xkXw4ng1WgkfIfc39Jba9RjSHskP7FK7D4b2I9tLgwAFbxp+KKxADu+73vUuMX\/+fPHWt751TAT3nnvuEZdddpn4zne+I17\/+teLa6+9VmC2PrdDBFITvIj+yHJi+K8sJ+ayAlpbfg75QAthc1uCt4toLviUJ6DhpRYpTRBPELy5bakIXhOhK9m6Cl55vLxmQiz8UvZ3atdjavbkdn3Q3m4JUPB2y3\/M2fHGu379erFy5UqxYsUKsXDhwuKBuHHjxiKii+2SSy4R5557Lhel0PgtJcGLqAUeDtjw\/xC5bU9AcxnaqT3QbNIZ0F8bwdVVNBfXNIQtWOOPmrKAPsRYWthlLLgcY8Pfpf2mY2yErmzrhUKIOZ6LrqiTTnG9I80hROnA1K7H1OxpGg\/8nQRUAhS8iY0HLC\/8rW99S1x++eXijW98o\/if\/\/kf8cMf\/lC87W1vK2rzxp4kkRgOK3O6FLzygbdv3z6xa9cuMW3aNHHCCScUQjfEg88KhMfOqT3QYgjerqO5iOBC9ELolldA61owegyd4tAu7XcRu7A5hOCV3OT1g7+HWF43tesxNXt8xyuPHxYBCt4E\/Q3Re+ONN4qLL764ELhYjOJlL3tZbxeiuPXWW4tottwQ4UZah+3WtuDV1cxFjWTkYc6YMaOzSgu23NT9U3ughRS8XQhd2xXQuhSMPuNGHtuF\/a5CN2SEt8xOLcGF9BTXl97UrsfU7AkxZtnGcAhQ8Hbsa3zaRFQQIre8bdmyRSxfvly8\/e1vL\/4gIoStTwtPQOwuXbpUrFu3TsycOVPceeedYtGiRWLJkiXWorctwaurmYt0BUy8yP2BkJr9voK3y7QFsMQqaNhw7R5\/\/PGNL0FdCMaQt8A27fcVujEFL9qG\/1HJAfcL12hvatdjavaEHLtsq\/8EKHg79rHpYhOqmX1ZeEL2fe7cuWLZsmWjLt5www1iw4YNYu3atVaVKGIKXl00FyK3XDM39wdCava7Ct6LLtLXR41ZUkwXzUXKwjHHHGM8Aa1NwUPb+HcAACAASURBVBjj1teG\/aGEbmzBK9v3WXAhtesxNXtijGG22V8CFLwd+9Z0sQnVzL4sPCGjuatXrxYQvXKr+vcmV8UQvLLKAtIUZJWFupq5uT8QUrPfVvB2kbYgVz9TJ6DJkmJNY7b8exuC0dYmm\/1j2h9a6LYleHEe1+V1U7seU7PHZmxyXxKg4OUY6IxAOZ2hLHht0xpCCd6qmrkmE9ByfyCkZr+p4G1b6FatgIYZ+jL1yOXCiikYXeyxPSaG\/TqhC7vUxSNs7VT3DzlprckO22hvatdjavY08ebvJKASoODteDwgcogqDK985SudSo3JCPELX\/jCLMpeqbirBK9MdViwYIFVHq+P4NWlLMhyYjaF\/3N\/IKRmf1Md3iqh++UvP2acRmB6C7CdgGbarrpfDMHoYofrMSHtjy1024zwqjxtor2pXY+p2eM6TnncMAlQ8Hbs9wMHDhQlyP7rv\/5LXHjhheLVr3610YxeHHfLLbeIz3\/+8+IFL3iBuPTSS4t80py2FASvbgIaJp+51szN\/YGQmv1VglcndJGfG1JwyWtJRnNtJ6C5XIsx7Hexw\/WYEPa3JXS7ErzyvCbR3tSux9TscR2nPG6YBCh4E\/C7XDIYwhfRzTe96U3it37rtwSitqqIffDBB8Vtt90m\/v3f\/13cdNNNxYQuCN1zzjlnzKpsCXTJyISuUhqkyMXNG3\/kCknlCWhGnSjtlPsDITX7y4K3SuhKN4QQXGgrxAQ0l\/ETyn6Xc4c4xsf+toVu14IX51ejvXjJLtdZT+16TM2eEGOWbQyHAAVvQr6WUVtUJ\/je975XadkrXvEKcd5554nXvOY1TmkQqXS5zUlrMmUBjNUJaKGX+c39gZCa\/VLwNgndUIIXQheR3BAT0FyuMx\/B6HK+0Me42N+V0JV9911aOARD3JNQwgwv3xC9Mo0qtesxNXtCsGcbwyFAwZuorxHN\/elPfyrw2Wvnzp3iOc95TlHL8YwzznD+3J5aV6tydZvKklXl6ur+vWoCGtIWYmy5PxBSs3\/cy3+9FJbirLrSYi6CS6Ys4Fj8kSug+U5AcxlfLva7nCfWMab2dy1y1f6nIHjL0V6MPSxYkdr1mJo9scYx2+0nAQrefvo1m15B3K5YsUKsX7++KE1muvAExC227du3j\/oqBW+oCWguEHN\/IKRi\/7hxv6Y\/55DgNamhayq4ZMqCLCkGkStLieG\/XW2m9ndlX9N5m+xPSeimFOFVucprUJZBRPQX4tdm8myTn1x\/T+X+4Go\/jxs2AQreYfs\/id77LC2sCl\/8P1ankykLeGAgZcF1ApoLnNwfCF3bPxK6Ev4cIf7u++aeaBJcVeXEJkyYYH6SiHs22R\/x1EGaBl98ucELBOYYyC1FoZuq4IVd8svU\/fffL\/C178wzz2z1PlY1GLq+PwQZpGxksAQoeAfr+n51XArfm2++uXgwmNTMjUEg9wdCV\/YfJnSFEFhtu6ksWdmHOsHYRjmxUGMpd8ELDrIPEL3vP+mQ6FUZhaqhG4J7KikNur7s3bu3mKiMvN5TTz21WL68y62r+0OXfea5+0OAgrc\/vhx8T3zq8IaCl\/sDoW37q4Su6g\/TxSdUsYV89yOOOGI0AQ2\/yZSFFD4NV423Pgjeqmgu+pyS0E05wittk9cjKshgwi2+WiG9Af\/tYmv7\/tBFH3nO\/hKg4O2vbwfXMwpef5e39UAzEbqyNzaCF\/bv27evELeI7HY5Ac3FGzkL3iqh+8l77ytQHH\/88SKV1BHpG7nKGv7e1ti3GReqTTgOlRyQ7oAXui6ivSkysuHJfYdNgIJ32P7vVe8peP3d2cYDzUbsokdNgledgLZ\/\/37x8MMPi1NOOaUQWF1OQHPxRo6C1yQ\/Fz5BuTe8gJxwwgleyy+7cK06JifBK79MmCxYEZKR2lYb94dYtrNdEqDg5RjoDQEKXn9Xxnyg2QrdpgivbgIa0hggGhEBy03sor+5CF6XtAX0DS8k8Fsq0d4cBS\/Gic3yxP53hUMtxLw\/hLSTbZGAjgAFL8dFbwhQ8Pq7MsYDzVXo6gRv0wS0XARjladSt98kmts0CtVoLyaXdvliIiesweYYY7+JRdPvTTa1He1tsqepP\/ydBLokQMHbJX2eOygBCl5\/nKEfaGWxi8oLthsqNdz+xBOFIMFncWz4NI4oYXkCWuqCsanvKdrvEs1t6ideXBDtRX+7jPbmLnjbjvaGvj80jRP+TgIhCVDwhqTJtjolQMHrjz\/UA803qisf5LAHwmjDUUeNJqAdc8wxlVHBFAWjjVdSsj9ENLep711He\/sgeCXjNqK9oe4PTeOCv5NADAIUvDGoss1OCFDw+mP3faCFELpS5MoFRBDF\/ebEiUYrTaUkGF280bX9MaK5TRy6ivbK\/F11+XHYMn369M7KfpVZ2V6Pam4v6pGjfm\/IzdaekOdmWyTgS4CC15dgpOOxWtHixYvFtm3bas8we\/ZssXbt2jGrGkUyKflmKXj9XeTzQPNJX1BFB\/4fdUbVBUSaKjXInnctGH090IX9XYhcHac2o70Qts87eFA8f9++IlVGLuOLnOKU6jS7Xo94WUQJM\/QLojdUn1zt8b0ueDwJhCBAwRuCIttIggAFr78bXB5orlFdCFsZzW0SHaYrrnUhGP2pH2qhTftTEboqPzXaC5EGARpyk5U9kAv+xw89NHqx6qKmrUm\/XK5H2a4a7QVLLFjhu\/nY43tuHk8CvgQoeH0J8vhkCFDw+rvC9oHmEtWV0Vwsm4oNUSgIDkR06zaTKG+bgtGf9uEtxLY\/RZGr4yhfhDA50beSQ1Vlj0VHHx0s8hljLKBN2+uxiuXu3buDiPsQ9sRixXZJoIkABW8ToY5\/v+GGG8SGDRvGpC3ceuutYunSpWLdunVi5syZHVuYzukpeP19YfNAsxG7umguBC6WTDX93ErB6+bfXERuuXdqtBely1Bb2WbD8YjkYkxDOKNGM8aaHG\/qhDWbdtvc1+Z6rLNLfclE\/5Hm4LI8cSh72mTIc5GAJEDBm\/BY0Ildae6dd94pFi1aJJYsWSLmz5+fcC\/aM42C15+1yQPNJoUBD1rMHlcnoEHkNkVzdT0xSWuIHSH1J1zfQij760QuLPg7h\/Jwsfte1b5kIpeJrlueWLcYybHHHluIXBwvN3XBia76ZXJek+vRpB25j2wPf3dZnji0PTa2c18S8CVAwetLMNLxmLT2wQ9+UHzoQx+qjOIi0nvVVVdx0tqvfUDB6z8Ymx5oJlHdpglorlZS8NaT65PILfdUzb1FtBdpDlLANi1GoqOWQ3QXdjddj67XklrCDLm9ptHeWPa49oPHkYANAQpeG1ot7msieBHlvfLKK8W1117LKg1CCApe\/wFa90CrE7syZeHAgQNjormhZ703pTWEipD6k3Rrwdb+PotcHUHJB78haos0BXUxEkRz6yLAss2hC14pplHJAdeuabSXgtftuuZRaRCg4E3DD1orEL3FtmzZMu3vSHnYvHmzWLVqVZELOfSNgtd\/BFQ90KrErm4CGtIVYs16b4ry2gpGf2JhWzCxf2giVyUso7kICGBBEgg1KXJNlyjOJZ0hZoRXZWqzYAUFb9jrna21S4CCt13eVmdTa\/GuXLlylKsLobtixQoxbdo0TlxTiFLwWg0v7c66B5oqdrE0sE05MX+LxrbQd8ELQYfrHp\/rJ0+ePOr8kEUuIOgmoIERXhDwm83yxLlEd9sSvDiPWsKsLtpLwRv6jsb22iRAwdsmbcdzIVd34cKFo6MpdPUgKXgdB5hyWPmBpordxx47NAENh8g6qS4T0HwsrUtrMImQ+pzb9tg\/UwB+DW8LBpvsw0Wn1K+SldPEM4NuH7aL6QQ02wUrKHirvdEU7aXgdRnJPCYVAhS8qXiCdngToOD1Rjhmksz48ceNGty1a\/doRSp1BTT\/M9q3kIPgVYVuuYd1wpeR3CeKqC2EFf6LKC5SFfByVZeyYLo8cU7pDG1GeNUxWhftpeC1v1\/xiHQIUPCm4wta4kmAgtcToDIr\/HnPe+6osbvv\/vkommtaM9ffkvoWqkRvChHeOrEreyVFb5PAxf6fvPe+4jCbz\/ax+YduX63CgLZNSpDpbGiK9uYU3e1K8EquumgvBW\/okc\/22iRAwdsmbYtz6ao0YBLbjBkzRrm8rNIwFigFr8UA0+yKyM79998vpk591ujXBx74VbQJaD7W5ix4x4nm1IZyuoIq5E444YQxNWV9OHZ5rEs5MRN7q6K9uUV3uxa8OL8a7ZVfdrBqG0qZpfLyazImuA8JgAAFb6LjgILX3jEUvPbMyhPQ1MiuYcqp\/UkDHJGq4K2K7jaJXJN8XESvUZnAdpJWANxBm2haAS3UycrRXiwlnNuWSkQVC8eghBnq9eKeQcGb20iivRS8CY8BCl5751DwmjPTlRNDKbGJE58xaiRlwQsjdaK365QGneCtErsmIlfn0abP9uajoL09TSeghbZIRnuf88gjYu7RRxdL6ua0pSJ41WgvxC\/KYE6fPj0nlLSVBBjhTXUMUPDae4aCt56ZFLl4iOIPojXqBDS1IsOBA48m\/8lSV6KsK8Erc3GfEuMOc4IqeNXfTas26LxqOknL\/ioKd0SslAUXCyF4p\/\/yl8WYh+jN5XN8SoIX3GHPjh07ChdggY+YNbdd\/MxjSKCOAFMaEh0fFLz2jqHgPZwZRK661C8e+LKcWPmhLwXv9u13ZfPJshzlbUvw1k0204nesmd8xK7aVorR3lAT0OzvAPojZO6umo+KsY\/P8qlvKQpe5PDipUFea2CJv5suT5w6c9rXXwIUvIn6loLX3jEUvIeYqSIX\/y+juVUroKnR3ZwEbznKG1PwmlRUgAfaFLw4XwrR3pSiueU7R7kygxSRTdeE\/R0o\/BGpCl6ZwyvtQ89NlycOT4ktkoAZAQpeM06t76WuslZ38tmzZ4u1a9eOWZWpdWMTOeHQBa\/PCmi5Cl4MPTXKG0rwmopbnF+Xi+tah9fnUuoi2tvWBDQfLrpSZGoOe8oRytQFr\/SLWsIMYpjRXp8Ry2NjEaDgjUWW7bZOYKiCt2oCms0KaDkLXjXK6yp4fQVu1WB3WWnN58JRo70ydcWnPd2xXU1Ac+lHU93d1COUuQhe+Aa2opID7keM9rqMVh4TmwAFb2zCDe3jzfjSSy8VkyZNEkuXLm2seYpavFhqmFHdw8EOSfDqorkQuJg97TIhJ2fBq0Z5TQSvjbitiuB2fNtoPD3EB0qYYQGHiRMn1q5S1tjYr9MmXFZAM2k75j5Ngjf1CGVOglfHMqcJgjHHIdtOgwAFb8d+UFMXkJ4AQTtr1qxKqyh4qx02BMErhQxKA8kJaBC5NtFcHcHcBa+M8uoE7xAEblUkFqIXTLAsL6JutltqE9Bs7DcVu7LNFCOUOQpe8KxbntjGh9yXBEISoOANSdOhLSl4x40bJ\/bt2ycefPBBsXLlSvG6171O4N\/KGwXv8ARv1QQ0RO5C5srlVpasPBI0l0vjFelaC7ex4YR2kC8Bpsv1pjwBzQarreBNMUKZq+BNkaXN2OG+\/SRAwduxX6XgnTt3rli4cKFYtmyZ+P73vy8uvvhisWjRosMEDQXvMASvzwQ0nyFdFo0pLj7hImzBZAjitsr3aqQW0V68LEEAq1sOE9BMx7ar2JXtdx2hVF9y4Rcs8hDy5daUY3k\/FwHeNUvXvvK4\/hGg4O3Yp6rghdhFhPeTn\/yk+OpXvyre\/va3i4suuqjI75UbBW+\/Ba9uAlqbxd3lA01dYlgSb1P8uopaaeuX\/\/fpz\/jcxhKQ0V786\/HHH1\/ke8PnBw8eLMqbySgw\/r0siHNh6St21X6q1Qdi56N29ZJr41cXwctorw1h7huTAAVvTLoGbZcFLw7Bje+GG24Ql19+uZgzZ06R1\/vsZz+7aI2Ct3+CN\/QENINhV7mL+kAbP\/4466Z0othXvNYZUT4f7P\/8\/v1FvioFr54chC3SpyDm8DINYQtWELm5M5OLTFgP3JoDYkcou37JtWHlI3jlsw2VHNAOXuRzW+rZhhX3TY8ABW\/HPtEJXpj01FNPiS1bthSVG7CtXr1avOpVrxJXX301qzRU+Cy3SWt40EF0yAlo+GSJT82+E9B8hrTugRZTsJraahpdBtOdO3eKjZMmsTZ1CW55AtqTTz5ZCFz8O6K9WCo29y1kdLfMImS0N4dorm4s+Ape2SbueRC+uS31nPv1MXT7KXg7HgFVglea9fOf\/7zI6926dav4+Mc\/Lu68807xve99j2XJNH7LQfDKaA4eHPgjV3sKPQHNdVibPtBCimBTMWvaJxnlNZ2kZdpujvuZTEDrYsGKGCxjil1prxrtdYlQpviSa+ML0\/uDSZuSJb4sVK0AadIO9yEBUwIUvKakIu3XJHhxWpnXu379+uLz4wte8AIK3owEr4zmHDhwYBTNlYsCuNTMjTQUi2ZDPtBi2tnUNph\/4aGHvEpyNZ0j5d9tJ6ClsDyxD88YqQx19thEKNuqsuLDz\/TYvtwfTPvL\/fpFgIK3Y39issiPf\/zjIpfuzDPPrLRGzes966yzKHgzELyI3qIOKh6O8J+M5qYczejbA23dY4+JvXv3jiZj9eGzfdVNIsQKaDlGe9sWu7poL15csaSu+pusmS2\/5KT6kmvzCOzb\/cGm79w3fwIUvBn5EHm9EFB4sJ1wwgniiCOOyMj6+KamkNIAcbtjx47CP1Lk5vSg69sDDQtS3P7EE0Xk+pFHHilyVnUlueKPzjhnMElZsD1zbtHeNlIZ6hjKawYvtFgEBl\/h8JKFLYeXXJvx0bf7g03fuW\/+BCh48\/dh0B6oK7+pDWMVOHU5Y3yeX758udi0aVOxW\/l3eSyqSqxZs6b467Rp08S6devEzJkzx9iMpZJRg1huWHhj\/vz51v3qUvDK3Lw9e\/YUkzFQVePkk0\/udAKaNcAepTSofS+vwobfcp+k1cYKaDlEe7sWuxhLMmXhl7\/8pcD9E5FemZeaWsqSyz1BPYaC15cgj++SAAVvl\/QTPDcmxWHBC1SFwGIYuk2KXQhYTKjDpiuXhn\/btWuXWLVqVRH5QKm16667bozohdhFJQophOX5lyxZYi162xa8ugloiOggygvBm+PDrq8PNCl6TRZgSPCyLEyKEc1t6qsa7ZVfKpqOaev3rlIZZP90E9DwBQFf4bChNF7K6Usufurr\/cGFBY\/JjwAFb34+i2oxBCiEqhrNLZ9QJ1xlZHjBggWFUJVRW0y0k8K5LJSrJuyh\/Q0bNljnKbcheJsmoOX+QMjd\/rqL42vKj1iAQaYHpR7ttZ2AFuMGIfNR8bkeKSEp1OvtIrprOgFNLWGGiG8Kq6SFGBd9vj+E4MM20iZAwZu2f1q3DmJz8+bNo6iszoBy5Fbuo\/77xo0btaJVFbP4\/K+LJptEmXV2xRS8VQ+6cgQn9wdC7vY3XTCq6MW+6md75F2nsrpYiAloTSxsf1ejvRC8iGB2tbUpdl1r5uJaQnoTju9LtLfv94euxjPP2w4BCt52OGdxFhmBvfvuu4sbNdIRsM2bN28kgHXpDKrgRWQX0eHrr79+TDqD3EeNDkPwqukMch\/XtIbQgtflQZf7AyF3+00utLLoTSXa20XKggmv8j5yeeKu6hy3JXZDrYAWcsEKF3+FPGYI94eQvNhWWgQoeNPyR6fWyBSD008\/\/TCBCxEMIYtcXExWU\/N3VTErUxGqBK+as1sleMvpEaZQQgle3YMOkVyTFdByfyDkbr\/pWCmL3i6jvW1MQDPlYrpfV7nQscWuy0uuCbPYyxOb2BBin6HcH0KwYhvpEaDgTc8nyVmkRlzPPffcXgpe3YMOAhcC32byWe4PhNztt7l4dKIXQm7fvn3Rl9vNJZrbxFNGe7Ff7FzomGK3rRXQco\/2Dun+0DT2+Xt+BCh48\/OZt8W60mPnn3\/+qOJC+QTq5LILLrigUvCqlRpySWmQk3FQWQETSyBuIXJNork6R+T+QMjdftuLQyd6y9HekJO0UpiAZsvIZP\/YudAxxK7pBDST\/tvsk3O0d2j3Bxu\/ct\/0CVDwpu+jzi0sV1PIfdJazAdd7g+E3O13uViqRG+oBRhSnIDmwqnpmBi50KFLj8VKWWhio\/s9x2jvEO8PLr7lMWkSoOBN0y+dWFU1WaxcNUFXNkxXlqw8Ia2qLJksZSY73VSWrCpXty6Ht60HXe4PhNztd71wqkSva7S3LykLLjxDLVgRUuyGmoDmwqPuGDXai69KU6ZMCX2KoO0N9f4QFCIb64wABW9n6NM8MaK3WD1NLgShm8im+7fywhNqxQdZ01dXvxf\/tmLFCiHr9ZpWaIC4xbZ9+\/YRSJ3g9ZmA5uKh3B8Iudvv4jN5jFycQteGabQ3xwloPsyqjjXlVXV8CLHb1ktuCH5IqUJlHKRVQfTazBsIcX7TNoZ8fzBlxP3SJUDBm65vOrNMilBpgC6\/N5WlhVWRK\/8\/1AQ0Fwfk\/kDI3X4Xn6nH1Ile7KeLXg45mtvE2yXa65uv29YEtKa+2\/6uRnsheLFgRWrb0O8PqfmD9tgRoOC148W9EyQgo70wbcuWLcXSviEmoLl0NfcHQu72u\/isfEyT6JXRy4MHDxaVCcAMNWmPOOKIIjKXanQuBBuXNmyiva5iN2ZevkuffY6R1yDuYUhzSGl5Yt4ffDzLY7smQMHbtQd4\/mAEIHxvvvnm4iGBmfVdLOeZ+wMhd\/uDDSYhRFXZMjCC2MWLFSp6YMN4o9Ctp18X7XVJYcgpZcF2XKqpWBhXSHPo4n5Wtpv3B1tPcv+UCFDwpuQN2uJFINTCEz5G5P5AwIN2586dxcN1+vTpPih6cSyivbc\/8YRABQL4Fv9FNBdL60KIIKq7f\/\/+4t9j16HtA1A12gt+eFGwFbupTkCL4R95P0HbKSxPnPv9LYaP2GY+BCh48\/EVLW0gQMEbZoik\/Ek1TA\/NWlGF1f\/zjGcUQvfYY48VEyZMOKwBl1xVMyv6uZesf\/0njz5qNEmrz9FcEw+rJcyQ29tVtJeC18Rb3CdVAhS8qXqGdlkToOC1RlZ5QKqfVMP1UN9SnbD6j+OOEz+pMUAXvYxtb67tI6o7+\/HHi8oEEFFVk7RynYAWwy\/gBF5g0na0V94PZDoPvv50JbpjsGWbwyBAwTsMPw+ilxS84d2c2ifV8D18ukUbYVVXsxdtMdpb7SVd+kL5iwLSHJAmghxp+EVO3uoqLz\/WmHNtt60FK+TLHyryqBOBma\/u6jke1zUBCt6uPcDzByNAwRsM5WENpfJJNWQPfWb2m1ZyQG6vzFUNaXtubTXl6cro5b333ismT548qrJCcVX9JUJGx0NHe4eUI53bdUR7\/QhQ8Prx49EJEaDgjeuMLj+phupZ6FzQJuErc1WR\/wvxhsluQ9vqxG5ZXI0bN0489dRTBaLQQq6P3ENFe0NfF31kzT7lT4CCN38fsge\/JkDB285QCPWQbcfap88SO2plukobBC+E3BC2KqFrIq76+EUhls\/VBStsXxLUVB7YJ79GoLQjNxLoGwEK3r55dMD9oeBtz\/k+D9m2rDQRVqFtqRO+Q4j2QuRim6MBa5MnjcP78EUh9Piqa8\/0RVSdgAbGzJFu00s8V5cEKHi7pM9zByVAwRsUp1Fjpg9Zo8YC7WQrrAKddkwzVcJXreTQp2hvXTTXdwJaimMsxpgJ0ab6IoooLRaswMYJaCHoso3cCVDw5u5B2j8iQMHbzWBIIdrrMwEtNjWd+MVktr1799bW9o1tV4j2dUI3RmQ9hTEWgldbbaCqAia1IScaaQoYb2rFi5SWK26LCc9DAhS8HAO9IUDB260r247ExRBWMQmWhS+ivfik\/MgjjxST2XLJ7W2K5kLIY5OfykOKq7bHWMzxEKtteV2A1e7du4vFUiZNmiSwYAWXv45Fne3mQICCNwcv0UYjAhS8Rpii7lT1STXkSWNPQAtpa1VbEL\/YsJBFDtFenwlooXky2qsnqrsu8LKBFw8I3xgvIKF9y\/ZIICYBCt6YdNl2qwQoeFvFXXsy+UkVD1nkEfpGlnKL5tp64v974gnxQyXaixJmKGXW1VY3+Qw2pZAnzWjvodxc5EnXTUBTxTCuRVyTXCmtq6uL5+2KAAVvV+R53uAEKHiDI\/VqUI3EVS0d23SCFIRVk40hf1dXtvvZySeLnRMnhmy+sq2mhSGkyPWdgBa6M218UQhtc4j2ZMUPdQW08ePHi6ZyYkNZOTEEY7bRPwIUvP3z6WB7RMGbpuvLS8c25XSmPAGtLcIyelkVjZMpEbb26MqF1bWRS2Q99BcFW65t7F91XTRdTzrbWOe4DY\/xHKkRoOBNzSO0x5kABa8zuugHNn1SzUVYRQelnKDLOrQ55kmH+KLQpn9NzhXzuuhyfMm+p3DPNvED9+kHAQrefviRvRBCpHDzlNFMzojWD8nyJ1XkquIzecyZ\/blfHG1F42KKqzZ9YPtFoU3bTM9VNQGtKWXBtH11vy5yoXGvxrZ9+3YXk3kMCTgRoOB1wsaDUiRAwZuiVw63CQ9zzBrHH3yOxWd7uaSp7+S2PAjYWymZ4b+2y8c2na2PedJNXxSamHTxu+6FAwIXubmxr4s2Kl9IkUuh28Xo4jlBgIKX46A3BFIQvHhw7Ny5s5gBPX369N6wDdGRsrBCHVpUIogh4kLYm2IboaJxQ8mTzmGSVvm6gLg1mYAWY3yGGl\/SNorcGF5im64EKHhdyfG45AikIHgBpQ+fVEM510RYhX7IhrI91XZco719SVlw8UtbaSGmtplcF6Zthd4vRLSXKQuhvcL2QhCg4A1BkW0kQSAVwQsYOX5SDeVEF2EV4iEbyv5c2jF9UchxAloMH3Q9ScvluojBwbRN0\/HFaK4pUe7XNQEK3q49wPMHI5CS4JWdyuGTaigHhBBWtg\/ZULbn2k7Vi0Ju4qpN\/m2PsRDXRZt81HOp4wv5xFiwQt2YstCVZ3heFwIUJvq\/CQAAIABJREFUvC7UeEySBFIUvBJUap9UQzkwhrBitNfeO3J8jRs3TkyYMEHIBQmQS45KGDFm99tbmc4RscdYjOuiS3rlOsdnnXVWYQ6rLHTpFZ7blgAFry0x7p8sgZQFL6B1\/Uk1pOPamNnfdiQuJJ8225IRxH379ok9e\/aIyZMnF5MmURqPy8fWeyL0GGvjumhzbDGa2xVtnjcGAQreGFTZZicEUhe8umgvPhHGLjkUyhldTLRp+qQaqm+5tVMXQTxw4EBR1xhiN6fx1ZUPfKO98rrACy3+gDsi6ois9+GFQzcBjZPSuhqtPK8PAQpeH3o8NikCuQheQPN9yLYFPpVPs0NYOtbEp6b5oLmML5M+t7WPTbRXXhd4uZDpI32qJW2am5vTPbetccTzpEuAgjdd39AySwI53nxtHrKWOLx2NxVWXiexPFgVcRAX+GQ\/hM3npSPV8ZWq35q+KFR95cACKrlvpiK33E9Ge3P3\/HDsp+Adjq9739McBW9K0V4fYdXm4BpKneNQ+aCM9tqPTvWLghSzWAJbpiz0MZrLCWj244RH5EWAgjcvf9HaGgK5Cl7Zpa6icaGEVZuDU41AQ3wgV7UP+ZIx86S7Gl9tjouQ54K4\/dnPflaI3EmTJo2Wv+5DxQvXaG5IvmyLBNomQMHbNnGeLxqB3AVvOdqrq3sZCl5MYRXKRpN2+lDnuM3IetMnexPmfd5HNwENL1QyTxfXZM7pC0w\/6PPoZd+aCFDwNhHi79kQ6IPglbBjTNJqU1i1PWhyrHPcZZ40o72HRij8oL4A4ktBOWUh5y8KjOa2fTfi+VIlQMGbqmdolzWBPgnecrTXZ5JWl8LK2okeB+RQ5zill46hR3tdJqDl8kWBItfjRsJDe0uAgre3rh1mx\/r4yc5lklZKwqrtkZhi9DLlPOkYXxPa9rnp+UJdF6l+Uejj\/c\/Ut9yPBJoIUPA2EeLvWRLoY7QXs8SxoEDdJK2UhVWbAymFygQ55Un3veSb7isHcnF9JqCl8kWB0dw27yw8V84EKHhz9h5tryXQx2iH7pNqTsKq7SHbdrQ3VASxbU7yfC5fE7qytem8Ol9A4I4fPz7o6oZtjzH0O1eR+8ADD4jFixeLZcuWiblz5za5kL+TQFACFLxBcbKxFAn0LdoLxvfdd5\/Ys2ePmDBhQjHhRjfRJkVfdGFTG9HePuVJ5x7thWjH1xB1BTSIXJ9obtO4bWOMqUI3x5q5WJVu+fLlYtOmTWL9+vUUvE2Dir8HJ0DBGxwpG0yRQF+ivaqwOnjwoDj22GML3CeeeGLW5ZLaGDOhI3G5R3ObmOcU7U3lK0foMZZzNFcdX3feeadYtGiR2LVrV\/HPFLxNVx9\/j0GAgjcGVbaZLIEco71NwirGQzZZB3oaFqIywZDypMvluFJazrnpuvAcKs6Hh4r29uUlXYrdefPmiTe\/+c2F8F29ejUjvM4jjAe6EqDgdSXH47IlkMuDxEZYhXrIZutUS8NtKxOkEkG07Gaw3VMqxxVjAlowUEpDLi+iuebmmvKT4peC15QY9wtJgII3JE22lRWBFIWvr7Byechm5bSAxjblqqYaQQyIwKqpcrS3zeWc25qAZgXEYGeTLwp9F7kqJgpeg0HDXaIRoOCNhpYN50Kg6zSH0MKK0V67kVfOVZ04cWIx6Qkl4LBhQmDuS8raEanfu81yXOWvHCjJF3sCWkhWsi3dF4UUX7hj9J2CNzZVtm9KgILXlBT36zWBLh4+sWf2M9prPmSliLv33nvF5MmTWfXCAF2sxRfkdQGf4I984cCLCP4\/1w39mjVr1sj8HCst+LJnhNeXII\/3IUDB60OPx\/aOQOxob+hobpMDTD6pNrXR59\/LEcRx48aJxx57rOgyK180ez5UtFdeFyhdpZYTg8hFVDf3rfxCHfs+kyovCt5UPTMMuyh4h+Fn9tKCQIxor80ENAtTjXe1naRl3HCGO5rkSceKXmaIy8hkV15VvsAqaLlvTbm5Me4zqTOj4E3dQ\/22j4K33\/5l7zwI+EZhTISVh3nWhzZN0rJuMKMDXCLroaKXGWHyMhWMd+\/eXSyEUhcdd\/GFl2EtHtwkcnWm+N5nWuye96koeL0RsgEPAhS8HvB4aP8J2EZhcniY57SggO8IC5EnzVxoOy9U8VK\/cqBFpCogZSHmCmh2lrvvbXufKJ9pKKKXgtd9jPFIfwIUvP4M2cIACDQ9kEIIqzYxdlliKnY\/Y7x0sPKFnddktBcvV6iqcOSRR\/ZqAhpouERz7ShybxIggZAEKHhD0mRbvSZQjuLEEFZtA0xpQQHfvreRJ81ob7OX1AlocinZSZMm9WISIEVus\/+5BwmkSoCCN1XP0K5kCeCht2XLltFscpRKyv3TrOuko66d1EWeNKO9eq\/X+eK+++4rIry5Vr7wTVno+jrh+UmABISg4OUoIAELAnjwbd26VcyZM6c46o477si6Nqja9VwmaaUSWWe0VxQT1DBusFCHrJkrc3PL5cRy48VorsWNkbuSQAYEKHgzcBJN7JZA3YOvKbe3W8vdzp6qMEkxT3qodY51vkApsaYJaDlExxnNdbtv8CgSSJ0ABW\/qHqJ9nREwffCZ7tdZRxxOnIowsYkgOnQz2CFDqHOs8wUELial2S4OkdpLFaO5wS4FNkQCyRKg4E3WNTQsNwJ9FL5dCZM2JqCFHl99rXMsUxbUFdAgcpuiuU18u36pyknkynJechLg+eefL5YtW1aL+IEHHhCLFy8W27ZtG7Pf7Nmzxdq1a4sltLmRwJAIUPAOydvsaysE+pbm0JYw6WICWowB0Yc6x236ou2XqtxeTMu1a6WQnTt3bq3oZc3bGFc328yZAAVvzt6j7ckSyO2hagIyhjDJJWXBhI+6T451jrv0ReyXqpyiueo4OnDggFi+fHnxT6tWrSrSR7DdeuutYunSpWLdunVi5syZ2uGJfa666ipGc20vXu7fWwIUvL11LTuWAoE+R3vxSXvKlClOmFOcgObUkYaDcqhz7DoBLQav0C9Vub94ymjuggULxPz580fIq\/5d9ckNN9wgNm\/ePEYox\/AZ2ySBXAhQ8ObiKdqZLYHcH7o68C6TtLqMIHY9eFKrcxxyAlpotr6VL3KN5uo4VqUlNKU1yMjw3XffLVADWeb+zps3jwI49IBle9kQoODNxlX9MdR0MoW8aW\/atKnofNVkC3y2W7NmTbHPtGnTtJ\/58Hlv4cKFI4grV64cEzFpg26fo72YpT916lQtxhwnoMUYDynUOS77An4LMQEtBi+bl6o+iVyVZZXglfdG3O90k9fkPfb0008fCVxVBHPSWowRyzZTJ0DBm7qHemifyWQK3Q0dwhbCVb1Z498QvZD5bfiMd911140RveV8N3n+JUuWdCJ64dLt27f3xrO6SVptTnrKDWToz\/ZN\/c\/ZF02VL\/r49SSE4K0aE13e+5rGKX8ngdgEKHhjE2b7hxEwmUyhE67lvDUZtV2\/fr3AjGVsZaFc9ekP7W\/YsKGzCR19jPZCyO3Zs0dMmDChWIELSy5Xrbo19Msi9iQtmbKA60EtJ4YlsG1r5qbgKzXae84554xM6tOLo46za0pDlc+aUiFS8DVtIIFYBCh4Y5Flu5UETCZTlCO3sjH13zdu3KgVraqYhQBbtGiRWL169UgUoy2TKHNsF\/YlOqVOepJCF+xOPPFEgdW3uFUTCB3t7etkQDVlQdLsu9hFP30mrelGHQUv70ZDJkDBO2Tvd9B3k8kUdflpalrD9ddfPyadQXZHjQ5D8OrK96T0aS\/HaG\/TBLTUJml1MNSNT+kb7W3yhbEhCe5Y91KY43Vji7jqXthUlqzq\/pbCi74tA+5PAqEIUPCGIsl2jAiYTKbAJBrUntRNyFCjt1WCV30YVAlek7I+Rh0KtFMu0V6bCWgpTNIK5J5WmrGN9qq+gIEyfcR3BbRWOltzEpsJaLlcNz5MZeqWnGhrGqVFcAATfmWtXt2918cuHksCuRGg4M3NYz21V41InHvuuYMTvNKtKUatfCc92Qq5ng5xo241leSSvsDLBP4gTxoCF7m5+P+cNx\/xmuJ1E9IXTUsLV0V0ESBYsWLFyBSTJYlD2s22SCAlAhS8KXljwLaoUYsLLrigUvD2MaWh7HafB3+oIRT6M7nvZ\/tQ\/cqlHXWSlsyD7ssENNUHNtHcXHxHO0mABNIkQMGbpl8GZ1X5M90QJq01ObkL4Rt70hOjvU1eP\/Q7Irg\/+9nPiioLkydPHkVzc58ISJFrPga4JwmQQDgCFLzhWLIlAwKmkyl0ZcN0ZcnKE9KqypKVl+bsuiyZAarRLrE\/14aO5jb1jdHeakJVvpClxZC+kKvg7eIFrmks8ncSIIHhEKDgHY6vk+mpyWQK3QSL8sITupWDdPV7ZR6brNebUoUGU6fEEAs2E9BM7bTZj9HeQ7R0kXUIWzkBTf0dk9OmTJmSRc4uo7k2VwT3JQESiEmAgjcmXbZdScBkMkXflhYOMRx8o72+E9BC9EFto2mSVujzpdSeywQ0uaod+pFqnWOK3JRGGW3RETCpBU9y\/SNAwds\/n7JHPSdgG+1tO2XBBb86SQvRyxxXAzPpN3yhvnS4rkaXYp1j23Fpwov7kEBoAvIL35w5c0ZL0oc+B9tLkwAFb5p+oVUk0EigKdobewJao4GWO6jRXgjeqVOnWraQ7u5VkXWffNwU6hznEs1tKuuV7sihZSEJ4KvhJz7xiaLJ\/fv3U\/CGhJtBWxS8GTiJJpJAFYFyVC2HaG6TN+Vne1lj1kcUNp0r5u9t+aKLXOicornl1cVMF26IOTbYdjcEkMqwY8cOMWPGDLF582YK3m7c0NlZKXg7Q88Tk0A4AhAgW7ZsKUpYQSjiDxYjyHXVrVwnacGjTRPQwnn9UEttVL7IJZqr8pXzAPBvq1atEljFEVvT0rwxfMQ2\/Qng5eWSSy4RV199tZg5c+aYBpui+HjRueyyy8RHP\/pRceONN1Lw+rsjuxYoeLNzGQ0mgbEEIES2bt0qkJOG7Y477shiBr+JH3OYpCVFLmzFZ1J1BTQIrDbzkUNHe3MUueq4qlpCPLWlxU2uhaHvI3123333jZZLlkxMovio8nPOOeeIuXPnCk5aG+ZoouAdpt\/Z68wJ1AmRptzeHLue4iQtcJQiV0bWIW4hcruMrIeI9uaUslA3nstCSO7LtIa87gKIyC9cuLAwetq0aWMEr0kU\/6STThKLFy8W27ZtG9PxefPmMa0hr6HgZS0Frxc+HkwC7RIwFSKm+7Vrvd\/ZUpikJaO5iORC5EJcylxjpJDg\/1PZbKO9uUdzddyrBG95gZpUfDYkOxBlxTZ\/\/vxRt+WLyLJly4pILDYpdleuXCme85zniPJiQy5RfEZ4hzTSDvWVgneYfmevB0Kg79HetkqYtTUBLfSwNKlz3MeXI8mRgjf0iArXnvQNoqwQuDqxWz6bLvfaJYpPwRvOjzm1RMGbk7doKwk4EOijoAnx2d4EZW6l3ar6VK5zfNZZZ4123b59uwmKLPdxEUNZdjRTo6V\/8OKK3NzVq1ePIru6LtkIXkbxMx0UEc2m4I0Il02TQEoE+ih8bT\/bm\/hDF81FTm7bE9BMbLXZp48pC039d\/nc3dQmfw9LwGapdwresOyH1hoF79A8zv4OnkDf0hxCRXvRDgR0X0q7yYGue9Hp48uP7sKuivKxLFkat8GYEV5OTEzDxylZQcGbkjdoCwm0RKCPgscl2itTFjAhTi0nltoENNthYRrN7dvLT9VncMzwx6QnTJCiELIdTXH2D5XDyyh+HP\/0sVUK3j56lX3KgoBaagcGywdym8b3TfCYTNKSKQuI\/qnlxCBy26yZG9rPpiK3fN4+vvyU+9i0KEFoX7C9ZgKYOIaqC7IaA45omrimi8wzit\/Mmns8TYCClyOBBDogUL5x2+SxhTa3j4KnPEkLQlZGc8vlxHJduliOgz76L\/QYZ3v9IFCViqKWLmMUvx++jtELCt4YVNkmCdQQqPqkiojHhg0bxNq1a8XkyZNbZ9jHaO8vfvGLIoo7YcKEUc1ciN+hRnNbH1Q8IQkEJFCXe80ofkDQPW2KgrenjmW30iVQVSqp6t\/b7ElfooXqBDSV34knnihyjuj2xT9tjmmeiwRIgASY0sAxQAIdEKiKUnSZ1lDGkGO0t24CGvqH1dH27t1b5Omi7mdKq6LVDUPX3NwOhjZPSQIkQALJEmCEN1nX0LC+EqgSvFWzjbvikEM00XYCGiox7N69u0CacrSXIrerUc\/zkgAJ9JUABW9fPct+JUsgF8ErAaYY7fWdgKaWMJs6dWoy0d4cXjKSvbBoGAmQAAnUEKDg5fAggZYJ5JDSUEaSghDTrYDmMwEN0V4sZ4p2u4z2Mprb8gXI05EACQySAAXvIN3OTndJIOVJa01cuoj2ymgu8m+xIfcWE8+w3G+IzWXBihDnTeElIkQ\/rrrqKrFmzZqiqWnTpol169aJmTNn1jaNiiQrVqw4bJ8ualGHYMA2SIAE0idAwZu+j2hhzwhU5ep2XZbMFHMbQk0XzYXAHT9+fJTFIUItT9zEsG\/RXIjdXbt2iVWrVhW+wRi+7rrrGkVv+bgmbvydBEiABHwJUPD6EuTxJOBAQEa41q9fX6w0lFKFBtPuxBC+SDNANQV1BTQIqVDR3Ka+xYj29k3kSoay2L8cw\/j3qlWvVO5yn7PPPrtY6pcbCZAACbRBgIK3Dco8BwloCKSwtHAIx\/imOVRNQMPiEF2UDgsV7Y3xQhDCX6HaqPoi0fSlomn52FD2sR0SIAESUAlQ8HI8kAAJeBOwFXehJ6B5d0DTgEu0t6\/RXB3fqrSEprQGOWkTtZC3bds2alqNFMfwJ9skARIYNgEK3mH7n70ngaAEmqK9sSegBe2MEEUFB1RyQKoF0iog0srbkESu2vcqwVu3\/CuOL6fz4N\/k1w5OWgs9gtkeCZCAJEDBy7FAAiQQlEA52tv2BLSgnfl1Y8gphvBFigVEL8qh2Ua1Y9jVZZuugrfKZrQH4bt27VoxefLkLrvGc5MACfSQAAVvD53KLpFACgQgCLds2dLZBLTQDCDcZ82aNWp2+\/btoU+RZHvlEmKy9Ng3vvGNMRUapPFNKQ1VnXQ9LkloNIoESCA5AhS8ybmEBpFA\/gQgdrdu3SrmzJlTdOaOO+7oZAJaKJLlaO7Qo7vg6jppjYI31KhkOyRAAjYEKHhtaHFfEiCBSgJ1uaxNub0pYjXJzc2xX6FY63J1TcqSVaUusDZvKM+wHRIgAR0BCl6OCxIgAS8CptFO0\/28jPE82ETklk+RQ788sWgPl+L27rvvHuXdmqQlyJrT8+bNE8uWLSva1k1ki2Ez2yQBEhguAQre4fqePSeBTgikGBUdqmgNMQCalhbWRXRlLV5Zlsx0SeIQ9rINEiCBYRKg4B2m39lrEuiUQAoC0yWa2yk0npwESIAESMCZAAWvMzoeSAIk4Eugi2hvCmLblxuPJwESIAESsCNAwWvHi3uTAAkEJtCGAGU0N7DT2BwJkAAJZEaAgjczh9FcEugrgdDClyK3ryOF\/SIBEiABewIUvPbMeAQJtEKgPLFHnnT27NljVqOSs+U3bdpU7FL+XR7XNLkI+8klXuUxXSz16pvmEFo4t+JsnoQESIAESCAqAQreqHjZOAm4E5Dlm1avXi3mzp2rbUhX91Q3K75c41RXPqpcV1Wef8mSJWL+\/PnuHXE40la0MprrAJmHkAAJkMCACFDwDsjZ7GpeBCBAIVTXrl0rJk+erDVeJ1xlZHjBggWFUJVR2\/Xr14+Ec1koy2MgrGVtVJywajWttkjWRXv7JnLBevPmzWLVqlVi\/PjxbSHmeUiABEhgEAQoeAfhZnYyRwImAqhqdSr13zdu3Cg2bNhwmHBWxeyePXvEokWLRDmabBJljs12CMv6ypcSLMZAwRt7RLF9EiCBIRKg4B2i19nn5Amoq1jdd999YteuXYXNqiCqW8ZVTWu4\/vrri+PLQkqNDkPwLl26VKxbt07MnDlzxKfLtIayk2zTHJJ38q8NVHOrKXhz8RrtJAESyI0ABW9uHqO9gyAgUwxOP\/30kVAtL+WKz97Lly8XWKVKTUMAIDV6WyV41ZzdKsFbTo8YBPwWOwmxi8mGeNH4xje+oX0xadEcnooESIAEekuAgre3rmXH+khAjbiee+65FLw9cnJVekqPusiukAAJkEBnBCh4O0PPE5OAELrSY+eff\/5hEVvJSp1cdsEFF1QK3j6mNPR9vFDw9t3D7B8JkECXBCh4u6TPc5OAJYFyNYUhTFqzRJTt7hS82bqOhpMACWRAgII3AyfRxOERqJosVq6aoCsbpitLVp6QVlWWTJYyk8S7LkvWB8+D4YoVK0ZdQc51eXIgfqTg7YO32QcSIIFUCVDwpuoZ2jV4AuqEJlRO0E1k0\/1beeGJ8mQ31PTV1e+VwkzW602pQsMQBgMF7xC8zD6SAAl0RYCCtyvyPC8JGBAoRwd1+b19W1rYAEsvdxmy4L3jjjvEX\/7lX4qzzjpLYDnrZzzjGWN8\/NBDD4mPfvSj4p577hGf+cxnxCmnnNL5GPjud79bvDheeeWV4rHHHhOLFy8uFnYpV0wJYejBgwfF5ZdfLl73uteJ3\/md3wnRJNsggcERoOAdnMvZYRIggRQJDFnwwh\/f\/va3xYUXXig+8IEPiPe85z1i3LhxhZsef\/xxcd1114m\/\/du\/FWvWrKlcZrtNn6I2NtKEzjvvPPFbv\/Vbo68vsQQv+nb77beLT3ziE8XiMKeddlqb3eW5SKAXBCh4e+FGdoIESCB3AkMXvBC2iJb+0z\/90xhhi0gqRPAll1wi\/uRP\/mQkhLv091e\/+lWxdevWUY3sqqW5Q9oo+Tzzmc8U73vf+5LgELJ\/bIsEYhOg4I1NmO2TAAmQAAkYEbj33nsFyu1hQ+rC3r17i1SHl7zkJeKyyy47LNXBqNHAO\/3iF78QSC1697vfLf7gD\/6gaF0neJGmAWF65plnFukITz31VJH2ANGORWOuvfZagUjx61\/\/evHhD39YHHnkkcXERUS6Tz75ZPHBD35QoNb2UUcdNeoBFov5q7\/6K\/GlL31JzJgxI3DP2BwJ9JsABW+\/\/St++MMfir\/4i78obqpVDwyZP5fKQwX5eh\/72MfE7\/\/+7xf5aroJViHdBkZf\/vKXi+jSxIkTQzbNtkiABCwJIKL73ve+V7zlLW8pBCEEJlIanvWsZ1m2FGf3f\/u3fyvuFai08exnP1sreKVwP+aYY4oUBOQcS1H84IMPihNOOEG8\/e1vF7jXfeELXxDTp08v8oDPPvts8dKXvlR861vfKv5glcRzzjln1JH777+\/ENp\/\/Md\/XBzPjQRIwJwABa85qyz3RFTh7\/\/+78VHPvKRIv+r\/ElQTga57bbbxBe\/+MUkogZYYvWWW24pHiqIhMQWvPJT4XOe85ziISJzB7N0OI0mgcwJ4Hr8m7\/5G\/HpT3+6eAHF\/yNPVt1kBZFdu3Yd1tvyxE5ETZH7KzdZhQR\/x71l8+bNo9SEJnSw7YorrijE66pVq8SECRMOE7yLFi0q8nv\/93\/\/d8wEOyl4EbFVBTxetj\/+8Y8Xk\/X+9E\/\/tLj\/7N69u4gGY5IaIr1yQ5uXXnppIY5xfzz22GObTObvJEACvyZAwTuAoSBFLSIna9euFb\/xG79R9BpiGLlon\/zkJ7UPlS7QIJqzZMkS8f73v3\/0kIsteNFPTAjBS8GnPvWpJER\/F+x5ThJIhcDNN99cpAMgMooX8TPOOOMwwYuc3quvvlqgZJ9uk9VL8BvEKV6edXWsbQQv7qUXXXRRYY9ajUGK2Re+8IUCFRW2b98uPve5z4mpU6eOTJP7vOIVryjqMssXa9zfkL6BiPGsWbOK\/R9++OFiFUVEfstVHz772c+Km266qbiXo8QgNxIgATMCFLxmnLLfS6YtnHrqqUXkZMqUKQL5YIiGvPOd7yw+Iaq5Yl11GAL8O9\/5jsBNfdKkSYUZZcErI0Bf+9rXiqj1G97whiKKvWHDhiIVAtEctIE+Ijryh3\/4h+Jf\/uVfxF\/\/9V+Ln\/\/856OcOXWms3w4Pv\/5z+eEkK6cz\/OSgBBFCgNeehHl3L9\/vzYdC8K1SfDi\/oboblkYqlHdjRs3WkV4pWhF5BWCXG7lJcIRmS5XlKia2KZ7oS8vDKMOjDYCAByIJNBHAhS8ffRqRZ9k2R+I2\/nz5xeRipNOOklb97ILLDJ68uIXv3jMw0S9wUOQ\/uM\/\/mMxCQQ5yW9605uKSImsh4mIx+\/+7u+Kl7\/85cXkjxtvvFG89rWvFb\/85S+LXOY9e\/YUOXOIciO3Tq33iTbwuRMPqlTyBbvwA89JAl0RgNDDlxbk1SNC+qMf\/UibjmUieCF2sdXVxbVNaSivYlgWvIjMIjUB95BHH320CC6gqgI2Ct6uRhXPSwJPE6DgHdBIQGQU6Qtf\/\/rXxQte8AKBCRCf\/\/znR5\/RJIpyzpv89\/KSqOVFEXCjh5BWb+542KA2pcmGtALkvyGtQM3Zk4IXuW7YBwXokcem5iNLW3R5cMizK+fMfeUrXynqej73uc8dmSYn+GEGND47ciMBEmiPAFKs8DKLaxsvtLi+ISBxvZfTsapyeOU9CqkASAnAJDB5T9L1xFXwVkV4ZR1e+fUMQQU5L4CCt72xxDORgI4ABe\/AxoWcPfyDH\/yg+MSPSgjlzSQyUl6+thz5kH+3EbyY\/YyHHXLZ1Lw8KXjx4MAEFnzKXLhw4ZjJZbo8uKqcOTw88TAsnwfpDhDcKCZf95Ac2JBhd0mgFQLyhRMpSLgPyAlZunSspgivTAkILXjlVyiUBFPzcMtiVgYXkDKBlAqUJgsleJnD28pw5El6SICCt4dObepSWazaCt6qh42aM4c2McvYRvBCtKIMj1ruB+2UI8lvfvObi5nSmIQiN11eW9UDBnZiFnVZ8LZRPL7JN\/ydBIb\/TMX\/AAAEnUlEQVRIQL6IP\/HEE9oSZDIdCzV5kd971113Nebwmry420Z4EYXGV6QdO3YUX8tkSpTu3oF9UEIMKVq4X0GE65YftsnhlVUaMEYQBUfZM24kQAJmBCh4zTj1ai9fwWvykHCJ8KJdTDzTTTJBXh9yjxEpwapLKNquRqcpeHs1RNmZARFQU610JciAQpYO\/Id\/+IdiQitSF1wnrakv5sjxt6nSAFvwJQqTZXGfkos\/VL0so8QiFpXA\/QrpDr6CV5Yre8c73sGvUAO6RtjVMAQoeMNwzKoVE8Gr1q2UnZP1LU0iJ66CtyrCixxcRGSRmwfx+5Of\/KQoVyQrLVDwZjUEaSwJeBFoSmlA47qyZOXUK5OX97KhWAwDFRr+7M\/+TMybN8+rH7YHo1wbosVcac2WHPcnAU5aG+QYMBG8AFM1uzmW4EXkBGJWrUcJO8pi9r\/\/+7+LSAmW3cSkEJRTCyF4mcM7yMuBnc6QQN3CExChsu4uuhZy4QmJCuUTERkuV3qJiRL1fXF\/xOpuENxcICcmbbbdRwKM8PbRqw198hW8JlERlwgvJqdg0hgiGFhSWG5lMYs8OpQWwx9Z6zKE4MWkGeTcIXKM5T25kQAJkICOAKK8mAeACa7lVeBiEUOFGpRivOaaa7g4TizIbLfXBCh4e+1efed8BW\/V50T131Hf13bSmpwB\/Zu\/+ZtjIhg6MYsHDlZjQ3QXtS7\/9V\/\/tZjsok5Es520hrJnWL8egvfkk08e4Mhgl0mABEwJoNqLrP+NhSZibojuYpLaq1\/9am1lnZjnZtsk0BcCFLx98WTAfpikLJRFc3llIJcIL7qAB8g3v\/nNovyYLNgesGuVTUFsI2Lzohe9iJ8L2wDOc5AACZAACZBAiwQoeFuEncupqhaegP1YiUwuJBF64Qm0j2VFUXrowgsvFOecc05ryLZu3Vrkx2F1JznzurWT80QkQAIkQAIkQAJRCVDwRsXLxl0IoJTPLbfcIq688soxtXZd2jI5RpY8mjJlinjPe97DySAm0LgPCZAACZAACWREgII3I2cNxVSkF3zsYx8rctXUyWux+o\/JalhiGZHtNtMoYvWH7ZIACZAACZAACYwlQMHLEUECJEACJEACJEACJNBrAhS8vXYvO0cCJEACJEACJEACJEDByzFAAiRAAiRAAiRAAiTQawIUvL12LztHAiRAAiRAAiRAAiRAwcsxQAIkQAIkQAIkQAIk0GsCFLy9di87RwIkQAIkQAIkQAIkQMHLMUACJEACJEACJEACJNBrAhS8vXYvO0cCJEACJEACJEACJEDByzFAAiRAAiRAAiRAAiTQawIUvL12LztHAiRAAiRAAiRAAiRAwcsxQAIkQAIkQAIkQAIk0GsCFLy9di87RwIkQAIkQAIkQAIkQMHLMUACJEACJEACJEACJNBrAhS8vXYvO0cCJEACJEACJEACJEDByzFAAiRAAiRAAiRAAiTQawIUvL12LztHAiRAAiRAAiRAAiRAwcsxQAIkQAIkQAIkQAIk0GsCFLy9di87RwIkQAIkQAIkQAIkQMHLMUACJEACJEACJEACJNBrAv8\/4Odgg\/n8L0IAAAAASUVORK5CYII=","height":337,"width":560}}
%---
%[output:4fdbba66]
%   data: {"dataType":"text","outputData":{"text":"\n--- Part (c) Results ---\n","truncated":false}}
%---
%[output:91609af9]
%   data: {"dataType":"text","outputData":{"text":"1. Local time at theta = theta_0 + pi\/2 is: 23:02:03 (next day if past midnight)\n","truncated":false}}
%---
%[output:60d354e5]
%   data: {"dataType":"text","outputData":{"text":"The number of iterations is 2","truncated":false}}
%---
%[output:3e42949e]
%   data: {"dataType":"text","outputData":{"text":"\n 2. State in ECI 6 hours later:\n","truncated":false}}
%---
%[output:3b1afd5c]
%   data: {"dataType":"text","outputData":{"text":" r0 = [2073.667, 8266.364, -3628.842] km\n","truncated":false}}
%---
%[output:936e1791]
%   data: {"dataType":"text","outputData":{"text":" v0 = [-6.11264, 1.67138, -1.26311] km\/s\n","truncated":false}}
%---
%[output:6b30419d]
%   data: {"dataType":"text","outputData":{"text":"3. State in SEZ 6 hours later:\n","truncated":false}}
%---
%[output:9e4e0fd0]
%   data: {"dataType":"text","outputData":{"text":" r1 = [-2909.901, 5079.220, -13550.614] km\n","truncated":false}}
%---
%[output:78f51728]
%   data: {"dataType":"text","outputData":{"text":" v1 = [-1.87996, -4.51006, -3.22634] km\/s\n","truncated":false}}
%---
%[output:6e1601cf]
%   data: {"dataType":"text","outputData":{"text":"\n--- Part (d) Results ---\n","truncated":false}}
%---
%[output:9d111e2f]
%   data: {"dataType":"text","outputData":{"text":"New Spacecraft Velocity after release:\n","truncated":false}}
%---
%[output:23b00857]
%   data: {"dataType":"text","outputData":{"text":" v_sc_new = [-6.15019, 1.22850, -1.15081] km\/s\n","truncated":false}}
%---
%[output:65f0f99f]
%   data: {"dataType":"text","outputData":{"text":"\n--- New Classical Orbital Elements (Post-Release) ---\n","truncated":false}}
%---
%[output:2405e50f]
%   data: {"dataType":"text","outputData":{"text":"Semi-major axis (a): 8778.7000 km\n","truncated":false}}
%---
%[output:9771e70c]
%   data: {"dataType":"text","outputData":{"text":"Eccentricity (e): 0.061268\n","truncated":false}}
%---
%[output:4334c999]
%   data: {"dataType":"text","outputData":{"text":"Inclination (i): 25.2827 deg\n","truncated":false}}
%---
%[output:67339525]
%   data: {"dataType":"text","outputData":{"text":"RAAN (Omega): -168.4358 deg\n","truncated":false}}
%---
%[output:25b18ce5]
%   data: {"dataType":"text","outputData":{"text":"Argument of periapsis (w): 90.8535 deg\n","truncated":false}}
%---
%[output:2b2c80a5]
%   data: {"dataType":"text","outputData":{"text":"True anomaly (nu): 155.6793 deg\n","truncated":false}}
%---
%[output:4c056766]
%   data: {"dataType":"text","outputData":{"text":"\n--- Part (e) Results ---\n","truncated":false}}
%---
%[output:044c00ed]
%   data: {"dataType":"text","outputData":{"text":"Radius of new pericenter (r_p): 8240.8451 km\n","truncated":false}}
%---
%[output:17d4a4ba]
%   data: {"dataType":"text","outputData":{"text":"Velocity at pericenter (current): 7.16466 km\/s\n","truncated":false}}
%---
%[output:32a6778c]
%   data: {"dataType":"text","outputData":{"text":"Velocity for circular orbit: 6.95477 km\/s\n","truncated":false}}
%---
%[output:9d4f7384]
%   data: {"dataType":"text","outputData":{"text":"Required Delta-V to circularize: 0.20989 km\/s\n","truncated":false}}
%---
