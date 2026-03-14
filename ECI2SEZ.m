function [r1, v1] = ECI2SEZ(r0, v0, rS, phi, t)
% Earth angular velocity (ignoring solar/sidereal difference -> 24h)
psi = 2*pi*t/24;
omega_E_vec = 2*pi/(24*3600) * [0; 0; 1]; 
    

ECI_R_ECEF = [cos(psi), -sin(psi), 0; 
              sin(psi),  cos(psi), 0;
              0 ,        0 ,       1];
ECEF_R_ECI = ECI_R_ECEF';
    
R = 6371 + rS;
lat = phi ;
lon = 0; % Assumed 0 as in part a)
    
% Station position in ECEF
rS_ECEF = R * [cos(lat)*cos(lon); cos(lat)*sin(lon); sin(lat)];
    
% Rotations to build ECEF to SEZ matrix
Rlon = [cos(lon), -sin(lon), 0;
        sin(lon),  cos(lon), 0;
        0,         0,        1];
Rlat = [sin(lat),  0, cos(lat);
        0,         1, 0;
        -cos(lat), 0, sin(lat)];
ECEF_R_SEZ = Rlon * Rlat;
SEZ_R_ECEF = ECEF_R_SEZ'; % Inverse is transpose
% --- Position Transformation ---
rP_ECEF = ECEF_R_ECI * r0; % Spacecraft in ECEF
rP_SEZ_rel = rP_ECEF - rS_ECEF; % Relative to station (in ECEF basis)
r1 = SEZ_R_ECEF * rP_SEZ_rel; % Position in SEZ
    
    

% --- Velocity Transformation ---
% v_ECI = v_rel_ECI + (omega_E x r_ECI)
% v_rel_ECI = v_ECI - (omega_E x r_ECI)
v_rel_ECI = v0 - cross(omega_E_vec, r0);
v_rel_ECEF = ECEF_R_ECI * v_rel_ECI;
v1 = SEZ_R_ECEF * v_rel_ECEF; % Velocity in SEZ
end