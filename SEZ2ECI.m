function [r0,v0] = SEZ2ECI(r1,v1,rS,phi,t)
%%%%% In r and v, 0 stands for S0 (ECI), 1 stands for S1(ECEF) and 2 stands for S2 (SEZ)
%%% Subindexes in r and v stand for ground station (S) or spacecraft (P)

%%% t to be included as input in hours
earth_rot_angle = 2*pi*t/24;

%%%% The statement does not say anything, so 0º longitude is assumed
station_lon = 0; 

%%% phi stands for latitude, and should be included as input in radians
station_lat = phi;

station_radius = 6371 + rS;

omega_earth = 2*pi/(24*3600)*[0;0;1]; %%% Earth angular velocity. 10 means that we change from S0 to S1

% Earth rotation matrices
matrix_ECEF_to_ECI = [cos(earth_rot_angle), -sin(earth_rot_angle), 0; 
                      sin(earth_rot_angle),  cos(earth_rot_angle), 0;
                      0 ,                    0 ,                   1];
                      
matrix_ECI_to_ECEF = matrix_ECEF_to_ECI'; %% Transpose

%Rotation from SEZ to ECEF
rot_Z_lon = [cos(station_lon), -sin(station_lon), 0;
             sin(station_lon),  cos(station_lon), 0;
             0,                 0,                1];  %%% Rotation about Z-axis of ECEF
             
rot_Y_lat = [sin(station_lat),   0,   cos(station_lat);
             0,                  1,   0    ;
            -cos(station_lat),   0,   sin(station_lat)];  %%%% Rotation about Y-axis of intermediate ref frame
            
matrix_SEZ_to_ECEF = rot_Z_lon * rot_Y_lat;
matrix_ECEF_to_SEZ = matrix_SEZ_to_ECEF';

%%%% Ground station position
pos_station_ECEF = station_radius * [cos(station_lat)*cos(station_lon); cos(station_lat)*sin(station_lon); sin(station_lat)];
pos_station_ECI = matrix_ECEF_to_ECI * pos_station_ECEF;

% Spacecraft relative coordinates
spacecraft_pos_SEZ = r1;
spacecraft_vel_SEZ = v1;

spacecraft_pos_rel_ECI = matrix_ECEF_to_ECI * matrix_SEZ_to_ECEF * spacecraft_pos_SEZ;
spacecraft_vel_rel_ECI = matrix_ECEF_to_ECI * matrix_SEZ_to_ECEF * spacecraft_vel_SEZ;

vel_P_1_ECI = spacecraft_vel_rel_ECI; %%%% because basis 1 and 2 share the same origin
%%%%% there is no angular velocity between ECEF and SEZ

pos_P_1_ECI = pos_station_ECI + spacecraft_pos_rel_ECI;

pos_P_0_ECI = pos_P_1_ECI;
vel_P_0_ECI = vel_P_1_ECI + cross(omega_earth, pos_P_1_ECI);

r0 = pos_P_0_ECI;
v0 = vel_P_0_ECI;
end