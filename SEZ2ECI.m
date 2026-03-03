function [r0,v0] = SEZ2ECI(r1,v1,rS,phi,t)
%%%%% In r and v, 0 stands for S0 (ECI), 1 stands for S1(ECEF) and 2 stands for S2 (SEZ)

%%% Subindexes in r and v stand for ground station (S) or spacecraft (P)

omega_10_ECI=2*pi/(24*3600)*[0;0;1]; %%% Earth angular velocity. 10 means that we change from S0 to S1

%%% t to be included as input in hours
psi=2*pi*t/24;

ECI_R_ECEF=[cos(psi), -sin(psi), 0; 
            sin(psi),  cos(psi), 0;
            0 ,         0 ,      1];
ECEF_R_ECI=ECI_R_ECEF';  %% Transpose

R=6371+rS;

%%% phi stands for latitude, and should be included as input in degrees
lat=phi*pi/180;
%%%% The statement does not say anything, so 0º longitude is assumed
lon= 0; 

%%%% Ground station position
rS_1_ECEF=R*[cos(lat)*cos(lon); cos(lat)*sin(lon); sin(lat)];
rS_0_ECI=ECI_R_ECEF*rS_1_ECEF;

%Rotation from SEZ to ECEF

Rlon=[cos(lon), -sin(lon), 0;
      sin(lon),   cos(lon),  0;
      0,          0,       1];  %%% About z-axis of ECEF

Rlat=[sin(lat),   0,   cos(lat);
       0,        1,        0    ;
       -cos(lat),  0,      sin(lat)];  %%%% About y-axis of intermediate ref frame

ECEF_R_SEZ=Rlon*Rlat;
SEZ_R_ECEF=ECEF_R_SEZ';

rP_2_SEZ=r1;
rP_2_ECI=ECI_R_ECEF*ECEF_R_SEZ*rP_2_SEZ;


vP_2_SEZ =v1;

vP_2_ECI=ECI_R_ECEF*ECEF_R_SEZ*vP_2_SEZ;
vP_1_ECI=vP_2_ECI;  %%%% because basis 1 and 2 share the same origin

%%%%% there is no angular velocity between ECEF and SEZ

rP_1_ECI=rS_0_ECI+rP_2_ECI;
rP_0_ECI=rP_1_ECI;
vP_0_ECI= vP_1_ECI+cross(omega_10_ECI,rP_1_ECI);

r0=rP_0_ECI;
v0=vP_0_ECI;
end