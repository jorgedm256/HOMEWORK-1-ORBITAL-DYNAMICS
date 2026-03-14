function [r_ECI, v_ECI] = coe2rv(mu, p, e, i, Omega, omega, nu)
% Computes r and v in ECI given classical orbital elements
    
% 1. State in perifocal frame (PQW)
r_pqw = (p / (1 + e*cos(nu))) * [cos(nu); sin(nu); 0];
v_pqw = sqrt(mu/p) * [-sin(nu); e + cos(nu); 0];
    
% 2. Rotation matrix from PQW to ECI (3-1-3 sequence)
cw = cos(omega); sw = sin(omega);
cO = cos(Omega); sO = sin(Omega);
ci = cos(i);     si = sin(i);
    
Q_pX = [cO*cw - sO*ci*sw,  -cO*sw - sO*ci*cw,  sO*si;
        sO*cw + cO*ci*sw,  -sO*sw + cO*ci*cw, -cO*si;
        si*sw,              si*cw,             ci];
            
% 3. Transform to ECI
r_ECI = Q_pX * r_pqw;
v_ECI = Q_pX * v_pqw;
end