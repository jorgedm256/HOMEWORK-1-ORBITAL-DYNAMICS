function E =theta2E(theta,e)

E=acos((e+cos(theta))/(1+e*cos(theta)));

if mod(theta,2*pi)> pi
    E=2*pi-E;
end

