function theta = E2theta (E,e)

theta=acos((cos(E)-e)/(1-e*cos(E)));
if mod(E,2*pi)> pi
    theta=2*pi-theta;
end