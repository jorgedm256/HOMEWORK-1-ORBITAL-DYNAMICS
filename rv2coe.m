function [a, e_mag, i, Omega, omega_arg, nu] = rv2coe(r_vec, v_vec, mu)
    % Extracts Classical Orbital Elements from state vectors
    
    r = norm(r_vec);
    v = norm(v_vec);
    
    h_vec = cross(r_vec, v_vec);
    h = norm(h_vec);
    
    i = acos(h_vec(3) / h);
    %%%% i is the angle that the K vector and h make
    
    K_vec = [0; 0; 1];
    n_vec = cross(K_vec, h_vec);
    n = norm(n_vec);
    
    if n ~= 0
        Omega = acos(n_vec(1) / n);
        if n_vec(2) < 0
            Omega = 2*pi - Omega;
        end
    else
        Omega = 0;
    end
    %%%% Omega (RAAN) is the angle between the I vector (vernal equinox in
    %%%% ECI) and n (line of nodes)


    e_vec = (1/mu) * ((v^2 - mu/r)*r_vec - dot(r_vec, v_vec)*v_vec);
    e_mag = norm(e_vec);

    %%%% e is got from the expression of slide 35 of Unit 1 presentation
    %%%% but simplified so that we do not have any cross product
    
    if n ~= 0 && e_mag > 1e-7
        omega_arg = acos(dot(n_vec, e_vec) / (n * e_mag));
        if e_vec(3) < 0
            omega_arg = 2*pi - omega_arg; % Because acos only gives values up to $\pi$, we must subtract the result from $2\pi$ to get the correct angle in the second half of the orbit
        end
    else
        omega_arg = 0;
    end
    
    if e_mag > 1e-7
        nu = acos(dot(e_vec, r_vec) / (e_mag * r));
        if dot(r_vec, v_vec) < 0
             nu = 2*pi - nu;  % Because acos only gives values up to $\pi$, we must subtract the result from $2\pi$ to get the correct angle in the second half of the orbit
        end
    else
        nu = 0;
    end
    % We use Vis-Viva equation
    xi = (v^2)/2 - mu/r;% Specific mechanical energy
    a = -mu / (2 * xi);
end
