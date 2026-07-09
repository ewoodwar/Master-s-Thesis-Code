function pdot = pachepsky2Pf_mu(t, p, mu_1, mu_2)
% Parameters 
rho_1 = 22;
rho_2 = 25;
alpha_1 = 10;
alpha_2 = 10;

% Patch 1 dynamics
pdot(1,:) = rho_1 * p(1) * (1 - p(1)) - alpha_1 * p(1) * p(2); % Resource 1
pdot(2,:) = -mu_1 * p(2);                                       % Consumer 1
pdot(3,:) = alpha_1 * p(1);                                      % Biomass 1

% Patch 2 dynamics
pdot(4,:) = rho_2 * p(4) * (1 - p(4)) - alpha_2 * p(4) * p(5);  % Resource 2
pdot(5,:) = -mu_2 * p(5);                                       % Consumer 2
pdot(6,:) = alpha_2 * p(4);                                      % Biomass 2

end