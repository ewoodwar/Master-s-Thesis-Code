function Rmean = run_model_mu_resource(p1, p2, years, burnin, mu1, mu2)

% Initial conditions [R1, C1, B1, R2, C2, B2]
p0 = [0.5, 1, 0, 0.5, 1, 0];

totalR = zeros(1, years);

for i = 1:years
    tspan = [0 1];
    
    % ODE step
    [~, p] = ode23s(@(t,p) pachepsky2Pf_mu(t, p, mu1, mu2), tspan, p0);
    
    p_end = p(end,:);
    
    % Dispersal step (consumers move, resource does not)
    C1 = (p_end(3)+1)*p1*p_end(2) + (p_end(6)+1)*(1-p2)*p_end(5);
    C2 = (p_end(6)+1)*p2*p_end(5) + (p_end(3)+1)*(1-p1)*p_end(2);
    
    p_end(2) = C1;
    p_end(5) = C2;
    
    % Reset biomass
    p_end(3) = 0;
    p_end(6) = 0;
    
    % Store TOTAL RESOURCE
    totalR(i) = p_end(1) + p_end(4);
    
    % Update state
    p0 = p_end;
end

% Return long-term average resource
Rmean = mean(totalR(burnin:end));

end