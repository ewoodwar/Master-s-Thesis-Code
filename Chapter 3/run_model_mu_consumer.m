function Nmean = run_model_mu(p1, p2, years, burnin, mu1, mu2)
% Initial conditions [R1, C1, B1, R2, C2, B2]
p0 = [0.5, 1, 0, 0.5, 1, 0];

totalC = zeros(1, years);  % preallocate

for i = 1:years
    tspan = [0 1];
    
    % Call ODE with mu1 and mu2
    [~, p] = ode23s(@(t,p) pachepsky2Pf_mu(t, p, mu1, mu2), tspan, p0);
    
    p_end = p(end,:);
    
    % Yearly dispersal step
    C1 = (p_end(3)+1)*p1*p_end(2) + (p_end(6)+1)*(1-p2)*p_end(5);
    C2 = (p_end(6)+1)*p2*p_end(5) + (p_end(3)+1)*(1-p1)*p_end(2);
    
    p_end(2) = C1;
    p_end(5) = C2;
    
    % Reset biomass for next year
    p_end(3) = 0;
    p_end(6) = 0;
    
    % Store total consumers
    totalC(i) = p_end(2) + p_end(5);
 
    
    % Update initial conditions for next year
    p0 = p_end;
end

% Remove burn-in transient and compute mean
Nmean = mean(totalC(burnin:end));

end