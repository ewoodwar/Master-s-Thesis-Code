% Params
rho1= 25;
rho2= 22;

alpha1= 10;
alpha2= 10;

p_small= 0.99;
p_none= 1.0;

years= 800;
burnin= 500;


% asymptotics
mu1_vals= linspace(1.5,2,50);
mu2_vals= linspace(1.5,2,50);
[MU1,MU2]= meshgrid(mu1_vals,mu2_vals);

NU1= exp(-MU1).*(1 + alpha1);
NU2= exp(-MU2).*(1 + alpha2);

ONE_MINUS_NU1= 1 - NU1;
ONE_MINUS_NU2= 1 - NU2;
D1= (exp(-MU1).*alpha1.^2 ./ (MU1*rho1)) .* (1 - exp(-MU1));
D2= (exp(-MU2).*alpha2.^2 ./ (MU2*rho2)) .* (1 - exp(-MU2));

W= (1./ONE_MINUS_NU1 - 1./ONE_MINUS_NU2) .*(ONE_MINUS_NU2./D2 - ONE_MINUS_NU1./D1);

mask= (NU1 > 1) & (NU2 > 1);
W(~mask)= NaN;


% plotting asymptotics
hfig1= figure;
contourf(MU1,MU2,sign(W),[-1 0 1],'LineStyle','none')
hold on
colormap([1 0.8 0.8; 0.8 0.9 1])

contour(MU1,MU2,W,[0 0],'k','LineWidth',2)
xlabel('$\mu_1$','Interpreter','latex')
ylabel('$\mu_2$','Interpreter','latex')
axis equal
set(gca,'FontSize',14,'TickLabelInterpreter','latex')


% simulation agreement
mu1_sim= linspace(1.5,2,10);
mu2_sim= linspace(1.5,2,10);

n1= length(mu1_sim);
n2= length(mu2_sim);

sim_benefit= NaN(n1,n2);
theory_benefit= NaN(n1,n2);


% both
for i= 1:n1
    for j= 1:n2
        mu1= mu1_sim(i);
        mu2= mu2_sim(j);

        nu1= exp(-mu1)*(1+alpha1);
        nu2= exp(-mu2)*(1+alpha2);

        if (nu1>1)&&(nu2>1)
            D1 = (exp(-mu1)*alpha1^2/(mu1*rho1))*(1-exp(-mu1));
            D2 = (exp(-mu2)*alpha2^2/(mu2*rho2))*(1-exp(-mu2));

            Wpoint = (1/(1-nu1)-1/(1-nu2))*((1-nu2)/D2-(1-nu1)/D1);
            theory_benefit(i,j) = Wpoint > 0;
        end
    end
end


% Sims w run_model_mu
for i= 1:n1
    for j= 1:n2
        if isnan(theory_benefit(i,j))
            continue
        end

        mu1= mu1_sim(i);
        mu2= mu2_sim(j);

        N_none= run_model_mu_consumer(p_none, p_none, years, burnin, mu1, mu2);
        N_small= run_model_mu_consumer(p_small, p_small, years, burnin, mu1, mu2);

        % Beneficial dispersal?
        sim_benefit(i,j)= (N_small - N_none) > 0;
    end
end


% OVERLAY: dots = beneficial, X = detrimental
for i = 1:n1
    for j = 1:n2
        if isnan(sim_benefit(i,j))
            continue
        end

        mu1 = mu1_sim(i);
        mu2 = mu2_sim(j);

        if sim_benefit(i,j)
            plot(mu1, mu2, 'k.', 'MarkerSize', 18)  % beneficial
        else
            plot(mu1, mu2, 'kx', 'MarkerSize', 12, 'LineWidth', 1.5)  % detrimental
        end
    end
end
% =========================

% EXPORT
set(hfig1,'Units','centimeters');
pos = get(hfig1,'Position');

set(hfig1,'PaperUnits','centimeters');
set(hfig1,'PaperSize',[pos(3) pos(4)]);
set(hfig1,'PaperPosition',[0 0 pos(3) pos(4)]);
set(hfig1,'PaperPositionMode','manual');

print(hfig1,'theory_vs_simulation_overlay_mu.png','-dpng','-r300');