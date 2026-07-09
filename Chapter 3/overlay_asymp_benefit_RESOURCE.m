clear; clc;

% Parameters
rho1 = 22;  rho2 = 25;
alpha1 = 10;  alpha2 = 10;

p_small = 0.99;
p_none = 1.0;
years  = 500;
burnin = 300;

% Grid
mu1_vals = linspace(1.501, 1.8, 50);
mu2_vals = linspace(1.501, 1.8, 50);
[MU1, MU2] = meshgrid(mu1_vals, mu2_vals);

% Numerical integration setup
tau = linspace(0, 1, 100);
sigma_fun = @(tau,rho,alpha,mu,w) exp(rho.*tau - (alpha*w./mu).*(1 - exp(-mu.*tau)));

% Precompute constants
NU1 = exp(-MU1).*(1 + alpha1);
NU2 = exp(-MU2).*(1 + alpha2);

ONE_MINUS_NU1 = 1 - NU1;
ONE_MINUS_NU2 = 1 - NU2;

B1 = max((exp(-MU1).*alpha1.^2 ./ (MU1.*rho1)) .* (1 - exp(-MU1)), 1e-12);
B2 = max((exp(-MU2).*alpha2.^2 ./ (MU2.*rho2)) .* (1 - exp(-MU2)), 1e-12);

% Compute baseline equilibrium w_{i,0}
W1_0 = -ONE_MINUS_NU1 ./ B1;
W2_0 = -ONE_MINUS_NU2 ./ B2;

%W1_0 = max(W1_0, 1e-12); % ensure positive
%W2_0 = max(W2_0, 1e-12);

% Compute G_i' at w_{i,0}
G1 = NaN(size(MU1));
G2 = NaN(size(MU2));

for i = 1:size(MU1,1)
    for j = 1:size(MU1,2)
        mu1 = MU1(i,j);
        mu2 = MU2(i,j);
        w1_0 = W1_0(i,j);
        w2_0 = W2_0(i,j);

        sigma1 = sigma_fun(tau, rho1, alpha1, mu1, w1_0);
        beta1 = trapz(tau, sigma1);
        I1 = trapz(tau, (1 - exp(-mu1*tau)) .* sigma1);
        sigma1_1 = sigma_fun(1, rho1, alpha1, mu1, w1_0);
        G1(i,j) = (-(alpha1/mu1)*(1 - exp(-mu1)) * sigma1_1 * beta1 + (alpha1/mu1)*(sigma1_1 - 1)*I1)/(rho1*beta1^2);

        sigma2 = sigma_fun(tau, rho2, alpha2, mu2, w2_0);
        beta2 = trapz(tau, sigma2);
        I2 = trapz(tau, (1 - exp(-mu2*tau)) .* sigma2);
        sigma2_1 = sigma_fun(1, rho2, alpha2, mu2, w2_0);
        G2(i,j) = (-(alpha2/mu2)*(1 - exp(-mu2)) * sigma2_1 * beta2 + (alpha2/mu2)*(sigma2_1 - 1)*I2)/(rho2*beta2^2);

    end
end

% Compute first-order dispersal effects w_{i,1}
w1_1 = (1 ./ ONE_MINUS_NU1) .* ((ONE_MINUS_NU2 ./ B2) - (ONE_MINUS_NU1 ./ B1));
w2_1 = (1 ./ ONE_MINUS_NU2) .* ((ONE_MINUS_NU1 ./ B1) - (ONE_MINUS_NU2 ./ B2));

% Compute true theoretical prediction V
V = G1 .* w1_1 + G2 .* w2_1;
V(~isfinite(V)) = 0; % remove any NaN/Inf

% Prepare plot
% Detect if all negative
if max(V(:)) <= 0
    % All negative → force all red
    V_plot = -ones(size(V));
elseif min(V(:)) >= 0
    % All positive → force all blue
    V_plot = ones(size(V));
else
    % Mixed values → use threshold
    tol = 1e-8;
    V_plot = zeros(size(V));
    V_plot(V > tol) = 1;
    V_plot(V < -tol) = -1;
end


hfig1 = figure;
contourf(MU1, MU2, V_plot, [-1 0 1], 'LineStyle', 'none');
hold on;
colormap([1 0.8 0.8; 0.8 0.9 1]); % red first, blue second
contour(MU1, MU2, V, [0 0], 'k', 'LineWidth', 2);

xlabel('$\mu_1$', 'Interpreter', 'latex');
ylabel('$\mu_2$', 'Interpreter', 'latex');
axis equal;
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex');

% Simulation overlay
mu1_sim = linspace(1.501, 1.8, 6);
mu2_sim = linspace(1.501, 1.8, 6);
sim_benefit = NaN(length(mu1_sim), length(mu2_sim));

for i = 1:length(mu1_sim)
    for j = 1:length(mu2_sim)
        mu1 = mu1_sim(i);
        mu2 = mu2_sim(j);
        N_none  = run_model_mu_resource(p_none, p_none, years, burnin, mu1, mu2);
        N_small = run_model_mu_resource(p_small, p_small, years, burnin, mu1, mu2);
        sim_benefit(i,j) = (N_small - N_none) > 0;
    end
end

for i = 1:length(mu1_sim)
    for j = 1:length(mu2_sim)
        if sim_benefit(i,j)
            plot(mu1_sim(i), mu2_sim(j), 'k.', 'MarkerSize', 18);
        else
            plot(mu1_sim(i), mu2_sim(j), 'kx', 'MarkerSize', 12, 'LineWidth', 1.5);
        end
    end
end

% Export figure
set(hfig1,'Units','centimeters');
pos = get(hfig1,'Position');
set(hfig1,'PaperUnits','centimeters');
set(hfig1,'PaperSize',[pos(3) pos(4)]);
set(hfig1,'PaperPosition',[0 0 pos(3) pos(4)]);
set(hfig1,'PaperPositionMode','manual');
print(hfig1,'overlay_resource_1.png','-dpng','-r300');
