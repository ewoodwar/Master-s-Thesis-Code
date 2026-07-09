rho1= 22;
rho2= 26;

alpha1= 10;
alpha2= 10;

% Grid
mu1_vals= linspace(1.5,2,150);
mu2_vals= linspace(1.5,2,150);
[MU1,MU2]= meshgrid(mu1_vals,mu2_vals);

NU1= exp(-MU1).*(1 + alpha1);
NU2= exp(-MU2).*(1 + alpha2);

ONE_MINUS_NU1= 1 - NU1;
ONE_MINUS_NU2= 1 - NU2;

D1= (exp(-MU1).*alpha1.^2 ./ (MU1*rho1)) .* (1 - exp(-MU1));
D2= (exp(-MU2).*alpha2.^2 ./ (MU2*rho2)) .* (1 - exp(-MU2));

% w1 and w2
w1 = (1 ./ ONE_MINUS_NU1) .* ((ONE_MINUS_NU2 ./ D2) - (ONE_MINUS_NU1 ./ D1));
w2 = (1 ./ ONE_MINUS_NU2) .* ((ONE_MINUS_NU1 ./ D1) - (ONE_MINUS_NU2 ./ D2));

% Mask invalid region
mask = (NU1 > 1) & (NU2 > 1);
w1(~mask) = NaN;
w2(~mask) = NaN;

% Colour scaling 
vals = abs([w1(:); w2(:)]);
vals = vals(~isnan(vals));

L = prctile(vals, 95);   % ignore extreme spikes

% Red–White–Blue colormap
n = 256;
half = n/2;

r = [linspace(1,1,half) linspace(1,0,half)]';
g = [linspace(0,1,half) linspace(1,0,half)]';
b = [linspace(0,1,half) linspace(1,1,half)]';

redblue = [r g b];

% Plot
figure;

%  w1
subplot(1,2,1)
contourf(MU1, MU2, w1, 30, 'LineStyle','none');
colormap(redblue)
caxis([-L L])
colorbar


hold on
contour(MU1, MU2, w1, [0 0], 'k', 'LineWidth', 1.5) % zero line
plot([min(mu1_vals) max(mu1_vals)], [min(mu2_vals) max(mu2_vals)], 'k--', 'LineWidth', 1);

title('$w_1$','Interpreter','latex')
xlabel('$\mu_{1,1}$','Interpreter','latex')
ylabel('$\mu_{2,1}$','Interpreter','latex')

axis equal
set(gca,'FontSize',14,'TickLabelInterpreter','latex')

% w2
subplot(1,2,2)
contourf(MU1, MU2, w2, 30, 'LineStyle','none');
colormap(redblue)
caxis([-L L])
colorbar


hold on
contour(MU1, MU2, w2, [0 0], 'k', 'LineWidth', 1.5) % zero line
plot([min(mu1_vals) max(mu1_vals)], [min(mu2_vals) max(mu2_vals)], 'k--', 'LineWidth', 1);

title('$w_2$','Interpreter','latex')
xlabel('$\mu_1$','Interpreter','latex')
ylabel('$\mu_2$','Interpreter','latex')

axis equal
set(gca,'FontSize',14,'TickLabelInterpreter','latex')


%export
set(gcf,'Units','centimeters');
pos = get(gcf,'Position');

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[pos(3) pos(4)]);
set(gcf,'PaperPosition',[0 0 pos(3) pos(4)]);
set(gcf,'PaperPositionMode','manual');


print(gcf,'w1_w2_plots_1.pdf','-dpdf','-painters');