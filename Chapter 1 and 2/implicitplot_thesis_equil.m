%set up initial grid
x = linspace(0, 5, 2000);
y = linspace(0, 3, 2000);
[X,Y] = meshgrid(x,y);
hfig1=figure;

%loops over cases
for k = 1:3
    subplot(1,3,k)
    hold on
    box on

    %set parameters for case k
    switch k
        case 1
            p1 = 0.9;
            p2 = 0.21;
            mu1 = 0.5;
            mu2 = 0.4;
            alpha1 =10;
            alpha2 = 10;
            rho1 = 20;
            rho2 = 20;
            nu1=exp(-mu1)*(alpha1+1);
            nu2=exp(-mu2)*(alpha2+1);

        case 2
           p1 = 0.9;
            p2 = 0.24;
            mu1 = 0.5;
            mu2 = 0.4;
            alpha1 =10;
            alpha2 = 10;
            rho1 = 20;
            rho2 = 20;
            nu1=exp(-mu1)*(alpha1+1);
            nu2=exp(-mu2)*(alpha2+1);

        case 3
           p1 = 0.9;
            p2 = 0.3;
            mu1 = 0.5;
            mu2 = 0.4;
            alpha1 =10;
            alpha2 = 10;
            rho1 = 20;
            rho2 = 20;
            nu1=exp(-mu1)*(alpha1+1);
            nu2=exp(-mu2)*(alpha2+1);

    end
%evaluate implicit functions
F = X-X.*p1.*( nu1 - exp(-mu1).*(alpha1.^2 .*X./(mu1*rho1)).*(1-exp(-mu1)) )-Y.*(1-p2).*(nu2-exp(-mu2).*(alpha2.^2.*Y./(mu2*rho2)).*(1-exp(-mu2)));
G = Y-Y.*p2.*( nu2 - exp(-mu2).*(alpha2.^2 .*Y./(mu2*rho2)).*(1-exp(-mu2)) )-X.*(1-p1).*(nu1-exp(-mu1).*(alpha1.^2.*X./(mu1*rho1)).*(1-exp(-mu1)));

%plot zero contours
contour(X,Y,F,[0 0],'color',[0.66, 0.33, 0.88],'LineWidth',2)
contour(X,Y,G,[0 0],'color',[0.5, 0.85, 1.0],'LineWidth',2)
xlabel('Consumer Patch 1')
ylabel('Consumer Patch 2')
axis equal
tol = 1e-2;
%Plotting the intersections to a tolerance of tol
idx = abs(F) < tol & abs(G) < tol;
plot(X(idx), Y(idx), 'ko', 'MarkerSize', 3, 'MarkerFaceColor','k')
end
drawnow
exportgraphics(gcf,'intersection_ellipse_4.pdf','ContentType','vector','BackgroundColor','none')
