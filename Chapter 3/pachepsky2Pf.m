function pdot=pachepsky2Pf(t,p)
%This is the ODE for the scaled pachepsky model on two patches
rho_1=20;
rho_2=20;
alpha_1=10;
alpha_2=10;
mu_1=2;
mu_2=2;

%Patch 1
pdot(1,:)=rho_1.*p(1).*(1-p(1))-alpha_1.*p(1).*p(2); %Resource Patch 1
pdot(2,:)=-mu_1.*p(2); %Consumer Patch 1
pdot(3,:)=alpha_1.*p(1); %Biomass Patch 1

%Patch 2
pdot(4,:)=rho_2.*p(4).*(1-p(4))-alpha_2.*p(4).*p(5); %Resource Patch 2
pdot(5,:)=-mu_2.*p(5); %Consumer Patch 2
pdot(6,:)=alpha_2.*p(4);% Biomass Patch 2