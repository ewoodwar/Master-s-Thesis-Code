function pdot=pachepsky2Pf(t,p)
%This is the ODE for the scaled pachepsky model on two patches
rho=20;
alpha=10;
mu=1.5;

%Patch 1
pdot(1,:)=rho.*p(1).*(1-p(1))-alpha.*p(1).*p(2); %Resource Patch 1
pdot(2,:)=-mu.*p(2); %Consumer Patch 1
pdot(3,:)=alpha.*p(1); %Biomass Patch 1

%Patch 2
pdot(4,:)=rho.*p(4).*(1-p(4))-alpha.*p(4).*p(5); %Resource Patch 2
pdot(5,:)=-mu.*p(5); %Consumer Patch 2
pdot(6,:)=alpha.*p(4);% Biomass Patch 2
