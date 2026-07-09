%This is a method to solve the Pachepsky equation on two patches. This plots the yearly map, using the 1+ (right after birth pulse), as the sampling time
clear all
%%params%%
rho=10;
alpha=15;
mu=0.1;

%% initial conditions
v=1;
w=0.75;
vv=1;
ww=0.27;
p11=0.1; %Probability of staying in patch 1
p22=0.67; %Probability of staying in patch 2
p12=1-p22; %Probability of moving from patch 2 to patch 1
p21=1-p11; %Probability of moving from patch 1 to patch 2
years=200; %number of years in the simulation
nodes=0:0.01:1; %select node spacing

%%% Calculate the integral parts of the resource equation %%%
for j= 1:years
  %Calculating Resources in both paches after reproduction
  f_vals_1 = exp(rho*nodes-(((alpha*w)/mu)*(1-exp(-mu*nodes)))); %Resource P1 integral part of equation
  f_vals_2 = exp(rho*nodes-(((alpha*ww)/mu)*(1-exp(-mu*nodes)))); %Resource P2 integral part of equation

  INT_f_1=trapz(nodes, f_vals_1);
  INT_f_2=trapz(nodes, f_vals_2);
  F_after_rep_1(j)=(v.*exp(rho-(((alpha.*w)./mu).*(1-exp(-mu)))))./((rho.*v.*INT_f_1)+1);
  F_after_rep_2(j)=(vv.*exp(rho-(((alpha.*ww)./mu).*(1-exp(-mu)))))./((rho.*vv.*INT_f_2)+1);

  %Calculating Consumers
  outer_int=int_eval_2P(rho, alpha, mu, w, v,ww,vv, nodes); %output of this function is the solved integral
  C_after_rep_1(j)= ((alpha.*v.*outer_int(1))+1)*p11*(w.*exp(-mu))+((alpha.*vv.*outer_int(2))+1)*p12*(ww.*exp(-mu)); %Calculate the Wt+1 (consumer density after rep)
  C_after_rep_2(j)= ((alpha.*vv.*outer_int(2))+1)*p22*(ww.*exp(-mu))+((alpha.*v.*outer_int(1))+1)*p21*(w.*exp(-mu));
  v=F_after_rep_1(j);
  vv=F_after_rep_2(j);
  w=C_after_rep_1(j);
  ww=C_after_rep_2(j);

endfor

%%% Plot everything %%%
k=1:years;
subplot(2,2,1)
plot(k, F_after_rep_1), xlabel('Years'), ylabel('Density'), title('Resource Patch 1')
subplot(2,2,2)
plot(k, C_after_rep_1), xlabel('Years'), ylabel('Density'), title('Consumer Patch 1')
subplot(2,2,3)
plot(k, F_after_rep_2), xlabel('Years'), ylabel('Density'), title('Resource Patch 2')
subplot(2,2,4)
plot(k, C_after_rep_2), xlabel('Years'), ylabel('Density'), title('Consumer Patch 2')

