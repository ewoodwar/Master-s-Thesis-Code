function [v_map_F_1,v_map_F_2, w_map_C_1,w_map_C_2] = map_twopatch(v_1_0,v_2_0,w_1_0,w_2_0,u_1,u_2,mu1,mu2, rho1,rho2, alpha1,alpha2,timing_1,timing_2,p1,p2)

%Map One takes initial state to state before lambda
nodes_1=0:0.1:timing_1; %select node spacing
nodes_2=0:0.1:timing_2; %select node spacing
f_vals_1 = exp(rho1*nodes_1-(((alpha1*w_1_0)/mu1)*(1-exp(-mu1*nodes_1))));
f_vals_2 = exp(rho2*nodes_2-(((alpha2*w_2_0)/mu2)*(1-exp(-mu2*nodes_2))));
INT_f_1=trapz(nodes_1, f_vals_1);
INT_f_2=trapz(nodes_2, f_vals_2);
F_1=v_1_0.*exp(rho1*timing_1-(((alpha1.*w_1_0)./mu1).*(1-exp(-mu1*timing_1))))./((rho1.*v_1_0.*INT_f_1)+1);
F_2=v_2_0.*exp(rho2*timing_2-(((alpha2.*w_2_0)./mu2).*(1-exp(-mu2*timing_2))))./((rho2.*v_2_0.*INT_f_2)+1);
%outer_int=int_eval_1(rho, alpha, mu, w_0, v_0, nodes_1); %output of this function is the solved integral
C_1_1=w_1_0.*exp(-mu1*timing_1);
C_2_1=w_2_0.*exp(-mu2*timing_2);
outer_int_1_1=int_eval_2p(rho1, alpha1, mu1, w_1_0, v_1_0, nodes_1);
outer_int_2_1=int_eval_2p(rho2, alpha2, mu2, w_2_0, v_2_0, nodes_2);
B_1_1= alpha1.*v_1_0.*outer_int_1_1;
B_2_1= alpha2.*v_2_0.*outer_int_2_1;


%second map takes state before spraying to immediately after spraying
  v_1_2=F_1;
  v_2_2=F_2;
  w_1_2=(1-u_1)*C_1_1;
  w_2_2=(1-u_2)*C_2_1;

%third map takes after spraying to just before the end of the year
nodes_1_2=0:0.1:(1-timing_1); %select node spacing for the time between lambda and the end of the year
nodes_2_2=0:0.1:(1-timing_2);
f_vals_1_2 = exp(rho1*nodes_1_2-(((alpha1*w_1_2)/mu1)*(1-exp(-mu1*nodes_1_2))));
f_vals_2_2 = exp(rho2*nodes_2_2-(((alpha2*w_2_2)/mu2)*(1-exp(-mu2*nodes_2_2))));
INT_f_1_2=trapz(nodes_1_2, f_vals_1_2);
INT_f_2_2=trapz(nodes_2_2, f_vals_2_2);
F_1_2=(v_1_2.*exp(rho1*(1-timing_1)-(((alpha1.*w_1_2)./mu1).*(1-exp(-mu1*(1-timing_1))))))./((rho1.*v_1_2.*INT_f_1_2)+1);
F_2_2=(v_2_2.*exp(rho2*(1-timing_2)-(((alpha2.*w_2_2)./mu2).*(1-exp(-mu2*(1-timing_2))))))./((rho2.*v_2_2.*INT_f_2_2)+1);

C_1_2=w_1_2.*exp(-mu1*(1-timing_1));
C_2_2=w_2_2.*exp(-mu2*(1-timing_2));
outer_int_1_2=int_eval_2p(rho1, alpha1, mu1, w_1_2, v_1_2, nodes_1_2);
outer_int_2_2=int_eval_2p(rho2, alpha2, mu2, w_2_2, v_2_2, nodes_2_2);
B_1_2= alpha1.*v_1_2.*outer_int_1_2;
B_2_2= alpha2.*v_2_2.*outer_int_2_2;
b1=B_1_2+B_1_1;
b2=B_2_2+B_2_1;


  %fourth map takes state before the end of the year to beginning of the next year (after birth pulse)
  v_1_4=F_1_2;
  v_2_4=F_2_2;
  w_1_4=p1*(1+b1)*C_1_2+(1-p2)*(1+b2)*C_2_2;
  w_2_4=p2*(1+b2)*C_2_2+(1-p1)*(1+b1)*C_1_2;

v_map_F_1=v_1_4;
v_map_F_2=v_2_4;
w_map_C_1=w_1_4;
w_map_C_2=w_2_4;
