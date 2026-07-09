function [deriv_h_u1] =dif_u_patch_1(v1,v2,w1,w2,u1,u2,mu1,mu2,rho1,rho2,alpha1,alpha2,adjoint_v_1,adjoint_v_2,adjoint_w_1,adjoint_w_2,timing_1,timing_2,p1,p2);
%in this section, I approximate the derivative of each of my function wrt that state

%Here I am using the 2-point difference method
h=0.0001; %step size
u_right_1 = min(0.9, u1 + h);
u_left_1 = max(0, u1 - h);
[v_right_1,v_right_2,w_right_1,w_right_2]=map_twopatch(v1,v2,w1,w2,u_right_1,u2,mu1,mu2,rho1,rho2,alpha1,alpha2,timing_1,timing_2,p1,p2);
[v_left_1,v_left_2,w_left_1,w_left_2]=map_twopatch(v1,v2,w1,w2,u_left_1,u2,mu1,mu2,rho1,rho2,alpha1,alpha2,timing_1,timing_2,p1,p2);

%This is the Hamiltonian- right now it includes ~minimizing cost only~ *since I want min on patch 1 I keep u2 the same
h_right = 5*u_right_1^2 + 30*u2^2 + w1^2 + w2^2 ...
          + adjoint_v_1 * v_right_1 ...
          + adjoint_v_2 * v_right_2 ...
          + adjoint_w_1 * w_right_1 ...
          + adjoint_w_2 * w_right_2;


h_left = 5*u_left_1^2 + 30*u2^2 + w1^2 + w2^2 ...
          + adjoint_v_1 * v_left_1 ...
          + adjoint_v_2 * v_left_2 ...
          + adjoint_w_1 * w_left_1 ...
          + adjoint_w_2 * w_left_2;

%Difference using divided difference formula
deriv_h_u1=(h_right-h_left)/(2*h);

