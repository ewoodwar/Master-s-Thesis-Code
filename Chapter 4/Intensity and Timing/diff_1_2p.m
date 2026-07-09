function Jac = diff_1_2p(v1,v2,w1,w2,u1,u2,u3,u4,mu1,mu2,rho1,rho2,alpha1,alpha2,p1,p2)


eps = 1e-6;
Jac = zeros(4,4);
x = [v1;
     v2;
     w1;
     w2];

for k = 1:4

    dx = zeros(4,1);
    dx(k) = eps;

    x_plus  = x + dx;
    x_minus = x - dx;

    F_plus = yearly_map( ...
        x_plus,...
        u3,u4,...
        u1,u2,...
        mu1,mu2,...
        rho1,rho2,...
        alpha1,alpha2,...
        p1,p2);

    F_minus = yearly_map( ...
        x_minus,...
        u3,u4,...
        u1,u2,...
        mu1,mu2,...
        rho1,rho2,...
        alpha1,alpha2,...
        p1,p2);

    Jac(:,k) = (F_plus - F_minus)/(2*eps);

end

end










##function [Jac] = diff_1_2p(v1,v2,w1,w2,u1,u2,u3,u4,mu1,mu2,rho1,rho2,alpha1,alpha2,p1,p2)
% Calculates the Jacobian of the map f(x, u) using central differences.
% x = [v1; v2; w1; w2]
% Jac(i,j) = ∂f_i / ∂x_j
##
##    h = 1e-5;
##
##    % Initialize Jacobian matrix
##    Jac = zeros(4,4);
##
##    % Baseline state vector
##    x = [v1; v2; w1; w2];
##
##    for j = 1:4
##        x_forward = x;
##        x_backward = x;
##
##        % Perturb j-th state variable
##        x_forward(j) = x(j) + h;
##        x_backward(j) = x(j) - h;
##
##        % Unpack perturbed variables
##        v1_f = x_forward(1); v2_f = x_forward(2);
##        w1_f = x_forward(3); w2_f = x_forward(4);
##
##        v1_b = x_backward(1); v2_b = x_backward(2);
##        w1_b = x_backward(3); w2_b = x_backward(4);
##
##        % Evaluate dynamics at x + h and x - h
##        [v1_plus, v2_plus, w1_plus, w2_plus] = map_2P_time(v1_f, v2_f, w1_f, w2_f, u1, u2,u3,u4, mu1, mu2, rho1, rho2, alpha1, alpha2, p1, p2);
##        [v1_minus, v2_minus, w1_minus, w2_minus] = map_2P_time(v1_b, v2_b, w1_b, w2_b, u1, u2,u3,u4, mu1, mu2, rho1, rho2, alpha1,alpha2, p1, p2);
##
##    Jac(1,j) = (v1_plus - v1_minus) / (2*h);  % f1 = v1_next
##    Jac(2,j) = (v2_plus - v2_minus) / (2*h);  % f2 = v2_next
##    Jac(3,j) = (w1_plus - w1_minus) / (2*h);  % f3 = w1_next
##    Jac(4,j) = (w2_plus - w2_minus) / (2*h);  % f4 = w2_next
##
##    end
##end
