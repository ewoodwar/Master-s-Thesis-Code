function [y,z] = FBS_2p_TI(mu1,mu2,rho1,rho2,alpha1,alpha2,init_c_1,init_c_2,init_r_1,init_r_2, p1,p2,T)

t = linspace(0,T,T+1);

test = Inf;
tol = 1e-3;
iter = 0;
maxIter = 200;

u_1_I = 0.1*ones(1,T);
u_2_I = 0.1*ones(1,T);

u_1_T = 0.8*ones(1,T);
u_2_T = 0.8*ones(1,T);

v_1 = zeros(1,T+1);
v_2 = zeros(1,T+1);
w_1 = zeros(1,T+1);
w_2 = zeros(1,T+1);

v_1(1) = init_r_1;
v_2(1) = init_r_2;
w_1(1) = init_c_1;
w_2(1) = init_c_2;

lambda = zeros(4,T+1);


while test > tol && iter < maxIter

    iter = iter + 1;

    % Store old values
    oldu_1_I = u_1_I;
    oldu_2_I = u_2_I;
    oldu_1_T = u_1_T;
    oldu_2_T = u_2_T;

    oldv_1 = v_1;
    oldv_2 = v_2;
    oldw_1 = w_1;
    oldw_2 = w_2;

    oldlambda = lambda;

    %% ==========================
    %% FORWARD SWEEP
    %% ==========================

    for i = 1:T

        x = [v_1(i);
             v_2(i);
             w_1(i);
             w_2(i)];

        F = yearly_map( ...
            x,...
            u_1_T(i),u_2_T(i),...
            u_1_I(i),u_2_I(i),...
            mu1,mu2,rho1,rho2,...
            alpha1,alpha2,p1,p2);

        v_1(i+1) = F(1);
        v_2(i+1) = F(2);
        w_1(i+1) = F(3);
        w_2(i+1) = F(4);

    end

    %% ==========================
    %% BACKWARD SWEEP
    %% ==========================

    lambda(:,T+1) = 0;

    for j = T:-1:1

        grad_L = [0;
                  0;
                  2*w_1(j);
                  2*w_2(j)];

        Jac = diff_1_2p( ...
            v_1(j),v_2(j),...
            w_1(j),w_2(j),...
            u_1_I(j),u_2_I(j),...
            u_1_T(j),u_2_T(j),...
            mu1,mu2,rho1,rho2,...
            alpha1,alpha2,p1,p2);

        lambda(:,j) = grad_L + Jac'*lambda(:,j+1);

end
%fprintf('lambda range = [%g,%g]\n', ...
    %min(lambda(:)), max(lambda(:)));

%fprintf('Iter %d: max state = %.3e, max lambda = %.3e\n', ...
    %iter, max(abs([v_1 v_2 w_1 w_2])), max(abs(lambda(:))));
  %fprintf('Iter %d: max state = %.3e, max lambda = %.3e\n', ...
    %iter, max(abs([v_1 v_2 w_1 w_2])), max(abs(lambda(:))));

    %% ==========================
    %% CONTROL UPDATE
    %% ==========================

    etaI = 1e-2;
    etaT = 1e-3;
    maxG1 = 0;
    maxG2 = 0;
    maxG3 = 0;
    maxG4 = 0;

    for i = 1:T

        g1 = dif_u_p1I( ...
            v_1(i),v_2(i),w_1(i),w_2(i),...
            u_1_I(i),u_2_I(i),...
            u_1_T(i),u_2_T(i),...
            mu1,mu2,rho1,rho2,...
            alpha1,alpha2,...
            lambda(1,i+1),lambda(2,i+1),...
            lambda(3,i+1),lambda(4,i+1),...
            p1,p2);

        g2 = dif_u_p2I( ...
            v_1(i),v_2(i),w_1(i),w_2(i),...
            u_1_I(i),u_2_I(i),...
            u_1_T(i),u_2_T(i),...
            mu1,mu2,rho1,rho2,...
            alpha1,alpha2,...
            lambda(1,i+1),lambda(2,i+1),...
            lambda(3,i+1),lambda(4,i+1),...
            p1,p2);

        g3 = dif_u_p1T( ...
            v_1(i),v_2(i),w_1(i),w_2(i),...
            u_1_I(i),u_2_I(i),...
            u_1_T(i),u_2_T(i),...
            mu1,mu2,rho1,rho2,...
            alpha1,alpha2,...
            lambda(1,i+1),lambda(2,i+1),...
            lambda(3,i+1),lambda(4,i+1),...
            p1,p2);

        g4 = dif_u_p2T( ...
            v_1(i),v_2(i),w_1(i),w_2(i),...
            u_1_I(i),u_2_I(i),...
            u_1_T(i),u_2_T(i),...
            mu1,mu2,rho1,rho2,...
            alpha1,alpha2,...
            lambda(1,i+1),lambda(2,i+1),...
            lambda(3,i+1),lambda(4,i+1),...
            p1,p2);
    maxG1 = max(maxG1,abs(g1));
    maxG2 = max(maxG2,abs(g2));
    maxG3 = max(maxG3,abs(g3));
    maxG4 = max(maxG4,abs(g4));

    theta = 0.9;   % averaging parameter

    u1I_new = u_1_I(i) - etaI*g1;
    u2I_new = u_2_I(i) - etaI*g2;

    u1T_new = u_1_T(i) - etaT*g3;
    u2T_new = u_2_T(i) - etaT*g4;

    u_1_I(i) = theta*u_1_I(i) + (1-theta)*u1I_new;
    u_2_I(i) = theta*u_2_I(i) + (1-theta)*u2I_new;

    u_1_T(i) = theta*u_1_T(i) + (1-theta)*u1T_new;
    u_2_T(i) = theta*u_2_T(i) + (1-theta)*u2T_new;

        u_1_I(i) = min(0.9,max(0,u_1_I(i)));
        u_2_I(i) = min(0.9,max(0,u_2_I(i)));

        u_1_T(i) = min(0.9,max(0,u_1_T(i)));
        u_2_T(i) = min(0.9,max(0,u_2_T(i)));




end
%fprintf('G1 %.3e  G2 %.3e  G3 %.3e  G4 %.3e\n',maxG1,maxG2,maxG3,maxG4);
%fprintf('maxGrad = %.3e\n',maxGrad);
%fprintf('u1I [%g,%g]  u2I [%g,%g]\n',min(u_1_I),max(u_1_I),min(u_2_I),max(u_2_I));
%fprintf('u1T [%g,%g]  u2T [%g,%g]\n',min(u_1_T),max(u_1_T),min(u_2_T),max(u_2_T));

    %% ==========================
    %% CONVERGENCE TEST

    err1 = max(abs(u_1_I - oldu_1_I));
    err2 = max(abs(u_2_I - oldu_2_I));

    err3 = max(abs(u_1_T - oldu_1_T));
    err4 = max(abs(u_2_T - oldu_2_T));

    err5 = max(abs(v_1 - oldv_1));
    err6 = max(abs(v_2 - oldv_2));

    err7 = max(abs(w_1 - oldw_1));
    err8 = max(abs(w_2 - oldw_2));

    err9 = max(max(abs(lambda-oldlambda)));
%fprintf('Controls = %.3e\n', max([err1 err2 err3 err4]));
%fprintf('States   = %.3e\n', max([err5 err6 err7 err8]));
%fprintf('Lambda   = %.3e\n', err9);

    test = max([err1 err2 err3 err4 err5 err6 err7 err8 err9]);

    fprintf('Iter %d: max error = %.8e\n',iter,test);
    %fprintf(['Iter %d | ' 'u1=%.2e u2=%.2e t1=%.2e t2=%.2e ' 'v1=%.2e v2=%.2e w1=%.2e w2=%.2e lam=%.2e\n'], iter,err1,err2,err3,err4, err5,err6,err7,err8,err9);

end



if iter >= maxIter
    fprintf('WARNING: Maximum iterations reached.\n');
else
    fprintf('Converged in %d iterations.\n',iter);
end

y(1,:) = t;
y(2,:) = v_1;
y(3,:) = v_2;
y(4,:) = w_1;
y(5,:) = w_2;

z(1,:) = t(1:T);
z(2,:) = u_1_I;
z(3,:) = u_2_I;
z(4,:) = u_1_T;
z(5,:) = u_2_T;
