function [y,z] = FBS_twopatch_JOINT(mu1,mu2,rho1,rho2,alpha1,alpha2,timing_1,timing_2, init_c_1, init_c_2,init_r_1,init_r_2,p1,p2,T)
  %this code is the forward backward sweep algorithm to approximate the solution to an optimal control problem
  %this code is modified from the code created by Suzanne Lenhart
  %this code represents the optimal control of a consumer resource population on one patch, where only the consumer is controlled

test = -1; % initial convergence flag (force at least one iteration)

delta = 0.1; %tolerance

t = linspace(0,T,T+1);

u_1 = zeros(1,T)+0.1; %Control variable patch 1
u_2 = zeros(1,T)+0.1; %Control variable patch 2
w_1 = zeros(1,T+1); %Setting initial conditions
w_2 = zeros(1,T+1); %Setting initial conditions
w_1(1) = init_c_1;
w_2(1) = init_c_2;
v_1 = zeros(1,T+1);
v_2 = zeros(1,T+1);
v_1(1) = init_r_1;
v_2(1) = init_r_2;


lambda = zeros(4, T+1);
lambda(:, T+1) = [0; 0; 1; 1];  % terminal condition


while(test < -1.0e-8)
% pick values
    oldu_1 = u_1;
    oldu_2 = u_2;
    oldw_1= w_1;
    oldw_2= w_2;
    oldv_1 = v_1;
    oldv_2 = v_2;
    oldlambda=lambda;

  %apply the map to get new state variables use the func "map_twopatch"
    for i = 1:T
      [v_1(i+1),v_2(i+1),w_1(i+1),w_2(i+1)] = map_twopatch(v_1(i),v_2(i),w_1(i),w_2(i),u_1(i),u_2(i),mu1,mu2,rho1,rho2,alpha1,alpha2,timing_1,timing_2,p1,p2);

  end
  %lambda(:, T+1) = [0; 0; 20*w_1(T); 20*w_2(T)];

  %calculate lambda from the adjoint equations
    for i= 1:T
        j = T+1-i;
        grad_L = [  %use this if obj function depends on state
    0;                 % if v1 not in cost
    0;                 % if v2 not in cost
    2*w_1(j);          % derivative of w1^2
    2*w_2(j)           % derivative of w2^2
];
    %grad_L = zeros(4,1);  % use this if obj depends only on cost

    % Compute Jacobian of the dynamics at time j
    Jac = diff_1_2p_joint(v_1(j),v_2(j),w_1(j),w_2(j),u_1(j),u_2(j),mu1,mu2,rho1,rho2,alpha1,alpha2,timing_1,timing_2,p1,p2);

    % Costate/Adjoint update:
    lambda(:, j) = grad_L + Jac' * lambda(:, j+1);
    end
utemp_1 = zeros(1,T);
utemp_2 = zeros(1,T);

 %Calculate new u from the characterization of the optimal control
   for i=1:T
      %x0 = u_1(i);
      lambda_v_1 = lambda(1, i+1);
      lambda_v_2 = lambda(2, i+1);
      lambda_w_1 = lambda(3, i+1);
      lambda_w_2 = lambda(4, i+1);
      %f1 = @(u1) myfun_1(v_1(i),v_2(i),w_1(i),w_2(i),u1,u_2(i),mu,rho,alpha,lambda_v_1, lambda_v_2, lambda_w_1, lambda_w_2,timing_1,timing_2,p1,p2);
      %f2 = @(u2) myfun_2(v_1(i),v_2(i),w_1(i),w_2(i),u_1(i),u2,mu,rho,alpha, lambda_v_1, lambda_v_2, lambda_w_1, lambda_w_2, timing_1,timing_2,p1,p2);
      %val = f1(x0);


      %zero_h_1=fzero(f1,x0);
      %zero_h_2=fzero(f2,x0);
      %utemp_1(i) = max(0,min(0.9,zero_h_1));
      %utemp_2(i) = max(0,min(0.9,zero_h_2));
      step = 0.01;
grad1 = dif_u_patch_1( v_1(i),v_2(i),w_1(i),w_2(i), u_1(i),u_2(i), mu1,mu2,rho1,rho2, alpha1,alpha2, lambda_v_1,lambda_v_2,lambda_w_1,lambda_w_2, timing_1,timing_2,p1,p2);

grad2 = dif_u_patch_2(v_1(i),v_2(i),w_1(i),w_2(i), u_1(i),u_2(i), mu1,mu2,rho1,rho2, alpha1,alpha2, lambda_v_1,lambda_v_2,lambda_w_1,lambda_w_2, timing_1,timing_2,p1,p2);

step = 0.01;

u_try_1 = u_1(i) - step * grad1;
u_try_2 = u_2(i) - step * grad2;

utemp_1(i) = min(0.9, max(0, u_try_1));
utemp_2(i) = min(0.9, max(0, u_try_2));
   end

  %determine a convex update for the u vector
  gamma = 0.5; % start small, maybe increase to 0.5
    u_1 = gamma * utemp_1 + (1 - gamma) * oldu_1;
    u_2 = gamma * utemp_2 + (1 - gamma) * oldu_2;


  %compare results using a convergence test
    %Here we test how much a variable changes according to it's size
 changes = [
        delta * sum(abs(oldu_1)) - sum(abs(oldu_1 - u_1));
        delta * sum(abs(oldu_2)) - sum(abs(oldu_2 - u_2));
        delta * sum(abs(oldw_1)) - sum(abs(oldw_1 - w_1));
        delta * sum(abs(oldw_2)) - sum(abs(oldw_2 - w_2));
        delta * sum(abs(oldv_1)) - sum(abs(oldv_1 - v_1));
        delta * sum(abs(oldv_2)) - sum(abs(oldv_2 - v_2));
        delta * sum(abs(oldlambda(1,:))) - sum(abs(oldlambda(1,:) - lambda(1,:)));  % v1 costate
        delta * sum(abs(oldlambda(2,:))) - sum(abs(oldlambda(2,:) - lambda(2,:)));  % v2 costate
        delta * sum(abs(oldlambda(3,:))) - sum(abs(oldlambda(3,:) - lambda(3,:)));  % w1 costate
        delta * sum(abs(oldlambda(4,:))) - sum(abs(oldlambda(4,:) - lambda(4,:)));  % w2 costate
    ];

    test = min(changes);  % track least converged variable
    fprintf('Min convergence gap: %.6f\n', test);
  end

%states and time
y(1,:) = t;
y(2,:) = v_1;
y(3,:) = v_2;
y(4,:) = w_1;
y(5,:) = w_2;

%control
z(1,:) = t(1:T);
z(2,:) = u_1;
z(3,:) = u_2;
