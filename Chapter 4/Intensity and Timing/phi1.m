function [v1_out,c1_out,b1_out,v2_out,c2_out,b2_out] = phi1( v1,v2,w1,w2,s,mu1,mu2,rho1,rho2,alpha1,alpha2)

N = 400;
t = linspace(0,s,N);

% =========================================================
% PATCH 1
% =========================================================

A1 = rho1*t - (alpha1*w1/mu1).*(1-exp(-mu1*t));
expA1 = exp(A1);

% Inner integral
inner1 = cumtrapz(t,expA1);

% Resource trajectory v1(t)
v1_vals = v1 .* expA1 ./ (1 + rho1*v1.*inner1);

% Resource at time s
v1_flow = v1_vals(end);

% Consumer at time s
c1_flow = w1*exp(-mu1*s);

% Birth accumulation
b1_flow = alpha1*trapz(t,v1_vals);

% =========================================================
% PATCH 2

A2 = rho2*t - (alpha2*w2/mu2).*(1-exp(-mu2*t));
expA2 = exp(A2);

% Inner integral
inner2 = cumtrapz(t,expA2);

% Resource trajectory v2(t)
v2_vals = v2 .* expA2 ./ (1 + rho2*v2.*inner2);

% Resource at time s
v2_flow = v2_vals(end);

% Consumer at time s
c2_flow = w2*exp(-mu2*s);

% Birth accumulation
b2_flow = alpha2*trapz(t,v2_vals);

% =========================================================
% OUTPUT

v1_out = v1_flow;
c1_out = c1_flow;
b1_out = b1_flow;

v2_out = v2_flow;
c2_out = c2_flow;
b2_out = b2_flow;



end
