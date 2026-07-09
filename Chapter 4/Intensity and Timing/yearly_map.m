function x_next = yearly_map(x,Lambda1,Lambda2,u1,u2,mu1,mu2,rho1,rho2,alpha1,alpha2,p1,p2)

% Determine order of spray events

if Lambda1 <= Lambda2

    firstTime  = Lambda1;
    secondTime = Lambda2;

    firstPatch  = 1;
    secondPatch = 2;

else

    firstTime  = Lambda2;
    secondTime = Lambda1;

    firstPatch  = 2;
    secondPatch = 1;

end

% Initial state
v1 = x(1);
v2 = x(2);
c1 = x(3);
c2 = x(4);


% Grow to first spray

[v1,c1,b1a,v2,c2,b2a] = phi1(v1,v2,c1,c2,firstTime,mu1,mu2,rho1,rho2,alpha1,alpha2);

%-------------------------------------------------
% First spray

if firstPatch == 1
    c1 = (1-u1)*c1;
else
    c2 = (1-u2)*c2;
end

%-------------------------------------------------
% Grow to second spray

[v1,c1,b1b,v2,c2,b2b] = phi1(v1,v2,c1,c2,secondTime-firstTime,mu1,mu2,rho1,rho2,alpha1,alpha2);

% Second spray

if secondPatch == 1
    c1 = (1-u1)*c1;
else
    c2 = (1-u2)*c2;
end

% Grow to end of year

[v1_end,c1_end,b1c,v2_end,c2_end,b2c] = phi1(v1,v2,c1,c2,1-secondTime,mu1,mu2,rho1,rho2,alpha1,alpha2);

b1 = b1a + b1b + b1c;
b2 = b2a + b2b + b2c;

% Reproduction
c1_rep = (1+b1)*c1_end;
c2_rep = (1+b2)*c2_end;


% Migration

c1_next = p1*c1_rep + (1-p2)*c2_rep;
c2_next = (1-p1)*c1_rep + p2*c2_rep;

% Output state
x_next = [v1_end;
          v2_end;
          c1_next;
          c2_next];

end
