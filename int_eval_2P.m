function outer_int=int_eval_2P(rho, alpha, mu, w, v, ww, vv,nodes)
  for i=1:length(nodes); %index through length of node array
    inner_1 = @(t) exp(rho*t-((alpha*w/mu)*(1-exp(-mu.*t))));
    inner_2 = @(t) exp(rho*t-((alpha*ww/mu)*(1-exp(-mu.*t))));
    inner_int_1(i)=integral(inner_1,0,nodes(i));
    inner_int_2(i)=integral(inner_2,0,nodes(i));
    outer_1(i)= (exp(rho*nodes(i)-(alpha*w/mu)*((1-exp(-mu.*nodes(i))))))./((rho.*v.*inner_int_1(i))+1);
    outer_2(i)= (exp(rho*nodes(i)-(alpha*ww/mu)*((1-exp(-mu.*nodes(i))))))./((rho.*vv.*inner_int_2(i))+1);
  endfor
%Calculates the trapezoidal rule for the outer integral
outer_int_1=trapz(nodes,outer_1);
outer_int_2=trapz(nodes,outer_2);
outer_int=[outer_int_1,outer_int_2];
