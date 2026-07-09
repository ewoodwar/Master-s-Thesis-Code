%This code runs the model by Pachepsky et al. piecing together continuous and discrete
clear
clc
%Parameters
p1= 0.6;
p2= 0.2; %moving from 1 to 2

%initialization
p0=[1,0.2,0, 1,0.2,0];
years=1; %the amount of years the simulation runs

for i=0:years
  t0=i;                %Section 01 runs the continuous portion
  tf=i+1;
  tspan=[t0 tf];
  [t,p]=ode23(@pachepsky2Pf,tspan,p0); %runs the simulation for one year
  p_end=[p(end,1),p(end,2),p(end,3),p(end,4),p(end,5),p(end,6)]; %collects the values at the end of the year
  year_end_R_1(i+1)=p(end,1); %puts the year end value in the ith spot of an array (matlab starts at 1 index)
  year_end_C_1(i+1)=p(end,2);
  year_end_R_2(i+1)=p(end,4); %puts the year end value in the ith spot of an array (matlab starts at 1 index)
  year_end_C_2(i+1)=p(end,5);

                  %Section 02 iterates the map once (this is different than the regular pachepsky as we add migration once per year)
  p_end(2)=(p_end(3)+1)*p1*p(end,2) + (p_end(6)+1)*(1-p2)*p_end(5); %consumers patch 1 moving to patch 1 + consumers in patch 2 moving to patch 1
  p_end(5)=(p_end(6)+1)*p2*p_end(5) + (p_end(3)+1)*(1-p1)*p(end,2); %consumers patch 2
  p_end(3)=p_end(6)=0;
  after_C1(1)=p0(2);
  after_C2(1)=p0(5);
  after_R1(i+1)=p(end,1);
  after_C1(i+1)=p_end(2);
  after_R2(i+1)=p(end,4);
  after_C2(i+1)=p_end(5);

  p0=p_end;
  %figure(1)
  %plot(t,p(:,1), 'color','#77AC30')%Resource plot
  %hold on
  %plot(t,p(:,4), 'color','#AAFF00'), title('Resources'), xlabel('Days'), ylabel('Density'), axis([0 5 0 4])
  %legend('Resources Patch 1', 'Resources Patch 2')
  %hold on

  %figure(2)
  %plot(t, p(:,2),'color','#96DED1'),  %Consumer Plot Patch 1
  %hold on
  %plot(t,p(:,5), 'color','#808000'), title('Consumers'), xlabel('Days'), ylabel('Density'), axis([0 5 0 15]) %Consumer Plot Patch 2
  %legend('Consumers Patch 1', 'Consumers Patch 2')
  %hold on
end




hfig1 = figure('visible','off');


% --- SIZE OF FIGURE WINDOW ---
hw_ratio = 0.65;
picturewidth = 28;
set(hfig1,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth]);

soft_burg = [153/255 98/255 30/255];
soft_blue = [44/255 45/255 88/255];

j = 0:years;

% ========== PLOTTING ==========
subplot(1,2,1)
plot(j, year_end_R_1, 'Color', soft_burg, 'LineWidth',1.5); hold on
plot(j, year_end_R_2, 'Color', soft_blue, 'LineWidth',1.5);
xlabel('Years');
ylabel('Density');
title('Yearly Resource Density');
lgd = legend('Patch 1','Patch 2','location','southoutside');

subplot(1,2,2)
plot(j, year_end_C_1, 'Color', soft_burg, 'LineWidth',1.5); hold on
plot(j, year_end_C_2, 'Color', soft_blue, 'LineWidth',1.5);
xlabel('Years');
ylabel('Density');
title('Yearly Consumer Density');
lgd = legend('Patch 1','Patch 2','location','southoutside');

set(findall(hfig1,'Type','axes'), 'FontSize', 17, ...
    'TickLabelInterpreter', 'latex', 'Box','off');

% ========== FREEZE LAYOUT BEFORE SAVING ==========
set(hfig1,'Units','centimeters');
pos = get(hfig1,'Position');

set(hfig1,'PaperUnits','centimeters');
set(hfig1,'PaperSize',[pos(3) pos(4)]);
set(hfig1,'PaperPosition',[0 0 pos(3) pos(4)]);
set(hfig1,'PaperPositionMode','manual');

print(hfig1, 'overcompensation_1.pdf', '-dpdf', '-r300'); %CHANGE NAME OF FILE HERE

##hfig2=figure;
##hw_ratio=0.65;
##picturewidth=20;
##set(findall(hfig2,'-property','FontSize'),'FontSize',17);
##set(findall(hfig2,'-property','Box'), 'Box','off');
##set(findall(hfig2,'-property','Interpreter'), 'Interpreter','latex');
##set(findall(hfig2,'-property','TickLabelInterpreter'), 'TickLabelInterpreter','latex');
##set(hfig2,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth]);
##pos=get(hfig2,'Position');
##set(hfig2,'PaperPositionMode','Auto','PaperUnits','centimeters','Papersize',[pos(3), pos(4)]);
##
##k=0:years;
##%figure(4)
##subplot(1,2,1)
##plot(k, after_R1)
##hold on
##plot(k, after_R2), title('Yearly Resource Density (after reproduction)'), xlabel('Years'), ylabel('Density'), legend('Patch 1', 'Patch 2')
##subplot(1,2,2)
##plot(k, after_C1)
##hold on
##plot(k, after_C2), title('Yearly Consumer Density (after reproduction)'), xlabel('Years'), ylabel('Density'), legend('Patch 1', 'Patch 2')

