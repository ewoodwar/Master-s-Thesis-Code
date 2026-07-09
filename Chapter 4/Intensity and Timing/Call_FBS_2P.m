tic
alpha1=10;
alpha2=11;
rho1=20;
rho2=25;
mu1=1.5;
mu2=1;
T=30; #number of years (finite time horizon)
p1=0.15;
p2=0.65;


#Call the function that calculates FBS with this command_line_path
#function [y,z] = FBS_twopatch_JOINT(mu1,mu2,rho1,rho2, alpha1,alpha2, init_c_1, init_c_2,init_r_1,init_r_2,p1,p2,T)


[y1,z1]=FBS_2p_TI(mu1,mu2, rho1,rho2, alpha1,alpha2,1,1,0.2,0.2,p1,p2,T);


years= 0:T-1; %creates a vector with the years
time=z1(1,:); %takes the output from the length of FBS

control_p1_I=z1(2,:); %vector of yearly intensities of control p1 (between 0-1)
control_p2_I=z1(3,:);
control_p1_T=z1(4,:); %vector of yearly control timing (between 0,1) p1 example (0.5,0.8,0.9, ..)
control_p2_T=z1(5,:);

%set up control timepoints with year starting from zero
timing_p1=years+control_p1_T; %sets up the actual timepoints of control, with the timing starting from zero ex. (0.5,1.8,2.9...)
timing_p2=years+control_p2_T;


%set up figure specifications
hfig=figure;
hw_ratio=0.65;
picturewidth=20;
set(findall(hfig,'-property','FontSize'),'FontSize',17);
set(findall(hfig,'-property','Box'), 'Box','off');
set(findall(hfig,'-property','Interpreter'), 'Interpreter','latex');
set(findall(hfig,'-property','TickLabelInterpreter'), 'TickLabelInterpreter','latex');
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth]);
pos=get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','Papersize',[pos(3), pos(4)]);
soft_blue = [0.4 0.6 0.8];
soft_green = [0.6 0.8 0.6];
soft_black = [0.3 0.3 0.3];


%Plotting C/R Dynamics (same as before)
subplot(2,2,1) %patch1
plot(y1(1,:), y1(2,:), '-o', 'Color', soft_blue, 'MarkerFaceColor', soft_blue);
hold on
plot(y1(1,:), y1(4,:), '-o', 'Color', soft_green, 'MarkerFaceColor', soft_green);
xlabel('Years')
ylabel('Density')
title('Population Dynamics Patch 1')
%legend('Resource','Consumer')
subplot(2,2,2)%patch2
plot(y1(1,:), y1(3,:), '-o', 'Color', soft_blue, 'MarkerFaceColor', soft_blue);
hold on
plot(y1(1,:), y1(5,:), '-o', 'Color', soft_green, 'MarkerFaceColor', soft_green);
xlabel('Years')
ylabel('Density')
title('Population Dynamics Patch 2')
%legend('Resource','Consumer')

%Plotting Control Intensity and Timing (changed structure slightly)

%patch 1
subplot(2,2,3)
picturewidth=20;
scatter(timing_p1, control_p1_I)
xlabel('Years')
ylabel('Control Intensity')
xlim([0, length(time)]);
title('Control Intensity')
hold on
for i = 1:length(timing_p1)
    plot([timing_p1(i), timing_p1(i)], [0, control_p1_I(i)], 'Color', soft_black, 'LineWidth', 2);
end

%patch 2
subplot(2,2,4)
picturewidth=20;
scatter(timing_p2, control_p2_I)
xlabel('Years')
ylabel('Control Intensity')
xlim([0, length(time)]);
title('Control Intensity')
hold on
for i = 1:length(timing_p2)
    plot([timing_p2(i), timing_p2(i)], [0, control_p2_I(i)], 'Color', soft_black, 'LineWidth', 2);
end
toc
