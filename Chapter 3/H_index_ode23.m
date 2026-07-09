
delta = 0.0:0.02:0.99;
years=300;
p_none=1;
orange  = [1, 0.33, 0];
bluee = [0.0117647, 0.0156863, 0.3686275];

% Baseline values (no dispersal)
[R_none, C_none] = run_model(p_none, years);

%Calculate H-index
for j = 1:length(delta)
prob=1-delta(j);
[R_delta,C_delta] = run_model(prob,years);

% compute H-index
H_R(j) = R_delta(end) - R_none(end);
H_C(j) = C_delta(end) - C_none(end);
end

plot(delta, H_C, 'Color', orange,  'LineWidth',2); hold on;
plot(delta, H_R, 'Color', bluee, 'LineWidth',2);
    yline(0,'k:');
    grid on;
lgd = legend('Consumers','Resources','Location', 'southwest');
%xlabel(t,'\delta (dispersal)');
%ylabel(t,'Benefit index H(\delta)');
xlabel('\delta (dispersal)');
ylabel('H-index')
ylim([-0.1 0.1])
exportgraphics(gcf,'H_index_f.pdf','ContentType','vector','BackgroundColor','none','Units', 'inches','Width', 6,'Height', 5)