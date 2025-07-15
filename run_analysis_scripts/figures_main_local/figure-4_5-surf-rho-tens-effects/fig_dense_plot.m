x = factor_l.*1000;
x_lab = 'Liquid density (kg m^{-3})'

figure('Position',[400 400 460 240]);

subplot(1,2,1)
box on

hold on
%plot(PA_space/1E5,total_work_M,'--',...
%    'LineWidth',1.25,'Color',	"#0072BD")
plot(x,work_done_M.*1E10,'-',...
    'LineWidth',1,'Color',	"#0072BD")% blue
plot(x,E_therm_M.*1E10,'-.',...
    'LineWidth',1,'Color',	"#D95319")% orange
% plot(PA_space/1E5,E_therm_H2O_M,'--',...
%     'LineWidth',1.25,'Color',	"#D95319")
plot(x,E_rxn_M.*1E10,'--',...
    'LineWidth',1,'Color',	"#77AC30")% green 
ylabel('Energy (J ×10^{10})')
xlabel(x_lab)
ax = gca; % Get the current axes
ax.YColor = 'k'; % Set the left y-axis color to blue
fontsize(gcf,10,"points")
fontname(gcf,"Arial")

xlim([700,1300])
hold off


subplot(1,2,2)
box on
hold on

plot(x,H2O_cons_M*1E-8,'-',...
    'LineWidth',1,'Color',	"#0072BD")
hold off
ylim([0,inf])
xlabel(x_lab)
ylabel('H_{2}O molecules consumed ({\epsilon}_{H2O} ×10^{-8})')

fontsize(gcf,10,"points")
fontname(gcf,"Arial")
xlim([700,1300])
ylim([0,2.5])
%%


f = gcf;
exportgraphics(f,'fig_4.png','Resolution',300)
