

figure('Position',[400 400 375/2 435/2]);


%%
box on

hold on
%plot(PA_space/1E5,total_work_M,'--',...
%    'LineWidth',1.25,'Color',	"#0072BD")
plot(PA_space/1E5,S_work*1E5*1E10,'-',...
    'LineWidth',1,'Color',	"#0072BD")
plot(PA_space/1E5,S_therm*1E5*1E10,'-.',...
    'LineWidth',1,'Color',	"#D95319")
% plot(PA_space/1E5,E_therm_H2O_M,'--',...
%     'LineWidth',1.25,'Color',	"#D95319")
plot(PA_space/1E5,S_rxn*1E5*1E10,'--',...
    'LineWidth',1,'Color',	"#77AC30")
plot(0:1:10,[0:1:10]*.0,'-',...
    'LineWidth',0.5,'Color',	'k')
ylabel('Sensitivity (J ×10^{10} bar^{-1})')
ax = gca; % Get the current axes
ax.YColor = 'k'; % Set the left y-axis color to blue
xlabel('Pressure amplitude (bar)')
hold off
fontsize(gcf,10,"points")
fontname(gcf,"Arial")
xticks(0:2:10)


f = gcf;
exportgraphics(f,'fig_2d_15May.png','Resolution',300)
