x = factor_l.*72;

figure('Position',[400 400 760 260]);

%%

subplot(1,3,1)

box on
%yyaxis left
hold on
plot(x,H2O_inits_M.*N_A.*1E-8,'-.',...
    'LineWidth',1,'Color',	"#D95319")
plot(x,N_A.*N_Ar_M(:).*1E-8,'-',...
    'LineWidth',1,'Color',	"#0072BD")

xlabel('Surface tension (mN m^{-1})')
ylabel('Molecules in bubble (×10^{-8})')
ylim([0,5.2E1])
ax = gca; % Get the current axes
ax.YColor = 'k'; % Set the left y-axis color to blue

hold off
ax = gca; % Get the current axes
ax.YColor = 'k'; % Set the left y-axis color to blue

hold off
fontsize(gcf,10,"points")
fontname(gcf,"Arial")

%%
subplot(1,3,2)
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
ax = gca; % Get the current axes
ax.YColor = 'k'; % Set the left y-axis color to blue
xlabel('Surface tension (mN m^{-1})')
hold off
fontsize(gcf,10,"points")
fontname(gcf,"Arial")

subplot(1,3,3)
box on
hold on

plot(x,H2O_cons_M.*1E-8,'-',...
    'LineWidth',1,'Color',	"#0072BD")
hold off
ylim([0, inf])
xlabel('Surface tension (mN m^{-1})')
ylabel('H_{2}O molecules consumed (\epsilon_{H2O} ×10^{-8})')

fontsize(gcf,10,"points")
fontname(gcf,"Arial")

f = gcf;
exportgraphics(f,'fig_5.png','Resolution',300)
