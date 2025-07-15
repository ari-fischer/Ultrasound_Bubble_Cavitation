nAr_0 = 8.1399e-15;

figure('Position',[400 400 500 260]);

subplot(1,2,1)
box on

hold on
plot(PA_space/1E5,E_therm_M,'-.',...
    'LineWidth',1,'Color',	"#D95319")

ylabel('{\Delta}E_{{\Delta}T} (J)')
ax = gca; % Get the current axes
ax.YColor = 'k'; % Set the left y-axis color to blue
xlabel('Pressure amplitude (bar)')
hold off

%%

subplot(1,2,2)
box on
hold on

plot(PA_space/1E5,H2O_inits_M.*0+nAr_0,'-.',...
    'LineWidth',1,'Color',	"#D95319")
plot(PA_space/1E5,H2O_inits_M,'-',...
    'LineWidth',1,'Color',	"#0072BD")
hold off
%ylim([0,16])
xlabel('Pressure amplitude (bar)')
ylabel('Molecules of Ar and H_{2}O')% or OH formed (N_{OH})')




f = gcf;
exportgraphics(f,'fig_S2_30Apr.png','Resolution',300)
