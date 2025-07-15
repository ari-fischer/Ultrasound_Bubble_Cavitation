%font = type
%font size = 

figure('Position',[400 400 460 500]);
subplot(2,2,2)
box on
hold on

plot(PA_space/1E5,OH_yields_M_nd.*n0_M*6.022E23.*1E-8,'-.',...
    'LineWidth',1,'Color',	"#D95319")
plot(PA_space/1E5,H2O_cons_M.*1E-8,'-',...
    'LineWidth',1,'Color',	"#0072BD")
hold off
%ylim([0,16])
xlabel('Pressure amplitude (bar)')
ylabel({'Molecules of H_{2}O consumed (\epsilon_{H2O} ×10^{-8})';'or OH formed (N_{OH} ×10^{-8}) '})

fontsize(gcf,10,"points")
fontname(gcf,"Arial")
xticks(0:2:10)

%%


%%
subplot(2,2,1)
box on

hold on
%plot(PA_space/1E5,total_work_M,'--',...
%    'LineWidth',1.25,'Color',	"#0072BD")
plot(PA_space/1E5,work_done_M.*1E10,'-',...
    'LineWidth',1,'Color',	"#0072BD")
plot(PA_space/1E5,E_therm_M.*1E10,'-.',...
    'LineWidth',1,'Color',	"#D95319")
% plot(PA_space/1E5,E_therm_H2O_M,'--',...
%     'LineWidth',1.25,'Color',	"#D95319")
plot(PA_space/1E5,E_rxn_M.*1E10,'--',...
    'LineWidth',1,'Color',	"#77AC30")
ylabel('Energy (J ×10^{10})')
ax = gca; % Get the current axes
ax.YColor = 'k'; % Set the left y-axis color to blue
xlabel('Pressure amplitude (bar)')
hold off

fontsize(gcf,10,"points")
fontname(gcf,"Arial")
xticks(0:2:10)


subplot(2,2,3)
x=E_rxn_M,
%y=OH_yields_M_nd.*n0_M*6.022E23
y = H2O_cons_M.*1E-8

A = sum(x .* y) / sum(x.^2);
A2 = sum(x .* H2O_cons_M) / sum(x.^2)
%-------------------------------
% Step 2: Generate fitted y values using the computed slope
y_fit = A * x;

box on
hold on
plot(x.*10.^10,H2O_cons_M.*1E-8,'o','Color',"#0072BD","MarkerFaceColor","#0072BD",'MarkerSize',6)
%plot(x,y,'o','Color',"#D95319","MarkerFaceColor","#D95319",'MarkerSize',6)

%plot(x, y_fit, '-','Color',"#D95319", 'LineWidth', 1);
plot(x.*10.^10, x.*A2.*1E-8, '-','Color',"#0072BD", 'LineWidth', 1);
xlabel('Reaction energy, {\Delta}E_{rxn} (J ×10^{10})')
ylabel('Molecules of H_{2}O consumed (\epsilon_{H2O} ×10^{-8})')
ax = gca; % Get the current axes
ax.YColor = 'k'; % Set the left y-axis color to blue
hold off

fontsize(gcf,10,"points")
fontname(gcf,"Arial")
xlim([0,5])
xticks(0:1:5)


y_mean = mean(y);                        % Mean of the observed y-values
SS_tot = sum((y - y_mean).^2);             % Total sum of squares
SS_res = sum((y - y_fit).^2);              % Sum of squared residuals from the fit
R_squared = 1 - (SS_res / SS_tot);     

%%



PA_Emax = PA_space(find(E_rxn_M==max(E_rxn_M)));




f = gcf;
exportgraphics(f,'fig_2_30Apr.png','Resolution',300)
