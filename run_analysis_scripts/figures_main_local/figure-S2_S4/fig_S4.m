%net energy consumed for each reversible reaction
OH_yields_rates_M_net = OH_yields_rates_M(:,1:2:end-1)+OH_yields_rates_M(:,2:2:end)
%average energy per E_rxn
OH_yields_ave = mean(OH_yields_rates_M_net,1)

[~, sortedIndices] = sort(abs(OH_yields_ave), 'descend');

OH_tot = sum(OH_yields_rates_M,2);

% Extract the top 5 indices
Indices = sortedIndices(1:6);

figure('Position',[400 400 260 300]);
box on
hold on
plot(PA_space/1E5,OH_tot ,'--',...
    'LineWidth',1.25,'Color',	"#77AC30")
plot(PA_space/1E5,OH_yields_rates_M_net(:,sortedIndices(1)),'-',...
    'LineWidth',1.25,'Color',	"#77AC30")
plot(PA_space/1E5,OH_yields_rates_M_net(:,sortedIndices(2)),'-',...
    'LineWidth',1.25,'Color',	"#77AC30")
plot(PA_space/1E5,OH_yields_rates_M_net(:,sortedIndices(3)),'-',...
    'LineWidth',1.25,'Color',	"#77AC30")
plot(PA_space/1E5,OH_yields_rates_M_net(:,sortedIndices(4)),'-',...
    'LineWidth',1.25,'Color',	"#77AC30")
plot(PA_space/1E5,OH_yields_rates_M_net(:,sortedIndices(5)),'-',...
    'LineWidth',1.25,'Color',	"#77AC30")
plot(PA_space/1E5,OH_yields_rates_M_net(:,sortedIndices(6)),'-',...
    'LineWidth',1.25,'Color',	"#77AC30")
yline(0)

%xlim([0,8])

hold off
xlabel('Pressure amplitude (bar)')
ylabel('OH yield (moles)')


f = gcf;
exportgraphics(f,'fig_S4.png','Resolution',300)
