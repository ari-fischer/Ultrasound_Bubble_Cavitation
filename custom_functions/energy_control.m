%net energy consumed for each reversible reaction
E_rxn_elem_M_net = E_rxn_elem_M(:,1:2:end-1)+E_rxn_elem_M(:,2:2:end)
%average energy per E_rxn
E_elem_ave = mean(E_rxn_elem_M_net,1)

[~, sortedIndices] = sort(abs(E_elem_ave), 'descend');

% Extract the top 5 indices
Indices = sortedIndices(1:6);

figure('Position',[400 400 260 300]);
box on
hold on
plot(PA_space/1E5,E_rxn_M ,'--',...
    'LineWidth',1.25,'Color',	"#77AC30")
plot(PA_space/1E5,E_rxn_elem_M_net(:,Indices(1)),'-',...
    'LineWidth',1.25,'Color',	"#77AC30")
plot(PA_space/1E5,E_rxn_elem_M_net(:,Indices(2)),'-',...
    'LineWidth',1.25,'Color',	"#77AC30")
plot(PA_space/1E5,E_rxn_elem_M_net(:,Indices(3)),'-',...
    'LineWidth',1.25,'Color',	"#77AC30")
plot(PA_space/1E5,E_rxn_elem_M_net(:,Indices(4)),'-',...
    'LineWidth',1.25,'Color',	"#77AC30")
yline(0)
%plot(PA_space/1E5,E_rxn_elem_M_net(:,Indices(5)),'-',...
%    'LineWidth',1.25,'Color',	"#77AC30")
%plot(PA_space/1E5,E_rxn_elem_M_net(:,Indices(6)),'-',...
%    'LineWidth',1.25,'Color',	"#77AC30")

%xlim([0,8])

hold off
%ylim([0,16])
xlabel('Pressure amplitude (bar)')
ylabel('Reaction energy (J)')


f = gcf;
exportgraphics(f,'fig_3.png','Resolution',300)


max(E_rxn_elem_M_net(2:end,Indices(1))./E_rxn_M(2:end))
max(E_rxn_elem_M_net(2:end,Indices(2))./E_rxn_M(2:end))
max(abs(E_rxn_elem_M_net(2:end,Indices(3)))./E_rxn_M(2:end))
max(E_rxn_elem_M_net(2:end,Indices(4))./E_rxn_M(2:end))