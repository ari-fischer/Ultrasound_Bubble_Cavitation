%font = type
%font size = 

figure('Position',[400 400 800 220]);
subplot(1,3,1)
box on
hold on
plot(time_profile(1:ad_switch_idx)*1E6,radial_profile(1:ad_switch_idx)*1E6,...
    'LineWidth',1,'Color',	"#0072BD")
plot(time_profile(ad_switch_idx:end)*1E6,radial_profile(ad_switch_idx:end)*1E6,...
    '--','LineWidth',1,'Color',	"#0072BD")
plot(time_profile(ad_switch_idx)*1E6,...
    radial_profile(ad_switch_idx)*1E6,'o','MarkerFaceColor','Black', ...
    'MarkerEdgeColor','Black')
%xlim([0,max(time_profile(1:ad_switch_idx))*1E6*1.2])
xlim([0,2.5])
set(gca,'XTick',(0:0.5:2.5))
xlabel('Time (\mus)')
ylabel('Radius (\mum)')
hold off

fontsize(gcf,10,"points")
fontname(gcf,"Arial")

% work and temperature plot

subplot(1,3,2)
box on
yyaxis left
hold on
plot(time_profile(1:ad_switch_idx+1)*1E6,...
    time_profile(1:ad_switch_idx+1).*0, ...
    '-','LineWidth',1,'Color',	"#0072BD")

plot(time_profile(ad_switch_idx+1:end)*1E6,cum_work.*1E10...
    ,'--','LineWidth',1,'Color',	"#0072BD")
ylabel('Cumulative adiabatic work (J ×10^{10})')
ax = gca; % Get the current axes
ax.YColor = 'k'; % Set the left y-axis color to blue
%ylim([0,16])
xlabel('Time (\mus)')
hold off

fontsize(gcf,10,"points")
fontname(gcf,"Arial")


subplot(1,3,2)

hold on
yyaxis right

plot(time_profile(1:ad_switch_idx)*1E6,T_profile(1:ad_switch_idx)...
    ,'-','LineWidth',1,'Color',"#D95319")

plot(time_profile(ad_switch_idx:end)*1E6,T_profile(ad_switch_idx:end)...
    ,'--','LineWidth',1,'Color',"#D95319")

hold off

ylabel('Temperature (K)')
ax = gca; % Get the current axes
ax.YColor = 'k'; % Set the left y-axis color to blue

xlim([time_profile(ad_switch_idx)*1E6*0.995,time_profile(end)*1E6*1.005])
set(gca,'XTick',(2.24:0.01:2.28))

fontsize(gcf,10,"points")
fontname(gcf,"Arial")

% H2O plot
subplot(1,3,3)
box on
hold on

plot(time_profile(1:ad_switch_idx)*1E6,mol_matrix_profile_nd(1:ad_switch_idx,1)...
    ,'-','LineWidth',1,'Color',	"#0072BD")

plot(time_profile(ad_switch_idx:end)*1E6,mol_matrix_profile_nd(ad_switch_idx:end,1)...
    ,'--','LineWidth',1,'Color',	"#0072BD")

ylim([0,max(mol_matrix_profile_nd(ad_switch_idx:end,1))*1.15])
xlabel('Time (\mus)')
ylabel('Moles H_{2}O (n_{H2O}/n_{0})')
ylim([0,0.25])
hold off

xlim([time_profile(ad_switch_idx)*1E6*0.995,time_profile(end)*1E6*1.005])
set(gca,'XTick',(2.24:0.01:2.28))

fontsize(gcf,10,"points")
fontname(gcf,"Arial")

f = gcf;
exportgraphics(f,'fig_collapse_rxn.png','Resolution',300)