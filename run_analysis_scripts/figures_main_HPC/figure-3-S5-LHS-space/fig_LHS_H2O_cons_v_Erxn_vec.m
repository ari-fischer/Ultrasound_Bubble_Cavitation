%read Latin hypercube sampling data
x = readmatrix("HPC-gen-data/E_rxn.csv");
y = readmatrix("HPC-gen-data/H2O_extent.csv");

%fit line
p = fitlm(x, y,  'y ~ x1 - 1');
slope = p.Coefficients.(1);
span = linspace(min(x(x~=0))/2,max(x)*2,20);

%Vary pressure amplitude (copy data from Fig 2)
H2O_cons_PA = readmatrix("out_data/H2O_cons_molec.csv");
E_rxns_PA = readmatrix("out_data/rxn-energy.csv");

%make figure
figure('Position',[400 400 240 240]);

box on
hold on
%linear fn
plot(span,span.*slope,'color','k','LineWidth',.75)
%fints
scat = scatter(x,y,60,'filled')
scat.MarkerFaceColor = "#0072BD"
scat.MarkerFaceAlpha = .7

scat = scatter(E_rxns_PA,H2O_cons_PA,120,'filled',"^")
scat.MarkerFaceColor = "#D95319"
scat.MarkerFaceAlpha = .7

hold off
set(gca,'YScale', 'log','YLim',[1e2 1e10],'YTick',10.^(2:10), ...
        'XScale', 'log','XLim',[1e-16 1e-8],'XTick',10.^(-16:-8))

ax=gca;
ax.FontSize=12;
xlabel('Total reaction energy ({\Delta}E_{rxn}; J)')
ylabel('H_{2}O molecules consumed (\epsilon_{H2O})')

fontsize(gcf,10,"points")
fontname(gcf,"Arial")

f = gcf;
exportgraphics(f,'fig_4_LHS.png','Resolution',300)
