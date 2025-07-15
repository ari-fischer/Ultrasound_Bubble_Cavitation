
%ideal gas constant
R_ig = 8.2057366080960E-5;

% Avogadros number
N_av = 6.022E23;

%import H2O consumed (# of molecules) in dimensionless form
H2O_cons = readmatrix("rho_l_H2O_cons.csv");

%import sensitivities from perterbations to each property
Y1= readmatrix("PA_z_value.csv");
Y2 = readmatrix("rho_l_z_value");
Y3 = readmatrix("sigma_l_z_value");

%import frequencies from the space and use them to calculate the initial
%moles of gas in the bubble
omegas = readmatrix("rho_l_omega_M.csv");
R_res = 3./omegas
Vs = 4./3.*pi.*(R_res./2.1).^3
n0s = Vs.*1.101325./R_ig./298;

%use H2O consumed and initial gas molecules to get absolute number of H2O
%molecules consumed
x = H2O_cons.*n0s.*N_av;

%remove all datapoints with an amount of H2O consumed that falls well below
%the continuum limit
x(x<1E4)=NaN;
Y1(x<1E4)=NaN;
Y2(x<1E4)=NaN;
Y3(x<1E4)=NaN;

%make the figure. 
figure('Position',[400 400 260 300]);

box on
hold on
%linear fn
scat = scatter(x(:),abs(Y1(:)),40,'filled')
scat.MarkerFaceColor = "#0072BD"
scat.MarkerFaceAlpha = .7
scat.Marker='o'


scat2 = scatter(x(:),abs(Y2(:)),40,'filled')
scat2.MarkerFaceColor = "#D95319"
scat2.MarkerFaceAlpha = .7
scat2.Marker='square'

scat3 = scatter(x(:),abs(Y3(:)),40,'filled')
scat3.MarkerFaceColor = "#77AC30"
scat3.MarkerFaceAlpha = .7
scat3.Marker='^'



hold off
set(gca,'YScale', 'log','YLim',[10^(-4) 1e2],'YTick',10.^(-4:2), ...
        'XScale', 'log','XLim',[1e4 1e10],'XTick',10.^(4:10))

ax=gca;
ax.FontSize=12;
xlabel('Molecules of H_{2}O consumed ({\epsilon}_{H2O})')
ylabel('Absolute normalized sensitivity (|S^{\wedge}_{H2O}|)')


fontsize(gcf,10,"points")
fontname(gcf,"Arial")

f = gcf;
exportgraphics(f,'fig_7_Erxn.png','Resolution',300)
