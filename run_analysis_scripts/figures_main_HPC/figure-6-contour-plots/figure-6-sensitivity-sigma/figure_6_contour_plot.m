close all
% Script plots out the contour plots of obtained from the simulations

% Dictionary mapping all the strings to their parameter names
var_dict=containers.Map(["Density","Viscosity","Surface Tension","Pressure Amplitude"],["rho_l","v_l","sigma_l","PA"]);



z=readmatrix("sigma_l_z_value");
r=readmatrix("sigma_l_Rsq.csv");
omega_M=readmatrix("sigma_l_omega_M.csv");
PA_M=readmatrix("sigma_l_PA_M.csv");
H2O_cons=readmatrix("sigma_l_H2O_cons.csv");
H2Os=readmatrix("sigma_l_H2O_init.csv");

omega_plot=omega_M/1000;
PA_plot=PA_M/10^5;

r(r<0.99)=NaN;
%z is omega columns, PA rows. need to transpose depending on shape we want.

%ideal gas constant
R_ig = 8.2057366080960E-5;

% Avogadros number
N_av = 6.022E23;
R_res = 3./omega_M
Vs = 4./3.*pi.*(R_res./2.1).^3
n0s = Vs.*1.101325./R_ig./298;

%use H2O consumed and initial gas molecules to get absolute number of H2O
%molecules consumed
x = H2O_cons.*n0s.*N_av;

z_values_plot=z;
%z_values_plot(convs<1E-6)
z_values_plot(x<1E4)=NaN;

% Contour plot


figure('Position',[400 400 260 240])
contourf(omega_plot,PA_plot,z_values_plot,15);
%colormap(flipud(autumn))
%caxis([0 max(max(z_values_plot))*1.1])  

colormap(summer)
caxis([-20 0])  

colorbar
xlabel('\omega (kHz)')
ylabel('P_A (bar)')

xticks(220:200:1020)
fontsize(gcf,10,"points")
fontname(gcf,"Arial")

f = gcf;
exportgraphics(f,'contour_sigma.png','Resolution',300)


