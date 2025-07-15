function [rho,v,beta,sigma]=temp_liq_properties(T,P)
% TEMP_LIQ_PROPERTIES are correlations obtained from literature relating
% the temperature (K) and pressure (kPa) to the density, viscosity, bulk
% modulus and surface tension of water

%density
rho=-0.004988*T^2+2.73*T+626.5;

%kinematic viscosity
v=4.038E-10*T^2-2.638E-7*T+4.364E-5;

%bulk modulus (P in kPa)
a=8.435281+4.585473*10^-5*P;
b=9.603812*10-3.264414*10^-2*P;
c=4.242841*10^5+7.846325*P;
d=-9.342021*10^7-5.877829*10^2*P;
beta=exp(a+b/T+c/T^2+d/T^3)*10^5;

%surface tension
sigma=-7.781E-7*T^2+0.0003108*T+0.04946;

%update to use temperatures at 298 K for manuscript:
rho = 997; %kg m-3
v = 0.889/rho/1000; %kinimatic viscosity, Pa s / rho
sigma = 72/1000 ; %N m-1