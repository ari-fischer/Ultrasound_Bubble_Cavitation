function [Cp_nd,Cv_nd,H_nd,S_nd,G_nd]=thermo_data(T_nd,thermo_param)
% THERMO_DATA takes in a non-dimensionalized temperature (scalar quantity)
% and outputs various thermodynamic properties of the species involved in
% the water reaction network.
%
% NASA Polynomials coefficients, universal gas constant, and T0 hould be
% imported as global variable 
%
% NASA Polynomials obtained from https://ntrs.nasa.gov/api/citations/20020085330/downloads/20020085330.pdf

R_ig = thermo_param{1};
T0 = thermo_param{2};
P0 = thermo_param{3};
n0 = thermo_param{4};
R0 = thermo_param{5};
NASA_table1 = thermo_param{6};
NASA_table2 = thermo_param{7};
NASA_table3 = thermo_param{8};

T_NASA=T_nd*T0;

if T_NASA<=1000
    NASA_coeffs=table2array(NASA_table1(:,2:end));
elseif T_NASA>1000 && T_NASA<=6000
    NASA_coeffs=table2array(NASA_table2(:,2:end));
else 
    NASA_coeffs=table2array(NASA_table3(:,2:end));
    1;
    %NASA_coeffs_LB=table2array(NASA_table2(:,2:end));
    %HT = 1;
end

% Cp/R = a1*T^-2 + a2*T^-1 + a3 + a4*T + a5*T^2 + a6*T^3 + a7*T^4
Cp_NASA=[T_NASA^-2,T_NASA^-1,1,T_NASA,T_NASA^2,T_NASA^3,T_NASA^4,0,0]';

Cp_0R=NASA_coeffs*Cp_NASA; %should give a column vector;
% has dimensions of Cp0/R
Cv_0R = Cp_0R-1;

% remove dimensions Cp0/R by dividing by n0

Cp_nd=Cp_0R; %non-dimensionalize
% Cp_0R (dimensionless) -> Cp_0R.*R_ig (J/mol/K) ->

Cv_nd=Cv_0R;

% H/RT = -a1*T^-2 + a2*ln(T)/T + a3 + a4*T/2 + a5*T^2/3 + a6*T^3/4
% + a7*T^4/5 + b1/T

H_NASA=[-T_NASA^-2,log(T_NASA)/T_NASA,1,T_NASA/2,T_NASA^2/3,T_NASA^3/4,T_NASA^4/5,1/T_NASA,0]';

H_0RT=NASA_coeffs*H_NASA; %likewise should give column vector. 
% dimeneions of H0/RT -> multiply bt T_NASA and divide by n0 T0
H_nd = H_0RT.*T_NASA./T0;
% S/R = -a1*T^-2/2 - a2*T^-1 + a3ln(T) + a4*T + a5*T^2/2 + a6*T^3/3
% + a7*T^4/4 + b2

S_NASA=[-T_NASA^-2/2,-T_NASA^-1,log(T_NASA),T_NASA,T_NASA^2/2,T_NASA^3/3,T_NASA^4/4,0,1]';

S_0R=NASA_coeffs*S_NASA; %also should give a column vector
% demsnsions of S/R
S_nd=S_0R;

% G=H-TS
G=H_0RT*R_ig*T_NASA-S_0R.*T_NASA.*R_ig; %units J/mol
G_nd=G./R_ig./T0; %units J/mol
%n_H2O;n_OH;n_H;n_H2;n_H2O2;n_O;n_HO2;n_O2;n_Ar
1;
