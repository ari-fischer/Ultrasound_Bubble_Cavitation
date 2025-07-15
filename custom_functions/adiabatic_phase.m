function dy=adiabatic_phase(t,y,adia_param,kin_param)
% Mod 2 Apr 2025

% ADIABATIC_PHASE is a function that describes momentum, energy and
% mole balances during the adiabatic phase of the bubble collapse in the
% S&S model

global t0 t_1 t_2 t_3 t_4 gamma_4 gamma_3 gamma_2 gamma_1 gamma0 ...
    pgdotterm10 pgdotterm1_1 pgdotterm1_2 pgdotterm1_3 pgdotterm1_4 

% collecting y variables
R=y(1); % non-dimensionalized bubble radius
dR=y(2); % derivative of the non-dimensionalized bubble radius
n_H2O=y(3); % non-dimensionalized number of moles of H2O
T_nd=y(4); % non-dimensionalized temperature
n_OH=y(5); % non-dimentionalized number of moles of OH
n_H=y(6); % non-dimentionalized number of moles of H
n_H2=y(7); % non-dimentionalized number of moles of H2
n_H2O2=y(8); % non-dimentionalized number of moles of H2O2
n_O=y(9); % non-dimentionalized number of moles of O
n_HO2=y(10); % non-dimentionalized number of moles of HO2
n_O2=y(11); % non-dimentionalized number of moles of O2


%n_Ar Mw_vector T0 R0 omega ...
%    Eu We Rel Mal PA_nd knd_mod
knd_mod = adia_param{1}; 
n_Ar = adia_param{2}; 
Mw_vector = adia_param{3}; 
T0 = adia_param{4}; 
R0 = adia_param{5};  
omega = adia_param{6};  
Eu = adia_param{7};  
We = adia_param{8};  
Rel = adia_param{9};  
Mal = adia_param{10};  
PA_nd = adia_param{11};  

R_ig = adia_param{12};  
P0 = adia_param{13}; 
n0 = adia_param{14}; 

NASA_table1 = adia_param{15}; 
NASA_table2 = adia_param{16}; 
NASA_table3 = adia_param{17}; 


%defining mol vector
mol_vector=[n_H2O;n_OH;n_H;n_H2;n_H2O2;n_O;n_HO2;n_O2;n_Ar];%non-dimensionalized mol column vector
mol_total=sum(mol_vector);
mol_frac=mol_vector./mol_total; %column vector of the mole fractions

% Pg = nRT/V
Pg = mol_total*T_nd/R^3;

Mw_ave=sum(mol_frac.*Mw_vector);


%thermo_param = R_ig T0 P0 n0 R0 NASA_table1 NASA_table2 NASA_table3 
thermo_param = {R_ig,T0,P0,n0,R0,...
    NASA_table1,NASA_table2,NASA_table3};
%get thermo data needed for energy balance
[Cp_nd,Cv_nd,H_nd,~,G_nd]=thermo_data(T_nd,thermo_param);
Cp_nd_bar=sum(mol_frac.*Cp_nd);
Cv_nd_bar=sum(mol_frac.*Cv_nd);
gamma=Cp_nd_bar./Cv_nd_bar;

cg=sqrt(gamma*R_ig*T_nd*T0/Mw_ave); 
Mag=R0*omega/(3*cg);


%get kinetics at reaction temperature
k_nd=get_k(T_nd,kin_param);%non-dimensionalized column vector of rate constants
%modify kinetics by toggleing with knd_mod
k_nd = k_nd.*knd_mod;
%calculate dimensionless concentrations by normalizing mols by concentrations
%defined as in S&S 2.4
c_nd=mol_vector./R.^3./(4./3.*pi);%non-dimensionalized concentrations
c_nd(end)=sum(c_nd); % the last element will be the sum of all the species as 'M' can be anything
%biggamma is the change in concentration of each species from kinetic rate
%laws
biggamma=rxn_kinetics(k_nd,c_nd,kin_param);
%change in number of moles of each species defined as in S&S 2.4
V = R.^3; %4./3.*pi.*
dmol=biggamma(1:end-1).*V;

%from energy balance S&S 2.6. 
% (1-gamma)sum(H_i-C_Pbar)Gamma_i) - 3gammaP_gRdot/T
% pgdotterm1 - (3.*gamma.*Pg.*dR./R)
pgdotterm1=pgdottermone(gamma,H_nd,Cv_nd_bar,Cp_nd_bar,T_nd,biggamma,V,G_nd);


Pgdot=pgdotterm1-3.*gamma.*Pg.*dR./R; 

% need to differentiate Pgdot to evaluate S&S 2.1. Cannot calculate 
% the derivative of all terms analytically. Instead approx some
% with backwards finite differences:
% gamma

x=[t0,t_1,t_2,t_3,t_4,gamma0,gamma_1,gamma_2,gamma_3,gamma_4];
gammadot=fobfd(x);

%can double check this part
Pgddotterm1=(3*gamma*Pg*dR.^2-3.*R.*gamma.*Pgdot.*dR-3.*R.*gammadot.*Pg.*dR)./R.^2;
Pgddotterm2=3*gamma*Pg/y(1);
x=[t0,t_1,t_2,t_3,t_4,pgdotterm10,pgdotterm1_1,pgdotterm1_2,pgdotterm1_3,pgdotterm1_4];
pgdotterm1dot=fobfd(x);

%liquid pressure 2.2
Pl=Pg+Mag.*R.*Pgdot-(1./Eu).*(1./(We.*R)+dR./(Rel.*R));
%derivative of 2.2
Pldotterm1=Mag*dR*Pgdot+dR/(Eu*We*R^2)+dR^2/(Eu*Rel*R^2);
Pldot1=Mag*R*pgdotterm1dot+Mag*R*Pgddotterm1+Pldotterm1;
Pldot2=Mag*R*Pgddotterm2+1/(Eu*Rel*R);

%momentum balance; 2.1
d2R=(Eu*(Pl-1+PA_nd*sin(2*pi*t)+Mal*y(1)*(Pgdot+Pldot1))-3*y(2)^2/2)/(y(1)+Eu*Mal*y(1)*Pldot2);

% T derivative from ideal gas law
dT= Pgdot*R^3/mol_total+3*R^2*Pg*dR/mol_total-sum(dmol)*Pg*R^3/mol_total/mol_total;

dy=[dR;d2R;dmol(1);dT;dmol(2:end)];

%this if statements will update the terms at each iteration which will be
%used in the fobfd method
if t~=t0
    pgdotterm1_4=pgdotterm1_3;
    pgdotterm1_3=pgdotterm1_2;
    pgdotterm1_2=pgdotterm1_1;
    pgdotterm1_1=pgdotterm10;
    pgdotterm10=pgdotterm1;

    gamma_4 = gamma_3;
    gamma_3 = gamma_2;
    gamma_2 = gamma_1;
    gamma_1 = gamma0;
    gamma0 = gamma;
    
    t_4=t_3;
    t_3=t_2;
    t_2=t_1;
    t_1=t0;
    t0=t;
end
end