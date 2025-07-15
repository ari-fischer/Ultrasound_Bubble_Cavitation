clc
clear
restoredefaultpath

%% 
% use this to generate the time-dependent radii, temperature, and other ... 
% properties plotted in figure 1
%store the following data 

add_paths()

%Cp toggle : 0 means use NASA, 1 means use Scho
use_NASA_Scho = 0;
%toggle kinetics
%0 is use full network; 1 is no network; 2 is only H2O dissociation
ks_mod_sel = 0;
% 0 means no phase change, 1 means normal
% tog_PC = 1;
tog_PC = 1;

Mw_vector=[18;17;1;2;34;16;33;32;40]./1000; % H2O;OH;H;H2;H2O2;O;HO2;O2;Ar kg/mol

R_ig = 8.314462618;

misc_param = {ks_mod_sel,tog_PC,Mw_vector,R_ig} ;
% test index:
% 1 - benchmark dataset storey and szeri (Table 1)
% SS standard conditions

%omega = 300.*10^3;
%PA = 9.*10^5;
omega = 355.*10^3;
PA = 2.*10^5;


R_res = 3/omega;
R0 = R_res/2.1;
R0 = 3.2*10^-6;

%rho, v_l, sigma_l
SA_param = [1,1,1]


[eta_metric,OH_yield,H2O_conversion,time_profile,radial_profile,...
Pg_profile,V_profile,T_profile,mol_matrix_profile,...
ad_switch_idx,n0,T0, P0,H2O_0,thermo_dat_params,kin_param,ie]=...
    code_driver_SA(omega, R0, PA,misc_param,SA_param);

%% extra energy for heating Ar

Cv_nd=zeros(length(T_profile),length(Mw_vector));
Cv_d=zeros(length(T_profile),length(Mw_vector));
H_nd=Cv_d;

[n_vector, A_vector,beta_exponent_vector,Ea_vector,...
    spec_exponent,rxn_extent,NASA_table1,NASA_table2,...
    NASA_table3,NASA_Tmax] = import_kin_thermo_param();

thermo_param = {thermo_dat_params{[2,1]},P0,n0,R0,...
    NASA_table1,NASA_table2,NASA_table3};

for i=1:length(T_profile)
    [Cp_nd,Cv_nd,Hnd_val,~,~]=thermo_data(T_profile(i)./T0,...
        thermo_param);
    H_nd(i,:)=Hnd_val';
    Cv_d(i,:)=Cv_nd'.*R_ig;
end

%get mass fraction Ar
mol_matrix_profile_nd = mol_matrix_profile./n0;
masses_out = mol_matrix_profile_nd(end,:).*Mw_vector';
mass_frac_Ar = masses_out(end)/sum(masses_out);

%% check thermo next
Ts = max(T_profile)
OHs = max(mol_matrix_profile(:,2))*100/n0;
Ars = mass_frac_Ar*100;

cum_work=zeros(1,length(time_profile(ad_switch_idx:end))-1);

for i=1:length(cum_work)
    cum_work(i) = -trapz(V_profile(ad_switch_idx:ad_switch_idx+i),...
        Pg_profile(ad_switch_idx:ad_switch_idx+i));
end


save('plot_figure1')
