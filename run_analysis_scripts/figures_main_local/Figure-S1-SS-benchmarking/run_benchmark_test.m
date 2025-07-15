clc
clear

%% run different tests
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
test_ind = 2;


R0s = [4.5;4.5;4.5;4.5;4.5;...
    2.0; 3.0; 6.55; 10.0;...
    4.5; 4.5].*10^-6; % initial radius (um)

omegas = [26.5; 26.5; 26.5; 26.5; 26.5;...
    26.5; 26.5; 26.5; 26.5;...
    38.3; 58.9].*10^3; % ultrasound frequency s-1

PAs = [1.1;1.15;1.2;1.25;1.30;...
    1.2;1.2;1.2;1.2;...
    1.2;1.2].*10^5; % pressure amplitude (bar)

Ts = [];
OHs = [];
Ars = [];

for i = 1:length(R0s)

    omega = omegas(i);
    R0 = R0s(i);
    PA = PAs(i);

    %rho, v_l, sigma_l
    SA_param = [1,1,1]
    
    
    [eta_metric,OH_yield,H2O_conversion,time_profile,radial_profile,...
    Pg_profile,V_profile,T_profile,mol_matrix_profile,...
    ad_switch_idx,n0,T0, P0,H2O_0,thermo_dat_params,kin_param,ie]=...
    code_driver_SA(omega, R0, PA,misc_param,SA_param);

    %get mass fraction Ar
    mol_matrix_profile_nd = mol_matrix_profile./n0;
    masses_out = mol_matrix_profile_nd(end,:).*Mw_vector';
    mass_frac_Ar = masses_out(end)/sum(masses_out);
    
    Ts = [Ts;max(T_profile)];
    OHs = [OHs;max(mol_matrix_profile(:,2))*100/n0];
    Ars = [Ars;mass_frac_Ar*100]

end

fprintf('Report:  ')

out=[[Ts(1:7);Ts(3);Ts(8:9);Ts(3);Ts(10:end)],...
    [OHs(1:7);OHs(3);OHs(8:9);OHs(3);OHs(10:end)],...
    [Ars(1:7);Ars(3);Ars(8:9);Ars(3);Ars(10:end)]]


save('benchmark_out.mat')



