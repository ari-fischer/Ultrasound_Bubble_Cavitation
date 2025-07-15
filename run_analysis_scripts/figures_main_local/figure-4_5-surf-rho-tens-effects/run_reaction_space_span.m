clc
clear
restoredefaultpath

add_paths

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
test_ind = 1;

M_size = 30;
omega=355.*1000;
PA=2*1E5;

factor_l = linspace(0.7,1.3,M_size);
%factor_l = linspace(20/72,1,M_size);

% rho, v, sig
param = 1;
%param = 3;

H2O_convs_M = [];
H2O_inits_M = [];
work_done_M = [];

n0_M = [];
E_rxn_M = [];

E_therm_M = [];
ie_M = [];

N_Ar_M = [];


tic

for i=1:length(factor_l)
    work_done_v = [];
    H2O_convs = [];
    H2O_inits = [];

    n0_v = [];
    E_rxn_v = [];

    E_therm_v = [];

    for j=1:length(omega)
        
        R0=3.2E-6; 
        
        % rho, v, sig
        SA_param=[1,1,1];
        SA_param(param) = factor_l(i);

       
        [eta_metric,OH_yield,H2O_conversion,time_profile,radial_profile,...
        Pg_profile,V_profile,T_profile,mol_matrix_profile,...
        ad_switch_idx,n0,T0, P0,thermo_dat_params,kin_param,H2O_0,ie]=...
            code_driver_SA(omega, R0, PA,misc_param,SA_param);

        H2O_inits = [H2O_inits, mol_matrix_profile(ad_switch_idx,1)];
        H2O_convs = [H2O_convs,H2O_conversion(end)];
        
        work_done = -trapz(V_profile(ad_switch_idx:end),Pg_profile(ad_switch_idx:end));
        work_done_v = [work_done_v,work_done];

        n0_v = [n0_v,n0];

        ind_max = find(radial_profile==max(radial_profile));
        total_work = -trapz(V_profile(ind_max:end),Pg_profile(ind_max:end));
        
        thermo_params = {T0 R_ig Mw_vector P0 n0 R0};
        
        [PdV_work,Urxn,Uheat,Ustate,UH2O,UAr,U_r1]=...
            AF_16Dec_thermo_analysis(T_profile(ad_switch_idx:end),...
                mol_matrix_profile(ad_switch_idx:end,:),...
                Pg_profile(ad_switch_idx:end),V_profile(ad_switch_idx:end),...
                radial_profile(ad_switch_idx:end),...
                time_profile(ad_switch_idx:end),thermo_params,thermo_dat_params);

        E_therm_v = [E_therm_v,Uheat];
        %% extract rxn energy
        u_rxn_int = get_rxn_energy(T_profile,Mw_vector,kin_param,T0,thermo_dat_params,...
    mol_matrix_profile,n0,V_profile,R0,omega,ad_switch_idx,time_profile);
        
        E_rxn_v = [E_rxn_v,Urxn];
  
        N_Ar_M = [N_Ar_M;mol_matrix_profile(1,end)];

        ie_M = [ie_M,ie]
    end
    H2O_convs_M = [H2O_convs_M; H2O_convs];
    H2O_inits_M = [H2O_inits_M; H2O_inits];
    work_done_M = [work_done_M; work_done_v];
    n0_M = [n0_M;n0_v];
    E_rxn_M = [E_rxn_M; E_rxn_v];

    E_therm_M = [E_therm_M;E_therm_v];
end

toc

H2O_cons_M = H2O_convs_M.*H2O_inits_M;
N_A = 6.022E23;
H2O_cons_M=H2O_cons_M*N_A;


save(['figure-5-plot-355kHz-density'])
