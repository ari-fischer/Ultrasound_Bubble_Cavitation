clc
clear

% add paths to custom functions
add_paths

%toggle kinetics
%0 is use full network; 1 is no network; 2 is only H2O dissociation
ks_mod_sel = 0;

%toggle phase change
% 0 means no phase change, 1 means normal
tog_PC = 1;

Mw_vector=[18;17;1;2;34;16;33;32;40]./1000; % H2O;OH;H;H2;H2O2;O;HO2;O2;Ar kg/mol

R_ig = 8.314462618;

misc_param = {ks_mod_sel,tog_PC,Mw_vector,R_ig} ;

% set rng to default settings for reproducibility 
rng default;
n_trials=200;
X = lhsdesign(n_trials,5); %5 variables

%span space of PA, freq, R0, rho, sigma
lb = [1.5E5, 200E3,  1/10,   0.7, 0.3];
ub = [10E5,  1100E3, 1/1.5, 1.2,   1];
diffs = ub-lb;

%param matrix
p_Mat = X.*diffs+lb;

H2O_convs = [];
H2O_inits = [];
n0_v = [];
E_rxn_v = [];
work_v = [];
rad_yields = [];

for i=1:n_trials


    omega = p_Mat(i,2);
    PA = p_Mat(i,1);

    R_res = 3/omega;
    %R0 = R_res/2.1;    
    R0 = R_res*p_Mat(i,3);    
    
    %rho, v_l, sigma_l
    SA_param = [p_Mat(i,4),1,p_Mat(i,5)]

    [~,~,H2O_conversion,time_profile,radial_profile,...
        Pg_profile,V_profile,T_profile,mol_matrix_profile,...
        ad_switch_idx,n0,T0, P0,thermo_dat_params,~,~]=...
            code_driver_SA(omega, R0, PA,misc_param,SA_param);

    H2O_inits = [H2O_inits, mol_matrix_profile(ad_switch_idx,1)];
    H2O_convs = [H2O_convs,H2O_conversion(end)];

    % n_H2O,n_OH,n_H,n_H2,n_H2O2,n_O,n_HO2,n_O2, Ar
    rad_yield = sum(mol_matrix_profile(end,[2,3,6,7]));

    rad_yields = [rad_yields,rad_yield];

    thermo_params = {T0 R_ig Mw_vector P0 n0 R0};

    [PdVwork,Urxn,~,~,~,~,~]=...
        AF_16Dec_thermo_analysis(T_profile(ad_switch_idx:end),...
            mol_matrix_profile(ad_switch_idx:end,:),...
            Pg_profile(ad_switch_idx:end),V_profile(ad_switch_idx:end),...
            radial_profile(ad_switch_idx:end),...
            time_profile(ad_switch_idx:end),thermo_params,thermo_dat_params);

    E_rxn_v = [E_rxn_v,Urxn];
    work_v = [work_v,PdVwork];
end

N_A = 6.022E23;
csvwrite("E_rxn.csv",E_rxn_v(:))
csvwrite("H2O_extent.csv",H2O_inits(:).*H2O_convs(:).*N_A)
csvwrite("H2O_conversions.csv",H2O_convs(:))
csvwrite("rad_yields.csv",rad_yields(:))
csvwrite("ad_work.csv",H2O_convs(:))

save('LHS_space')

