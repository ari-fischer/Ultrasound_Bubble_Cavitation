function [eta_metric,OH_yield,H2O_conversion,time_profile,radial_profile,...
    Pg_profile,V_profile,T_profile,mol_matrix_profile,...
    adiabatic_switch_idx,n0,T0,P0,thermo_dat_params,kin_param,H2O_0,ie] = ...
            code_driver_SA(omega, R0, PA,misc_param,SA_param)

    % misc_param = {ks_mod_sel,tog_PC,Mw_vector,R_ig};
    

    % Extract parameters from NASA polynomials and kinetic data
    [n_vector, A_vector,beta_exponent_vector,Ea_vector,...
        spec_exponent,rxn_extent,NASA_table1,NASA_table2,...
        NASA_table3,NASA_Tmax] = import_kin_thermo_param();

    %toggle kinetics
    %0 is use full network; 1 is no network; 2 is only H2O dissociation
    %ks_mod_sel = 0;
    ks_mod_sel = misc_param{1};
    R_ig = misc_param{4};

    knd_mod = ones(length(Ea_vector),1);
    if ks_mod_sel == 1
        knd_mod = knd_mod.*0;
    elseif ks_mod_sel == 2
        knd_mod(3:end) = knd_mod(3:end).*0;
    end
    % 0 means no phase change, 1 means normal
    tog_PC = misc_param{2};
    
    %specify parameters 
    
    P0=101325; % bulk pressure 
    T0=298; % bulk temperature
    sigmahat=0.35*tog_PC;
    
    [rho_l,v_l,beta_l,sigma_l]=temp_liq_properties(T0,P0/1000);
    
    rho_l = rho_l*SA_param(1) ;
    v_l = v_l*SA_param(2) ;
    sigma_l = sigma_l*SA_param(3) ;


    params=[omega,R0,PA,T0,rho_l,v_l,beta_l,sigma_l,sigmahat,P0,misc_param{4}];

    misc_param2 = {knd_mod,n_vector, A_vector,beta_exponent_vector,Ea_vector,...
        spec_exponent,rxn_extent,NASA_table1,NASA_table2,...
        NASA_table3,NASA_Tmax};

    %run the bubble dynamics simulation
    [eta_metric,OH_yield,H2O_conversion,time_profile,radial_profile,...
        Pg_profile,V_profile,T_profile,mol_matrix_profile,...
        adiabatic_switch_idx,n0,kin_param,H2O_0,ie]=full_model(params,misc_param,misc_param2);

    thermo_dat_params = {R_ig,T0,P0,n0,R0,NASA_table1,NASA_table2,NASA_table3};

end

