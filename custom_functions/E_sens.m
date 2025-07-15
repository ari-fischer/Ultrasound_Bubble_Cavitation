function [S_work, S_therm, S_rxn]=E_sens(omega_space,R0_space,PA_space,factor_ls,s_param,...
        misc_param)

z=zeros(length(omega_space),length(PA_space)); %this will be the eta gradients
r=z; % this will be the R2 scores of the straight lines
outs=z;
outs_work=z;
H2Os=z;
ie_M=z;


for i=1:length(PA_space)
    for j=1:length(omega_space)
        works=zeros(1,length(factor_ls));
        E_therms = works;
        E_rxns = works;
        for k=1:length(factor_ls)

        if s_param == "PA"
            % turning the values into % changes
            factor=factor_ls*PA_space(i);%-100;  
    
                [eta_metric,OH_yield,H2O_conversion,time_profile,radial_profile,...
            Pg_profile,V_profile,T_profile,mol_matrix_profile,...
            ad_switch_idx,n0,T0, P0,thermo_dat_params,kin_param,H2O_0,ie]=code_driver_SA(omega_space(j),R0_space(j),PA_space(i)*factor_ls(k),...
                    misc_param,...
                    [1,1,1]);
        elseif s_param == "rho_l"
            factor=factor_ls*997;  
    
                [eta_metric,OH_yield,H2O_conversion,time_profile,radial_profile,...
            Pg_profile,V_profile,T_profile,mol_matrix_profile,...
            ad_switch_idx,n0,T0, P0,thermo_dat_params,kin_param,H2O_0,ie]=code_driver_SA(omega_space(j),R0_space(j),PA_space(i),...
                    misc_param,...
                    [1*factor_ls(k),1,1]);
        elseif s_param == "sigma_l"
            factor=factor_ls*72/1000; 
    
                [eta_metric,OH_yield,H2O_conversion,time_profile,radial_profile,...
            Pg_profile,V_profile,T_profile,mol_matrix_profile,...
            ad_switch_idx,n0,T0, P0,thermo_dat_params,kin_param,H2O_0,ie]=code_driver_SA(omega_space(j),R0_space(j),PA_space(i),...
                    misc_param,...
                    [1,1,1*factor_ls(k)]);
        end

        thermo_params = {T0 misc_param{4} misc_param{3} P0 n0 R0_space};

        [PdV_work,Urxn,Uheat,Ustate,UH2O,UAr,U_r1]=...
            AF_16Dec_thermo_analysis(T_profile(ad_switch_idx:end),...
            mol_matrix_profile(ad_switch_idx:end,:),...
            Pg_profile(ad_switch_idx:end),V_profile(ad_switch_idx:end),...
            radial_profile(ad_switch_idx:end),...
            time_profile(ad_switch_idx:end),thermo_params,thermo_dat_params);

            works(k)=-trapz(V_profile(ad_switch_idx:end),Pg_profile(ad_switch_idx:end));
            E_therms(k) = Uheat
            E_rxns(k) = Urxn
        end

        % turning the values into % changes

        %fitting the linear model and getting the gradient and R2
        %scores
        mdl=fitlm(factor,works);
        S_work(j,i)=mdl.Coefficients.Estimate(2);
        
        mdl=fitlm(factor,E_therms);
        S_therm(j,i)=mdl.Coefficients.Estimate(2);

        mdl=fitlm(factor,E_rxns);
        S_rxn(j,i)=mdl.Coefficients.Estimate(2);

    end
end
            