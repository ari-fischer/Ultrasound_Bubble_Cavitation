function [z,r,outs,outs_work,H2Os,ie_M]=zr_generator_H2O(omega_space,R0_space,PA_space,factor_ls,s_param,...
        misc_param)

% ZR_GENERATOR generates the sensitivities (z) and the corresponding R2
% scores (r). Inputs the frequency, R0 and PA space of interest as well as
% some of the bulk liquid properties. Factor_ls is the amount of
% perturbation we will give the specific parameter. It is important to note
% that s_param in this code must either be the strings "rho_l", "v_l",
% "sigma_l" or "PA" which are the liquid density, liquid viscosity, 
% surface tension or pressure ampltitude respectively.
%
% The rows of z and r correspond to the frequency, while the columns
% correspond to the pressure ampltitude.

% ie specifies the time that the stopping event was executed in the code
% use this to check if the simulation reached a collapse endpoint

z=zeros(length(omega_space),length(PA_space)); %this will be the eta gradients
r=z; % this will be the R2 scores of the straight lines
outs=z;
outs_work=z;
H2Os=z;
ie_M=z;


for i=1:length(PA_space)
    for j=1:length(omega_space)
        H2O_cons=zeros(1,length(factor_ls));
        dPg_work=zeros(1,length(factor_ls));
        H2Os_0=zeros(1,length(factor_ls));
        for k=1:length(factor_ls)
            %[factor_ls(k),1,1] -> rho_l, v_l, sigma_l
            if s_param=="rho_l"
                [~,OH_yield,X_H2O,~,~,...
                 Pg_profile,~,~,~,...
                 ~,n0,~,~,~,~,H2O_0,ie]=code_driver_SA(omega_space(j),R0_space(j),PA_space(i),...
                    misc_param,...
                    [factor_ls(k),1,1]);
            elseif s_param=="v_l"
                1;
                [~,OH_yield,X_H2O,~,~,...
                 Pg_profile,~,~,~,...
                 ~,n0,~,~,~,~,H2O_0,ie]=code_driver_SA(omega_space(j),R0_space(j),PA_space(i),...
                    misc_param,...
                    [1,factor_ls(k),1]);
            elseif s_param=="sigma_l"
                [~,OH_yield,X_H2O,~,~,...
                 Pg_profile,~,~,~,...
                 ~,n0,~,~,~,~,H2O_0,ie]=code_driver_SA(omega_space(j),R0_space(j),PA_space(i),...
                    misc_param,...
                    [1,1,factor_ls(k)]);
            elseif s_param=="PA"
                [~,OH_yield,X_H2O,~,~,...
                 Pg_profile,~,~,~,...
                 ~,n0,~,~,~,~,H2O_0,ie]=code_driver_SA(omega_space(j),R0_space(j),PA_space(i)*factor_ls(k),...
                    misc_param,...
                    [1,1,1]);
            end

            
            %get Pg work
            H2O_cons(k) = X_H2O(end)*H2O_0;
            dPg_work(k)=trapz(Pg_profile,Pg_profile.*0+1);

            H2Os_0(k) = H2O_0;
        end

        % turning the values into % changes
        H2O_change=H2O_cons./H2O_cons(1);%.*100-100;
        factor=factor_ls;%-100;

        % add these codes in to inspect linearity
        %figure
        %plot(factor_ls,OH_yields_change,'o')   

        %fitting the linear model and getting the gradient and R2
        %scores
        mdl=fitlm(factor,H2O_change);
        z(j,i)=mdl.Coefficients.Estimate(2);
        r(j,i)=mdl.Rsquared.Ordinary;
        %reportes output values for the unperturbed collapse
        outs(j,i) = H2O_cons(3);
        outs_work(j,i) = dPg_work(3);
        H2Os(j,i) = H2Os_0(3);
        ie_M(j,i) = ie;
    end
end
            