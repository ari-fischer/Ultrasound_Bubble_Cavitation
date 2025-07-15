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

%define senstivity analysis space
omega_space=355.*1000;
PA_space=2.*10^5;

% Perturbation amount
factor_ls=linspace(1.00000,1.01,5);

% Sensitivity Parameter of interest
% rho_l, v_l, sigma_l, PA
s_param="PA";
fprintf("Running sensitivity analysis for "+s_param+"\n");

% Minnaert Resonance using approximation by physical review letters
% According to Physical Review Letters, R_res (m)=3/omega (Hz)

R0_space=3.2E-6;

% z is the slope of the sensitivities
% R is the Rsq
% outs is the amount of H2O consumed
[z,r,outs,outs_work,H2Os,ie]=zr_generator_H2O(omega_space,R0_space,PA_space,factor_ls,s_param,...
    misc_param);


z