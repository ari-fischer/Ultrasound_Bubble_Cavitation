function [OH_yields,OH_yield] = get_OH_yield(T_profile,Mw_vector,kin_param,T0,thermo_dat_params,...
    mol_matrix_profile,n0,V_profile,R0,omega,ad_switch_idx,time_profile)

%initialize the heat capacities

R_ig = thermo_dat_params{1};

% get rate constants and reaction energy for each elementary step
k_nds = zeros(length(T_profile),length(kin_param{4}));

rxn_extent = kin_param{11};
n_vector = kin_param{6};

for i=1:length(T_profile)
    k_nds(i,:)=get_k(T_profile(i)./T0,kin_param);
end

%%
% calculate the reaction energy for each step
%k * C_A * C_B 

%calculate dimensionless concentrations by normalizing mols by concentrations
%defined as in S&S 2.4
c_nd=mol_matrix_profile./n0./V_profile.*(R0.^3);%non-dimensionalized concentrations
c_nd(:,end)=sum(c_nd,2); % the last element will be the sum of all the species as 'M' can be anything
%get rates of each reaction step
rates = zeros(length(T_profile),length(kin_param{4}));
for i=1:length(T_profile)
   rates(i,:)=get_rates(k_nds(i,:)',c_nd(i,:)',kin_param)';
end

rates = rates.*n0.*V_profile./(4./3.*pi.*R0.^3).*omega;
%integrate rates that depend on OH ()
for i=1:length(kin_param{4})
    rates_int(i) = trapz(time_profile(ad_switch_idx:end),...
        rates(ad_switch_idx:end,i));
end

OH_yields = rxn_extent(:,2).*rates_int';
OH_yield = sum(OH_yields);

end