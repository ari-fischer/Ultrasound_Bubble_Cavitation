function [u_rxn_int,delu_max,delu_min] = get_rxn_energy(T_profile,Mw_vector,kin_param,T0,thermo_dat_params,...
    mol_matrix_profile,n0,V_profile,R0,omega,ad_switch_idx,time_profile)

%initialize the heat capacities
Cv_nd=zeros(length(T_profile),length(Mw_vector));
Cv_d=zeros(length(T_profile),length(Mw_vector));
H_nd=Cv_d;
U_d=Cv_d;


R_ig = thermo_dat_params{1};

% get rate constants and reaction energy for each elementary step
k_nds = zeros(length(T_profile),length(kin_param{4}));
u_rxns = zeros(length(T_profile),length(kin_param{4}));

rxn_extent = kin_param{11};
n_vector = kin_param{6};

for i=1:length(T_profile)
    [~,Cv_nd,Hnd_val,~,~]=thermo_data(T_profile(i)./T0,thermo_dat_params);
    H_nd(i,:)=Hnd_val';
    Cv_d(i,:)=Cv_nd'.*R_ig;
    U_d(i,:)=Hnd_val'.*T0.*R_ig-R_ig*T_profile(i);

    u_rxns(i,:)=rxn_extent*U_d(i,:)';
    k_nds(i,:)=get_k(T_profile(i)./T0,kin_param);
end

%%
% calculate the reaction energy for each step
%k * C_A * C_B * delU

%calculate dimensionless concentrations by normalizing mols by concentrations
%defined as in S&S 2.4
c_nd=mol_matrix_profile./n0./V_profile.*(R0.^3);%non-dimensionalized concentrations
c_nd(:,end)=sum(c_nd,2); % the last element will be the sum of all the species as 'M' can be anything
%get rates of each reaction step
rates = zeros(length(T_profile),length(kin_param{4}));
for i=1:length(T_profile)
   rates(i,:)=get_rates(k_nds(i,:)',c_nd(i,:)',kin_param)';
end
% get the delta U rxn for each step
u_rxns_t = u_rxns.*rates.*n0.*V_profile./(4./3.*pi.*R0.^3).*omega;

%initialize u_rxn_int
u_rxn_int = zeros(1,length(kin_param{4}));

for i=1:length(kin_param{4})
    u_rxn_int(i) = trapz(time_profile(ad_switch_idx:end),u_rxns_t(ad_switch_idx:end,i));
end

sum(u_rxn_int);
delu_min = min(u_rxns(:,1));
delu_max = max(u_rxns(:,1));

end