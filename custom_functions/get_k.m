function k_nd_vector=get_k(T_nd,kin_param)
% GET_K will output a column vector of the non-dimensionalized rate
% constants for every reaction
% Input will be the non-dimensionalized temperature at each point during
% the reaction phase
%
% This reaction rate constant is obtained using a modified Arrhenius
% equation k=A*T^beta*exp(-Ea/RT)

T0 = kin_param{1};
R_ig = kin_param{2};
A_vector = kin_param{3};
beta_exponent_vector = kin_param{4};
Ea_vector = kin_param{5};
n_vector = kin_param{6};
omega = kin_param{7};
R0 = kin_param{8};
n0 = kin_param{9};

T_d=T_nd*T0;
% Modified Arrhenius Equation: k=A*T^beta*exp(-Ea/RT)
k_vector=A_vector.*(T_d.^beta_exponent_vector).*exp(-Ea_vector./(R_ig.*T_d)); %cm^(3n-3)/mol^(n-1) s
k_nd_vector=k_vector.*n0.^(n_vector-1)./((100.*R0).^(3.*n_vector-3).*omega);
end




