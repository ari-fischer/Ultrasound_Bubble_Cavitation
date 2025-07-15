function rate=get_rates(k_nd,c_nd,kin_param) 
% RXN_KINETICS is a function desribing the reaction kinetics in the S&S
% model. The input values are non-dimensionalized rate constants and
% concentrations, and the output will be the change in concentration. All
% inputs and outputs should be column vectors
spec_exponent = kin_param{10}; 
rxn_extent = kin_param{11};

c_row=c_nd'; % change it to row vector
[m,n]=size(spec_exponent);
mc_row=zeros(m,n);%preallocation
for row=1:m
    mc_row(row,:)=c_row; %change all rows to c_row
end
mc_row_power=mc_row.^spec_exponent;% species involved raise to their stoichiometric coefficient
mc_row_power_prod=prod(mc_row_power,2);%product sum of each row
rate=k_nd.*mc_row_power_prod;%column vector of all the reaction rates

end





