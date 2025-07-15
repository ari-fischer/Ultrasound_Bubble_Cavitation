function [n_vector, A_vector,beta_exponent_vector,Ea_vector,...
    spec_exponent,rxn_extent,NASA_table1,NASA_table2,NASA_table3,NASA_Tmax] ...
    = import_kin_thermo_param()

%read data from NASA polynomials
NASA_table1=readtable("WaterArgon200to1000.csv");
NASA_table2=readtable("WaterArgon1000to6000.csv");
NASA_table3=readtable("WaterArgon6000to20000.csv");
NASA_Tmax=readtable("WaterArgonT_max.csv");

%read data from kinetics parameterization
rate_constant_table=readtable("hydrogen_rxn_rate_constant.csv");
spec_exponent_table=readtable('spec_exponent.csv');
rxn_extent_table=readtable('hydrogen_rxn_extent.csv');

%% convert the data into useful vectors and matrices
% total stoichiometry of reactants in that step
n_vector=table2array(rate_constant_table(:,3));
%preexponential factor
A_vector=table2array(rate_constant_table(:,4));
beta_exponent_vector=table2array(rate_constant_table(:,5));
%activation energies
Ea_vector_cal=table2array(rate_constant_table(:,6));% units in cal/mol, need to convert to J/mol
Ea_vector=Ea_vector_cal./1000.*4184;%convert into J/mol
%stoichiometric coefficient for each species in that step
spec_exponent=table2array(spec_exponent_table(:,3:end));
%extent of reaction for each species in that step
rxn_extent=table2array(rxn_extent_table(:,3:end));
end