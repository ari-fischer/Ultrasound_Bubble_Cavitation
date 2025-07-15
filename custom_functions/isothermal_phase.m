function dy=isothermal_phase(t,y,iso_params)
% y = [R,dR,n_H2O,T,n_OH,n_H,n_H2,n_H2O2,n_O,n_HO2,n_O2]
% ISOTHERMAL_PHASE describes the momentum and energy balances during the
% isothermal phase of the bubble collapse as desribed in the S&S model


Mw_vector = iso_params{1}; R_ig = iso_params{2}; T0 = iso_params{3};
R0 = iso_params{4}; omega = iso_params{5}; nncd = iso_params{6}; tauc = iso_params{7};
n_Ar = iso_params{8}; Eu = iso_params{9}; Mal = iso_params{10}; We = iso_params{11};
Rel = iso_params{12}; PA_nd = iso_params{13}; Psat_nd= iso_params{14};
P0 = iso_params{15}; n0 = iso_params{16}; NASA_table1 = iso_params{17};
NASA_table2 = iso_params{18}; NASA_table3 = iso_params{19};

%% initialize
% access variables from y vector
R=y(1); % non-dimensionalized bubble radius
dR=y(2); % derivative of the non-dimensionalized bubble radius
n_H2O=y(3); % non-dimensionalized number of moles of H2O
%T_nd=y(4); % non-dimensionalized temperature
n_OH=y(5); % non-dimentionalized number of moles of OH
n_H=y(6); % non-dimentionalized number of moles of./ H
n_H2=y(7); % non-dimentionalized number of moles of H2
n_H2O2=y(8); % non-dimentionalized number of moles of H2O2
n_O=y(9); % non-dimentionalized number of moles of O
n_HO2=y(10); % non-dimentionalized number of moles of HO2
n_O2=y(11); % non-dimentionalized number of moles of O2


T_nd = 1;
%non-dimensionalized mol column vector
mol_vector=[n_H2O;n_OH;n_H;n_H2;n_H2O2;n_O;n_HO2;n_O2;n_Ar];
mol_total=sum(mol_vector);
mol_frac=mol_vector./mol_total; %column vector of the mole fractions

thermo_param = {R_ig,T0,P0,n0,R0,...
    NASA_table1,NASA_table2,NASA_table3};

%get heat capacities
[Cp_nd,Cv_nd,~,~,~]=thermo_data(T_nd,thermo_param); % some of the values here will be unused
Cp_nd_bar=sum(mol_frac.*Cp_nd);
Cv_nd_bar=sum(mol_frac.*Cv_nd);
%specific heat ratio
gamma=Cp_nd_bar./Cv_nd_bar; % this gamma will be used in the model

%molecular weight weighted average
Mw_bar=sum(mol_frac.*Mw_vector); 

%gas speed of sound in bubble
cg=sqrt(gamma.*R_ig.*T_nd.*T0./Mw_bar);
%mach number in gas phase
Mag=R0*omega./(3.*cg);

%% calculate variable changes
%change in the moles of water from evaporation
dn_H2O=-tauc./R.*(n_H2O-Psat_nd.*y(1).^3);

%gas pressure
Pg=(n_H2O+nncd)./R.^3;

%change in the gas pressure by differentiating {Pg=(n_H2O+nncd)./R.^3}
dPg=dn_H2O/R^3-3*(n_H2O+nncd)*dR/R^4;
%change in the rate of H2O evaporation by differentiating {dn_H2O}
ddn_H2O = -tauc.*(dn_H2O./R-n_H2O./R./R.*dR-2.*Psat_nd.*R.*dR);

%from differentiating the equation 2.1 (https://www.jstor.org/stable/3067224)
%  d2R = 1./(R./Eu+3.*Mag.*Mal.*(n_H2O+nncd)./R./R+Mal./Eu./Rel).*(-3./2.*dR.^2./Eu ...
%      +(Pg+Mag.*R.*(dn_H2O./R.^3-3.*dR./R.^4.*(n_H2O+nncd)) - 1./Eu.*(1./We./R+dR./R./Rel)-1 ...
%      +PA_nd.*sin(2*pi*t)+Mal.*R.*((dn_H2O./R.^3 - 3.*dR./R.^4.*(n_H2O+nncd)).*(1+Mag.*dR)  ...
%      +Mag.*R.*(ddn_H2O./R.^3-3.*dn_H2O.*dR./R.^4+12.*dR.^2./R.^5.*(n_H2O+nncd)-3.*dR.*dn_H2O./R.^4)  ...
%      +dR./Eu./We./R.^2+dR.^2./Rel./Eu./R.^2)));

%remove redundant -3.*dR.*dn_H2O./R.^4?
d2R = 1./(R./Eu+3.*Mag.*Mal.*(n_H2O+nncd)./R./R+Mal./Eu./Rel).*(-3./2.*dR.^2./Eu ...
     +(Pg+Mag.*R.*(dn_H2O./R.^3-3.*dR./R.^4.*(n_H2O+nncd)) - 1./Eu.*(1./We./R+dR./R./Rel)-1 ...
     +PA_nd.*sin(2*pi*t)+Mal.*R.*((dn_H2O./R.^3 - 3.*dR./R.^4.*(n_H2O+nncd)).*(1+Mag.*dR)  ...
     +Mag.*R.*(ddn_H2O./R.^3-3.*dn_H2O.*dR./R.^4+12.*dR.^2./R.^5.*(n_H2O+nncd))  ...
     +dR./Eu./We./R.^2+dR.^2./Rel./Eu./R.^2)));

%isothermal-- temperature constant
dT=0;

%output variables- water is the only species with change in moles
dy=[dR;d2R;dn_H2O;dT;0.*mol_vector(2:end-1)];
end