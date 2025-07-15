function [PdVwork,Urxn_int,Uheat_int,Ustate,Urxn,UAr,UH2O]=...
            AF_16Dec_thermo_analysis(T,mol,Pg,V,R,t,thermo_params,thermo_params2)
% THERMO_ANALYSIS takes in the temperature, mol, gas pressure and volume
% profiles of the bubble during the adiabatic phase. Various thermodynamic
% properties are then calculated as the output. These thermodynamic
% properties are of the system, not per molar values.
% H2O, OH, H, H2, H2O2, O, HO2, O2, Ar


%% initialize
T0 = thermo_params{1};
R_ig = thermo_params{2};
Mw_vector = thermo_params{3};
P0 = thermo_params{4};
n0 = thermo_params{5};
R0= thermo_params{6};

Tnd=T./T0;

mol_frac = mol./sum(mol,2);


%% state functions

%% Enthalpy
% Reference state is the initial state during the adiabatic phase
% Initial State
[~,~,Hnd_1,~,Gnd_1]=thermo_data(Tnd(1),thermo_params2);
Hd_1=Hnd_1.*T0.*R_ig;% units of J/mol
Ud_1=Hd_1-R_ig*T(1);
H_1=mol(1,:)*Hd_1;
U_1=mol(1,:)*Ud_1;

% 2nd state after reaction completes at the temperature of initial state
[~,~,Hnd_2,~,Gnd_2]=thermo_data(Tnd(end),thermo_params2);
Hd_2=Hnd_2.*T0.*R_ig;% units of J/mol
Ud_2=Hd_2-R_ig*T(end);
H_2=mol(1,:)*Hd_2;
U_2=mol(1,:)*Ud_2;

% Final state is heating up all the products to the final temperature
[~,~,Hnd_3,~,Gnd_3]=thermo_data(Tnd(end),thermo_params2);
Hd_3=Hnd_3.*T0.*R_ig;% units of J/mol;
Ud_3=Hd_3-R_ig*T(end);
H_3=mol(end,:)*Hd_3;
U_3=mol(end,:)*Ud_3;

Hrxn=H_2-H_1;
Hheat=H_3-H_2;
Hstate=H_3-H_1;

VdPwork=trapz(Pg,V);
%% Internal Energy
% Make use of the relationship, U=H-PV
% Reference state is the initial state during the adiabatic phase
% Initial State



Urxn=U_3-U_2;
Uheat=U_2-U_1;
Ustate=U_3-U_1;

PdVwork=-trapz(V,Pg);

UH2O_1 = mol(1,1).*(Hd_1(1)-R_ig*T(1));
UH2O_2 = mol(1,1).*(Hd_2(1)-R_ig*T(end));

UH2O =  UH2O_2-UH2O_1;

UAr_1 = mol(1,end).*(Hd_1(end)-R_ig*T(1));
UAr_2 = mol(1,end).*(Hd_2(end)-R_ig*T(end));
UAr =  UAr_2-UAr_1;
1;
%    ...
%    *Hd_2(end)*mol(1,1)*R_ig*T(end)-mol(1,1)*Hd_1(1)*mol(1,1)*R_ig*T(1));

%energy going to just driving the H2O dissociation reaction

H_nd2=zeros(length(Tnd),length(Mw_vector));
Cv_nds=zeros(length(Tnd),length(Mw_vector));
Cp_nds=zeros(length(Tnd),length(Mw_vector));
for i=1:length(Tnd)
    [Cp_nd,Cv_nd,Hnd_val,~,~]=thermo_data(Tnd(i),thermo_params2);
    H_nd2(i,:)=Hnd_val';
    Cv_nds(i,:)=Cv_nd';
    Cp_nds(i,:)=Cp_nd';
end

H_d2=H_nd2.*T0.*R_ig; % J/mol

U_d2 = H_d2-R_ig.*T;

%integrate heating for all species

Uheat_int = 0;
for i=1:length(Mw_vector)
    Uheat_int = Uheat_int + ...
        trapz(T,mol(:,i).*Cv_nds(:,i).*R_ig);
end

Urxn_int = 0;
for i=1:length(Mw_vector)
    Urxn_int = Urxn_int + ...
        trapz(mol(:,i),H_d2(:,i)-R_ig.*T);
end


"Uheat_int/Uheat"
Uheat_int/Uheat
"delU/PdVwork"
(Uheat_int+Urxn_int)/PdVwork
"Ustate/PdVwork"
Ustate/PdVwork

1;
end





