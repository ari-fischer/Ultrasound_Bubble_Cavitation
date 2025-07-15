function [eta_metric,OH_yield,H2O_conversion,time_profile,radial_profile,...
    Pg_profile,V_profile,T_profile,mol_matrix_profile,...
    adiabatic_switch_idx,n0,kin_param,H2O_0,ie]...
    =full_model(params,misc_param,misc_param2)

% FULL_MODEL Computes various outputs that may be of interest to
% user. They are:
% 1) eta_metric: is the metric computed in get_eta_metric. 
% 2) OH radical yield as a result of the bubble collapse
% 3) H2O conversion due to the bubble collapse
% 4) the time, radial, gas pressure, volume, temperature and mol profiles
% of the model
% 5) the index at which the adiabatic switch occurs, which can be used for
% plotting the profiles, and used for obtaining the beginning of the
% adiabatic phase, where some thermodynamic properties may be needed in
% thermo_analysis.
%
% Requires a row vector known as params as the input. The following is the
% order of the inputs required in params:
% 1) omega (the frequency in Hz)
% 2) R0 (initial bubble radius in m)
% 3) PA (Pressure Amplitude in Pa)
% 4) T0 (Initial Temperature in K)
% 5) rho_l (Bulk Liquid Density in kg/m3)
% 6) v_l (Bulk liquid viscosity in m2/s)
% 7) beta_l (Liquid Bulk Modulus in Pa)
% 8) sigma_l (Surface tension in N/m)
% 9) sigmahat (accommodation coefficient, dimensionless)
%
% Ensure all global variables required are imported in the script before
% using this function. These are:
% 1) NASA Tables
% 2) Reaction Kinetic tables and tables derived from them
% 3) Molecular Weights
% 4) Ideal gas constant
% 5) Thermal and Diffusion Coefficients

% See also GET_ETA_METRIC, THERMO_ANALYSIS

    % misc_param = {ks_mod_sel,tog_PC,Mw_vector};

%    misc_params2 = {knd_mod,n_vector, A_vector,beta_exponent_vector,Ea_vector,...
%        spec_exponent,rxn_extent,NASA_table1,NASA_table2,...
%        NASA_table3,NASA_Tmax};


      
global t0 t_1 t_2 t_3 t_4 gamma_4 gamma_3 gamma_2 gamma_1 gamma0 ...
    pgdotterm10 pgdotterm1_1 pgdotterm1_2 pgdotterm1_3 pgdotterm1_4

% Defining the variable in params
omega=params(1); R0=params(2); PA=params(3); T0=params(4);
rho_l=params(5); v_l=params(6); beta_l=params(7); sigma_l=params(8);
sigmahat=params(9); P0=params(10); R_ig=params(11);

A_vec = misc_param2{3}; beta_vec = misc_param2{4}; Ea_vec = misc_param2{5};
n_vec = misc_param2{2}; spec_exponent = misc_param2{6}; 
rxn_extent = misc_param2{7};


Mw_vector = misc_param{3};

Mw_H2O=Mw_vector(1); % Molecular weight of H2O
Mw_Ar=Mw_vector(end); % mol weight Ar

%% dimensionless quantities that go into bubble model
%speed of sound liquid state
cl=sqrt(beta_l/rho_l);

% Non-dimensionalized Saturated Pressure of Water at T0, which will be
% 298K, requires antoine equation parameters
Psat=10^(6.20963-2354.731/(T0+7.559))*10^5;%antoine eqeuation paramters taken from NIST
Psat_nd=Psat/P0;

%non-dimensionalized pressure amplitude
PA_nd=PA/P0;

%Dimensionless Numbers, note some of the values needed will be defined at a later time
Eu=P0./(omega.^2.*R0.^2.*rho_l); %Euler's number
We=(rho_l*R0^3*omega^2)/(2*sigma_l); %Weber number
Mal=R0*omega/cl; %liquid mach number
Rel=(omega*R0^2)/(4*v_l); %Reynold number

%Information on the number of moles
n0=P0*(4/3*pi*R0^3)/R_ig/T0;%non-dimensionalization is done using this value
n_init=(P0+2*sigma_l/R0)*(4/3*pi*R0^3)/R_ig/T0;%initial total number of moles in the bubble
n_H2O_init=0;%initially a dry bubble
nnc=n_init-n_H2O_init;%number of moles of non-condensible gas
nncd=nnc/n0;%non-dimensionalized number of moles of non-condensible gas
n_H2O_init_nd=n_H2O_init/n0;%non-dimensionalized initial number of moles of water vapor
tauc=(3*sigmahat/(R0.*omega))*(R_ig*T0/(2*pi*Mw_H2O))^0.5;
n_Ar=1*nncd;%for now the whole bubble will be argon, at a later time, it could be eg 90% Ar etc


%% run the isothermal simulation
%specify a long time range that encompasses collapse
t_nd=50*10^-6*omega;%non-dimensional time we will try to simulate the isothermal phase until
tspan1=linspace(0,t_nd,1000);%time points

%initial values for R,dR,n_H2O,T,n_OH,n_H,n_H2,n_H2O2,n_O,n_HO2,n_O2
y0=[1;0;n_H2O_init_nd;1;0;0;0;0;0;0;0]; 

%pass important parameters for isothermal model
iso_params = {Mw_vector R_ig T0 R0 omega nncd tauc n_Ar Eu ...
    Mal We Rel PA_nd Psat_nd P0 n0 ... 
    misc_param2{8} misc_param2{9} misc_param2{10}};

%integrate the isothermal model
options=odeset('Events',@(t,y) isothermalstoppingevent(t,y,[omega,R0]),'RelTol',10^-6,'AbsTol',10^-6);
[t1,y1]=ode23(@(t,y) isothermal_phase(t,y,iso_params),tspan1,y0,options);


% run simulation to the adiabatic stop condition to get approximate tf
options=odeset('Events',@adiabaticstoppingevent,'RelTol',10^-6,'AbsTol',10^-6);
[t1f,~]=ode23(@(t,y) isothermal_phase(t,y,iso_params),tspan1,y0,options);

t_fin = (t1f(end)-t1(end))+t1f(end);

%% run adiabatic simulation


thermo_param = {R_ig,T0,P0,n0,R0,...
    misc_param2{8}, misc_param2{9}, misc_param2{10}};


%initialize BFD method 
%to calculate gamma, i need the last 5 values of the isothermal phase
y1last5=y1(end-4:end,[3,5:11]);
y1last5(:,end+1)=n_Ar;
y1last5molfrac=y1last5./sum(y1last5,2);
[Cp_nd,Cv_nd,~,~,~]=thermo_data(1,thermo_param);
Cp_nd_bar=y1last5molfrac*Cp_nd;
Cv_nd_bar=y1last5molfrac*Cv_nd;
gamma_vector=Cp_nd_bar./Cv_nd_bar;

%get a time-span that with surpasse the collapse
tinitial=t1(end);
tspan2=linspace(tinitial,t_fin,20000);
y0=y1(end,:); %updating the initial values for the adiabatic phase

% for the fourth order backward finite difference in the adiabatic phase
t0=t1(end,1);
t_1=t1(end-1,1);
t_2=t1(end-2,1);
t_3=t1(end-3,1);
t_4=t1(end-4,1);
pgdotterm10=0;
pgdotterm1_1=0;
pgdotterm1_2=0;
pgdotterm1_3=0;
pgdotterm1_4=0;
gamma0=gamma_vector(end);
gamma_1=gamma_vector(end-1);
gamma_2=gamma_vector(end-2);
gamma_3=gamma_vector(end-3);
gamma_4=gamma_vector(end-4);


%    misc_params2 = {knd_mod,n_vector, A_vector,beta_exponent_vector,Ea_vector,...
%        spec_exponent,rxn_extent,NASA_table1,NASA_table2,...
%        NASA_table3,NASA_Tmax};

% 1 is ks_mod
adia_param = {misc_param2{1},n_Ar, Mw_vector, T0, R0, omega,...
    Eu, We, Rel, Mal, PA_nd,R_ig,P0,n0,...
    misc_param2{8},misc_param2{9},misc_param2{10}};



kin_param = {T0, R_ig, A_vec, beta_vec, Ea_vec,...
    n_vec, omega, R0, n0, spec_exponent, rxn_extent};

% integrate forward using the adiabatic model

%testing with no reaction to get energy balance correct.
options=odeset('Events',@adiabaticstoppingevent,'RelTol',10^-6,'AbsTol',10^-6,'NonNegative',[3,5:11]);
[t2,y,ie]=ode23(@(t,y) adiabatic_phase(t,y,adia_param,kin_param),tspan2,y0,options);
% ie returns the time that the stopping event was executed

%% collected output data
if length(ie)~=1
    ie = 0;
end
OH_yield=y(end,5);% OH yield at the end of the adiabatic phase %nd

%extracting important data profiles
time_profile_nd=zeros(length(t1)+length(t2)-1,1);
time_profile_nd(1:length(t1))=t1;
time_profile_nd(length(t1)+1:end)=t2(2:end);
time_profile=time_profile_nd./omega;

radial_profile_nd=time_profile_nd;
radial_profile_nd(1:length(t1))=y1(:,1);
radial_profile_nd(length(t1)+1:end)=y(2:end,1);
radial_profile=radial_profile_nd.*R0;

V_profile=(4/3).*pi.*radial_profile.^3;

% R,dR,n_H2O,T,n_OH,n_H,n_H2,n_H2O2,n_O,n_HO2,n_O2
mol_matrix_profile_nd=zeros(length(t1)+length(t2)-1,length(Mw_vector)-1);
mol_matrix_profile_nd(1:length(t1),:)=y1(:,[3,5:end]);
mol_matrix_profile_nd(length(t1)+1:end,:)=y(2:end,[3,5:end]);
mol_matrix_profile_nd(:,end+1)=n_Ar;
mol_matrix_profile=mol_matrix_profile_nd.*n0;

T_profile_nd=time_profile_nd;
T_profile_nd(1:length(t1))=y1(:,4);
T_profile_nd(length(t1)+1:end)=y(2:end,4);
T_profile=T_profile_nd.*T0;

%from ideal gas EOS
% P = nRT/V
Pg_profile=sum(mol_matrix_profile,2).*R_ig.*T_profile...
    ./V_profile;

adiabatic_switch_idx=length(t1);

%final values
H2O_conversion=(y(1,3)-y(:,3))/y(1,3);
H2O_0 = y(1,3);
%radical flux
%rad_matrix=y(end,5)*n0; % just the OH radical
%radius=y(end,1)*R0;
%[eta_metric]=get_eta_metric(rad_matrix,radius);
eta_metric = 1;
end