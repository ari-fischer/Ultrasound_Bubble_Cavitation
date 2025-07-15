function [val,isterminal,dir]=isothermalstoppingevent(t,y,params)
% ISOTHERMALSTOPPINGEVENT is the stopping criteria for the isothermal
% phase. Beyond which, the adiabatic phase will take over. This stopping
% event is as described by the relative momentum and transport time scales
% in the S&S model
omega = params(1);
R0 = params(2);
a=1.77*10^-5;% thermal diffusivity of Ar
Dv=0.2570/100/100;% mass diffusivity of Ar-- Fuller equation

%update 11/22. need to add in the dimensions to y(1)
taudiffh=R0.^2.*omega./(4.*a.*y(1)); 
taudiffm=R0.^2.*omega./(4.*Dv.*y(1));
taudiff=(taudiffh+taudiffm)./2;
taudyn=y(1)./abs(y(2));
val=(y(2)<0 && taudyn<taudiff);
isterminal=1;
dir=1;
end