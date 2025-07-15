function [val,isterminal,dir]=adiabaticstoppingevent(t,y)
% ADIABATICSTOPPINGEVENT is the stopping criteria during the adiabatic
% phase. We stop it in the first instant the first derivative becomes
% positive.
val=y(2)>0;
isterminal=1;
dir=1;
end