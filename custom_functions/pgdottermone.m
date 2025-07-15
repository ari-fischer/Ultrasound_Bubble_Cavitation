function pgdotterm1=pgdottermone(gamma,H_nd,Cv_nd_bar,Cp_nd_bar,T_nd,biggamma,V,G_nd)

% PGDOTTERMONE is a function that desribes the first term in equation 2.6
% in https://www.jstor.org/stable/3067224
% 
% Inputs are non-dimensionalized 

%SS derivation
sum_term=sum((H_nd-T_nd.*(Cp_nd_bar)).*biggamma);
pgdotterm1=(1-gamma).*sum_term;