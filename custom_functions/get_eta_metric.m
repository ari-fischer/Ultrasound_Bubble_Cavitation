function eta=get_eta_metric(rad_matrix,radius)
% GET_ETA_METRIC calculates a metric which represents the rate of dissolution
% of radical species from the interior of a bubble into the bulk liquid.
% According to https://doi.org/10.1063/1.1925607, we can find that the rate
% of dissolution is proportional to the concentration of the radical species
% in the bubble, divided by the root of their molecular weight and
% multiplied by the radius of the bubble
% Simplifying, this becomes n_rad/(root(Mw)*R). In this study, we only care
% about the OH radical, so the input would be the rad_matrix containing
% just the OH radical and the radius of the bubble at the end of the
% collapse.

global Mw_vector 

%Molecular weights of radical species
Mw_rad=Mw_vector(2); %OH;H;O;HO2 column vector
rad_rootMw=sum(rad_matrix./(sqrt(Mw_rad)));

%First, at the bubble collapse
eta=rad_rootMw/radius;

%Next, calculate the rebound bubble radius using ideal gas law
% R_rebound=((3*total_mol*R_ig*T0)/(4*P0*pi))^(1/3);
% chi_rebound=rad_rootMw/R_rebound;