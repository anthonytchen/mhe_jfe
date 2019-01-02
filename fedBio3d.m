function dx = fedBio3d (t, x, p, u)
% Differential equations for a three-dimensional fed-batch bioreactor
% 
% Ref: 
%	K. Versyck, J. Van Impe, Chemical Engineering Commmunications 172 (1999) 107-124
% args:
%	t - time
%	x - state variables
%   p - model parameters
%   u - input, including feed flow rate and other process variables
%
% outputs:
%	dx - RHS of the ODEs
%

% fixed values
Cs_in = 500; % g/L
m = 0.29; % g/g
Y_xs = 0.47; % g/g
mu_max = 2.1;

F_in = u;

x(x<0)=0;

Cs = x(1);
Cx = x(2);
V = x(3);

% assign model parameters
Ks=exp(p(1)); Ki=exp(p(2)); % log transform to make sure estimated p is positive


% calculations
muu = mu_max * Cs / (Ks + Cs + Cs*Cs/Ki); % specific growth rate, Haldane kinetics
sig = muu / Y_xs + m; % substrate consumption rate, a linear law
u_by_V = F_in/V;

% RHS of ODEs
dx = zeros(size(x));
dx(1) =  -sig*Cx + u_by_V*Cs_in - u_by_V*Cs;
dx(2) =  muu*Cx - u_by_V*Cx;
dx(3) = u;

return;