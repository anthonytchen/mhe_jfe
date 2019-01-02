%
% Dynamic model of fedBio with three states
%
% args:
%   x0 - states at current time
%   p - model parameters
%   u - input at current time
%   dt - simulation interval, from t to t+dt
%
% returns:
%   x - current state
%

function x = fedBio_state (x0, p, u, dt)
  
%opt = odeset('NonNegative', [1 2 3]);
opt = [];

x0(x0<0) = 0; % make sure they are non-negative
[tt1, tt2] = ode23(@(t,x) fedBio3d(t,x,p,u), [0 dt], x0, opt);
x = tt2(end, :);
x(x<0) = 0; % make sure they are non-negative



