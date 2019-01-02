% measurement function of fedBio
%
%   x - states
%   p - model parameters

function z = fedBio_mea (x, p)

idm=[1 3]; % only the 1st and 3rd states are measured
z = x(idm);
