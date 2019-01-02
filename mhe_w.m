function [x, p, f_obj] = mhe_w(xseq, useq, yseq, p, fn_state, fn_obsv, Cov_x, Cov_y, Cov_p, t_sample, bSnP)
% Moving horizon estimation
% Optimisation wrt x_{t-n} and w_{t-n:t}, not the state sequence x_{t-n:t} directly
%
% Inputs:
%   xseq -- current estimate of state sequence, nt by nx; 
%           used as the initial value and its size tells the time horizon for MHE
%   useq -- sequence of inputs, nt by nu
%   yseq -- sequence of measurements, nt by ny
%   p -- model parameters [or their current estimate], np by 1
%   fn_state, fn_obsv -- name of state and observation functions
%       these functions take the following syntax:
%           fn_state(x, p, u, dt)
%           fn_obsv(x, p)
%   Cov_x, Cov_y, Cov_p -- vectors, the diagonal of covariance matrices of 
%                           state, measurement and parameters, respectively
%   t_sample -- sampling interval of the system
%   bSnP -- integer; 
%       1: both states and parameters; 0: only states to be updated; -1: only parameters to be updated
%
% Outputs:
%   x -- estimated state sequence; nt by nx
%   p -- estimated model parameters; np by 1
%   f_obj -- value of objective function

if nargin < 11
    bSnP = 0;
end

% Optimiser's options; use the 2nd line for no display of information
opt = optimset('Display', 'iter', 'MaxFunEvals', 5e4, 'MaxIter', 5e4, 'TolX', 1e-8);    
%opt = optimset('Display', 'off', 'MaxFunEvals', 5e4, 'MaxIter', 5e4, 'TolX', 1e-8);    

[nt, nx] = size(xseq);
np = length(p);

w0 = zeros(nt-1,nx); %initially assume no state function noise

if bSnP == 1 % estimate both state and parameters   

    % organise states and parameters into a single vector for optimisation    
    x_w_p0 = [xseq(1,:)'; reshape(w0, (nt-1)*nx, 1); p'];
    nxwp = length(x_w_p0);    

    x0 = xseq(1,:);
    x0Q = Cov_x; % this is a crude estimate; probably should use EKF to approximate the prior distribution
    p0 = p;
        
    [x_w_p1, f_obj] = fminunc(@(x_w_p) mhe_obj(x_w_p, nx, x0, x0Q, p0, Cov_p, useq, yseq, fn_state, Cov_x, ...
        fn_obsv, Cov_y, t_sample), x_w_p0, opt);
    
    x = xseq;
    x(1,:) = x_w_p1(1:nx)';
    w = reshape(x_w_p1(nx+1:nx*nt), nt-1, nx);
    p = x_w_p1(nx*nt+1 : nxwp);
    for i = 1 : nt-1,
        ui = useq(i, :);
        xx = feval(fn_state, x(i,:), p, ui, t_sample) + w(i,:);
        %xx(xx<0)=0;
        x(i+1,:) = xx;
    end    
    
elseif bSnP == 0 % estimate states only

    % organise states into a single vector for optimisation
    x_w0 = [xseq(1,:)'; reshape(w0, (nt-1)*nx,1)];   

    x0 = xseq(1,:);
    x0Q = Cov_x; % this is a crude estimate; probably should use EKF to approximate the prior distribution

    [x_w1, f_obj] = fminunc(@(x_w) mhe_obj_states(x_w, nx, x0, x0Q, p, useq, yseq, fn_state, Cov_x, ...
        fn_obsv, Cov_y, t_sample), x_w0, opt);
    
    x = xseq;
    x(1,:) = x_w1(1:nx)';
    w = reshape(x_w1(nx+1:nx*nt), nt-1, nx);
    for i = 1 : nt-1,
        ui = useq(i, :);
        xx = feval(fn_state, x(i,:), p, ui, t_sample) + w(i,:);   
        x(i+1,:) = xx;
    end
    
elseif bSnP == -1 % estimate parameters only
 
    % organise states and parameters into a single vector for optimisation    
    p0 = p';
    x0 = xseq(1,:);    
        
    [p1, f_obj] = fminunc(@(p) mhe_obj_paras(p, nx, x0, useq, yseq, fn_state, fn_obsv, Cov_y, t_sample), p0, opt);
    
    x = xseq;    
    p = p1;
    for i = 1 : nt-1,
        ui = useq(i, :);
        xx = feval(fn_state, x(i,:), p, ui, t_sample);
        x(i+1,:) = xx;
    end
        
else % error message
    error('bSnP must be 0, 1, or -1');
end



%% Internal function to calculate objective function to be minimised
%   by adjusting both states and parameters
function f_mhe = mhe_obj(x_w_p, nx, x0, Cov_x0, p0, Cov_p, useq, yseq, fn_state, Cov_x, fn_obsv, Cov_y, t_sample)
% The first nx columns of x_n_p contain the states, whilst the
%   rest (nxp-nx) columns contain the parameters.
%   Parameters are fixed across the time span of nt

[nt, ny] = size(yseq);
nxwp = length(x_w_p);
x00 = x_w_p(1:nx)';
w = reshape(x_w_p(nx+1:nx*nt), nt-1, nx);
p = x_w_p(nx*nt+1 : nxwp)';
np = length(p);

log_2pi = log(2*pi);

% part 1: arriving cost
er = x00 - x0;
c_arrv = nx*log_2pi + sum(log(Cov_x0)) + sum( (er.^2) ./ Cov_x0 );
c_arrv = c_arrv/2;

   % arriving cost for parameters
   er_p = p0 - p;
   c_arrv_p = np*log_2pi + sum(log(Cov_p)) + sum( (er_p.^2) ./ Cov_p );
   c_arrv = c_arrv + c_arrv_p/2;
   

% part 2: cost regarding process noise
c_proc = 0;
for i = 1 : nt-1
    er = w(i,:);
    c_proc = c_proc + sum( (er.^2) ./ Cov_x );
end
c_proc = c_proc + nx*(nt-1)*log_2pi + (nt-1)*sum(log(Cov_x));
c_proc = c_proc/2;

   
% part 3: cost regarding measurement
c_obsv = 0;
x = x00;
for i = 1 : nt    
    ypred = feval(fn_obsv, x, p);
    er = ypred - yseq(i,:);
    c_obsv = c_obsv + sum( (er.^2) ./ Cov_y );
    
    if i == nt
        break;
    end
    
    % then move on to the next time step    
    ui = useq(i, :);
    x = feval(fn_state, x, p, ui, t_sample) + w(i,:);
end
c_obsv = c_obsv + ny*nt*log_2pi + nt*sum(log(Cov_y));
c_obsv = c_obsv/2;


f_mhe = c_arrv + c_proc + c_obsv;

if ~isreal(f_mhe)
    error('not real value of the objective function');
end

return;



%% Internal function to calculate objective function to be minimised
%   by adjusting only the states
function f_mhe = mhe_obj_states(x_w, nx, x0, Cov_x0, p, useq, yseq, fn_state, Cov_x, fn_obsv, Cov_y, t_sample)

[nt, ny] = size(yseq);
x00 = x_w(1:nx)';
w = reshape(x_w(nx+1:nx*nt), nt-1, nx);

log_2pi = log(2*pi);

% part 1: arriving cost
er = x00 - x0;
c_arrv = nx*log_2pi + sum(log(Cov_x0)) + sum( (er.^2) ./ Cov_x0 );
c_arrv = c_arrv/2;

% part 2: cost regarding process noise
c_proc = 0;
for i = 1 : nt-1    
    er = w(i,:);
    c_proc = c_proc + sum( (er.^2) ./ Cov_x );
end
c_proc = c_proc + nx*(nt-1)*log_2pi + (nt-1)*sum(log(Cov_x));
c_proc = c_proc/2;


% part 3: cost regarding measurement
c_obsv = 0;
x = x00;
for i = 1 : nt    
    ypred = feval(fn_obsv, x, p);
    er = ypred - yseq(i,:);
    c_obsv = c_obsv + sum( (er.^2) ./ Cov_y );
    
    if i == nt
        break;
    end
    
    % then move on to the next time step    
    ui = useq(i, :);
    x = feval(fn_state, x, p, ui, t_sample) + w(i,:);    
end
c_obsv = c_obsv + ny*nt*log_2pi + nt*sum(log(Cov_y));
c_obsv = c_obsv/2;


f_mhe = c_arrv + c_proc + c_obsv;

if ~isreal(f_mhe)
    error('not real value of the objective function');
end

return;


%% Internal function to calculate objective function to be minimised
%   by adjusting only the parameters
function f_mhe = mhe_obj_paras(p, nx, x0, useq, yseq, fn_state, fn_obsv, Cov_y, t_sample)
                     
[nt, ny] = size(yseq);

log_2pi = log(2*pi);

% cost regarding measurement
c_obsv = 0;
x = x0;
for i = 1 : nt    
    ypred = feval(fn_obsv, x, p);
    er = ypred - yseq(i,:);
    c_obsv = c_obsv + sum( (er.^2) ./ Cov_y );
    
    if i == nt
        break;
    end
    
    % then move on to the next time step    
    ui = useq(i, :);
    x = feval(fn_state, x, p, ui, t_sample);
end
c_obsv = c_obsv + ny*nt*log_2pi + nt*sum(log(Cov_y));
c_obsv = c_obsv/2;

f_mhe = c_obsv;

if ~isreal(f_mhe)
    error('not real value of the objective function');
end

return;