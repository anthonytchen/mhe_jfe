% The test script to demonstrate the use of moving horizon estimation
%   code developed at the University of Surrey
%
% Author: Tao Chen (t.chen@surrey.ac.uk)
% Revision: 13 March 2018
%
% returns:
%
% Refs:
%

%function Demo_fedBio()


% +++ set initial states +++
x0 = [1 10.5/7 7]; % [Cs Cx V]
nx = length(x0);
y0 = fedBio_mea(x0);

% +++ initliazation for simulation +++
t_sample = 2; % sampling interval, hours
t_simu = 40; % simulation time, hours
n_sample = t_simu / t_sample;
n_mhe = 10; % time horizon for moving horizon estimation

std_x = [0.001 0.001 0.001]; % [Cs Cx V], standard deviation for state noise
std_y = [0.1 0.05]; % [Cs V], standard deviation for measurement noise
nz = length(std_y);

% model parameters   
Ks = 50; Ki = 0.5;
p = log([ Ks, Ki ]); % using log transform to ensure estimated parameters are positive


u_nom = ones(n_sample,1) * 0.075; % constant feed rate, 0.075 L/hr
np = length(p);


% +++ initialisation for MHE +++

% MHE for state & parameter estimation
x_mhe_sp = zeros(n_sample+1, nx);
x_mhe_sp(1,:) = x0;
x_mhe_sp_smooth = x_mhe_sp; % smoothed version of states in the past

% MHE for state estimation only
x_mhe_s = zeros(n_sample+1, nx);
x_mhe_s(1,:) = x0;
x_mhe_s_smooth = x_mhe_s; % smoothed version of states in the past


Cov_x = std_x.^2;
Cov_y = std_y.^2;
Cov_p = ones(size(p)) * 1e10; % a very large covariance for parameter, meaning no prior for parameters

paras_mhe = zeros(n_sample+1, np); 
paras_mhe(1,:) = [log(2) log(2)] + p; % Initial guess of the parameters are twice of the true values


%% simulation, and MHE state/parameter estimation

for i = 1 : n_sample
        
    fprintf('Iter %d\n', i);
          
    %% 1. simulation to generate 'measured' data
    
    if i == 1
		xc = x0;
	else
		xc = simu.x(i-1,:); 
    end
    u = u_nom(i);
    
    xc = fedBio_state(xc, p, u, t_sample);
    simu.xtrue(i,:) = xc;
    xc = xc + randn(1, nx) .* std_x; xc(xc<0)=0;
	simu.x(i,:) = xc;
    
    z = fedBio_mea(xc); % 'z' and 'y' used interchangeably
    simu.ztrue(i,:) = z; 
    z = z + randn(1, nz) .* std_y;   z(z<0)=0;
    simu.z(i,:) = z;
    	    
    
    %% 2. MHE for state and parameter estimation    
    
    % preparing data for MHE
    xt = x_mhe_sp(i,:); % MHE estimate of the current states
    pt = paras_mhe(i,:); % MHE estimate of the parameters
    x = fedBio_state(xt, pt, u, t_sample);
    if i < n_mhe+0.5 % m_mhe is the window size of the moving horizon
        xseq = [x_mhe_sp_smooth(1:i,:); x]; % this is the initial value to be optimised by MHE
        yseq = [y0; simu.z(1:i,:)];
        useq = u_nom(1:i,:);
    else
        xseq = [x_mhe_sp_smooth(i-n_mhe+1:i,:); x]; % this is the initial value to be optimised by MHE
        yseq = simu.z(i-n_mhe:i,:);
        useq = u_nom( i-n_mhe+1:i, : );
    end           
    
    % call MHE function
    [x, pa, f_obj] = mhe_w(xseq, useq, yseq, pt, 'fedBio_state', 'fedBio_mea', Cov_x, Cov_y, Cov_p, t_sample, 1);    
    paras_mhe(i+1,:) = pa';
    
    if i < n_mhe
        x_mhe_sp_smooth(1:i+1,:) = x;
        x_mhe_sp(i+1,:) = x(i+1,:);
    else
        x_mhe_sp_smooth(i-n_mhe+1:i+1,:) = x;
        x_mhe_sp(i+1,:) = x(size(x,1),:);
    end
    
    
    
    %% 3. MHE for state estimation only
    
    % preparing data for MHE
    xt = x_mhe_s(i,:); % MHE estimate of the current states
    pt = paras_mhe(1,:); % taking the initial guess for the parameters
    x = fedBio_state(xt, pt, u, t_sample);
    if i < n_mhe+0.5 % m_mhe is the window size of the moving horizon
        xseq = [x_mhe_s_smooth(1:i,:); x]; % this is the initial value to be optimised by MHE
        yseq = [y0; simu.z(1:i,:)];
        useq = u_nom(1:i,:);
    else
        xseq = [x_mhe_s_smooth(i-n_mhe+1:i,:); x]; % this is the initial value to be optimised by MHE
        yseq = simu.z(i-n_mhe:i,:);
        useq = u_nom( i-n_mhe+1:i, : );
    end           
    
    % call MHE function
    [x, pa, f_obj] = mhe_w(xseq, useq, yseq, pt, 'fedBio_state', 'fedBio_mea', Cov_x, Cov_y, Cov_p, t_sample, 0);    
    
    if i < n_mhe
        x_mhe_s_smooth(1:i+1,:) = x;
        x_mhe_s(i+1,:) = x(i+1,:);
    else
        x_mhe_s_smooth(i-n_mhe+1:i+1,:) = x;
        x_mhe_s(i+1,:) = x(size(x,1),:);
    end
end

%% Plot the results

t_seq = 0:t_sample:t_simu;
plot(t_seq, [x0(2); simu.xtrue(:,2)], 'r-', t_seq, x_mhe_s(:,2), 'k-o', t_seq, x_mhe_sp(:,2), 'b-x');
legend('Simu', 'State Est.', 'State & Para. Est.');
xlabel('Time (hr)');
ylabel('X [g/L]');





