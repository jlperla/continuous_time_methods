% Modification of Ben Moll's: http://www.princeton.edu/~moll/HACTproject/option_simple_LCP.m
addpath('../../lib/'); %Adds library files used by the method.


%% Example 1 with the LCP based method. 
%Parameters
tic;
mu_bar = -0.01; %Drift.  Sign changes the upwind direction.
sigma_bar = 0.01; %Variance
S_bar = 10.0; %the value of stopping
gamma = 0.5; %u(x) = x^gamma

%Passing on to parameters.
p.rho = 0.05; %Discount rate
p.x_min = 0.01; %Reflecting barrier at x_min.  i.e. v'(x_min) = 0 as a boundary value
p.x_max = 1.0; %Reflecting barrier at x_max.  i.e. v'(x_max) = 0 as a boundary value

%Relevant functions for u(x), S(x), mu(x) and sigma(x) for a general diffusion dx_t = mu(x) dt + sigma(x) dW_t, for W_t brownian motion
p.u_x = @(x) x.^gamma; %u(x) = x^gamma in this example
p.S_x = @(x) S_bar.*ones(numel(x),1); %S(x) = S_bar in this example
p.mu_x = @(x) mu_bar * ones(numel(x),1); %i.e. mu(x) = mu_bar
p.sigma_2_x = @(x) (sigma_bar*x).^2; %i.e. sigma(x) = sigma_bar x
%TODO: In the old version, why does the sigma_bar multiple the grid x, but not the drift?  Is this a mistake in missing out the drift of the brownian motion term?

%Settings for the solution method
settings.I = 1000; %number of grid variables for x
settings.display = false;
settings.error_tolerance = 10^(-6);

%Create uniform grid and determine step sizes.
results = simple_optimal_stopping_diffusion(p, settings);
toc;
p2=plot(results.x,results.v,results.x,results.S,'--','LineWidth',2);
set(gca,'FontSize',16);
legend('v(x)','S(x)','Location','NorthWest');
xlabel('x');

