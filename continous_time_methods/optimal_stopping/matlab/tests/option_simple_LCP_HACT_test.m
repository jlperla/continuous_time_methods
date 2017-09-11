% Generates the output data
% call original code
option_simple_LCP();

% Define an absolute tolerance
test_tol = 1e-10;

%% Test 1: baseline HACT comparison
%load results from the baseline setup.
addpath('../');
load output_option_simple_LCP_HACT_raw.mat 
test_output_simple_HACT = test_output;  

%Run the new code.
mu_bar = -0.01; %Drift.  Sign changes the upwind direction.
sigma_bar = 0.01; %Variance
S_bar = 10.0; %the value of stopping
x_min = 0.0; %TODO: Is there something crucial here to check? v(x_min) < S_bar for example for a boundary value?  Or is there even a binding value?
x_max = 1.0; %Reflecting barrier at x_max.  i.e. v'(x_max) = 0 as a boundary value
gamma = 0.5; %u(x) = x^gamma

%Passing on to parameters.
p.rho = 0.05; %Discount rate
p.x_min = 0.0; %TODO: Is there something crucial here to check? v(x_min) < S_bar for example for a boundary value?  Or is there even a binding value?
p.x_max = 1.0; %Reflecting barrier at x_max.  i.e. v'(x_max) = 0 as a boundary value

%Relevant functions for u(x), S(x), mu(x) and sigma(x) for a general diffusion dx_t = mu(x) dt + sigma(x) dW_t, for W_t brownian motion
p.u_x = @(x) x.^gamma; %u(x) = x^gamma in this example
p.S_x = @(x) S_bar.*ones(numel(x),1); %S(x) = S_bar in this example
p.mu_x = @(x) mu_bar * ones(numel(x),1); %i.e. mu(x) = mu_bar
p.sigma_2_x = @(x) (sigma_bar*x).^2; %i.e. sigma(x) = sigma_bar x
%TODO: In the old version, why does the sigma_bar multiple the grid x, but not the drift?  Is this a mistake in missing out the drift of the brownian motion term?

%Settings for the solution method
settings.N_x = 1000; %number of grid variables for x
settings.display = false;
settings.error_tolerance = 10^(-6);

%Create uniform grid and determine step sizes.
results = simple_optimal_stopping_diffusion(p, settings);

%Check all values
assert(max(abs(results.x - test_output_simple_HACT.x)) < test_tol, 'grid different');
assert(max(abs(results.v - test_output_simple_HACT.v)) < test_tol, 'value solution different');
assert(max(abs(results.S - test_output_simple_HACT.S)) < test_tol, 'S solution different');

