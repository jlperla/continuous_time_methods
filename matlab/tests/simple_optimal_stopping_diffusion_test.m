addpath('../lib/');
main_script_tested = 'simple_optimal_stopping_diffusion';

% Define an absolute tolerance for floating point comparisons
test_tol = 1e-9;

%% Test 1: baseline HACT comparison
%option_simple_LCP(); %A minimally modified version of the HACT code for comparison.  The main difference is generality and the boundary value at 0.  See /graveyard
%In this test, the parameters need to be maintained to match the HACT code
mu_bar = -0.01; %Drift.  Sign changes the upwind direction.
sigma_bar = 0.01; %Variance
S_bar = 10.0; %the value of stopping
gamma = 0.5; %u(x) = x^gamma

%Relevant functions for u(x), S(x), mu(x) and sigma(x) for a general diffusion dx_t = mu(x) dt + sigma(x) dW_t, for W_t brownian motion
p.rho = 0.05; %Discount rate
p.x_min = 0.1; %Reflecting barrier at x_min.  i.e. v'(x_min) = 0 as a boundary value
p.x_max = 1.0; %Reflecting barrier at x_max.  i.e. v'(x_max) = 0 as a boundary value

p.u_x = @(x) x.^gamma; %u(x) = x^gamma in this example
p.S_x = @(x) S_bar.*ones(numel(x),1); %S(x) = S_bar in this example
p.mu_x = @(x) mu_bar * ones(numel(x),1); %i.e. mu(x) = mu_bar
p.sigma_2_x = @(x) (sigma_bar*x).^2; %i.e. sigma(x) = sigma_bar x

%Settings for the solution method
settings.I = 1000; %number of grid variables for x
settings.display = false; %Optional
settings.error_tolerance = 10^(-6); %Optional
settings.method = 'Yuval'; %Optional, defaults to `Yuval'

%Create uniform grid and determine step sizes.
results = simple_optimal_stopping_diffusion(p, settings);
v = results.v;

%Check all values
v_old = dlmread(strcat(main_script_tested,'_1_v_output.csv')); %Loads old value, asserts identical.  Note that the precision of floating points in the .csv matters, and can't be lower than test_tol.
assert(max(abs(v - v_old)) < test_tol, 'Value of solution no longer matches HACT example');

%% Test 2: Check default values for settings
mu_bar = -0.01; %Drift.  Sign changes the upwind direction.
sigma_bar = 0.01; %Variance
S_bar = 10.0; %the value of stopping
x_min = 0.01; %TODO: Is there something crucial here to check? v(x_min) < S_bar for example for a boundary value?  Or is there even a binding value?
x_max = 1.0; %Reflecting barrier at x_max.  i.e. v'(x_max) = 0 as a boundary value
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

%Settings for the solution method
settings.I = 1000; %number of grid variables for x

%Create uniform grid and determine step sizes.
results = simple_optimal_stopping_diffusion(p, settings);
v = results.v;
%dlmwrite(strcat(main_script_tested,'_2_v_output.csv'), v, 'precision', default_csv_precision); %To save results again

v_old = dlmread(strcat(main_script_tested,'_2_v_output.csv')); %Loads old value, asserts identical.  Note that the precision of floating points in the .csv matters, and can't be lower than test_tol.
assert(max(abs(v - v_old)) < test_tol, 'Value of solution no longer matches default value');

