%Using the unit testing framework in matlab.  See https://www.mathworks.com/help/matlab/matlab_prog/write-simple-test-case-with-functions.html
%To run tests:
% runtests %would run all of them in the current directory
% runtests('my_test') %runs just the my_test.m file
% runtests('my_test/my_function_test') %runs only `my_function_test function in `my_test'.

function tests = simple_optimal_stopping_diffusion_test
    tests = functiontests(localfunctions);
end

%This is run at the beginning of the test.  Not required.
function setupOnce(testCase)
    addpath('../lib/');
    
    testCase.TestData.tolerances.test_tol = 1e-9;    
    testCase.TestData.tolerances.default_csv_precision = '%.10f'; %Should be higher precision than test_tol    
end
%To add in cleanup code, add here
%function teardownOnce(testCase)
%end

%This runs code prior to every test.  Not required
function setup(testCase)
    %Setup defaults.
    mu_bar = -0.01; %Drift.  Sign changes the upwind direction.
    sigma_bar = 0.01; %Variance
    S_bar = 10.0; %the value of stopping
    gamma = 0.5; %u(x) = x^gamma

    parameters.rho = 0.05; %Discount rate
    parameters.x_min = 0.01; %Reflecting barrier at x_min.  i.e. v'(x_min) = 0 as a boundary value
    parameters.x_max = 1.0; %Reflecting barrier at x_max.  i.e. v'(x_max) = 0 as a boundary value

    parameters.u_x = @(x) x.^gamma; %u(x) = x^gamma in this example
    parameters.S_x = @(x) S_bar.*ones(numel(x),1); %S(x) = S_bar in this example
    
    %Baseline test is GBM
    parameters.mu_x = @(x) mu_bar * x; %i.e. mu(x) = mu_bar * x
    parameters.sigma_2_x = @(x) (sigma_bar*x).^2; %i.e. sigma(x) = sigma_bar x
    settings.I = 1000; %number of grid variables for x
    settings.display = false; %Optional
    settings.error_tolerance = 10^(-6); %Optional    
    settings.method = 'Yuval'; %Optional, defaults to `Yuval'
    
    %These will be overwritten as required.
    testCase.TestData.baseline_parameters = parameters;
    testCase.TestData.baseline_settings = settings;
end

%This unpacks everything stored in the testCase
function [settings, parameters, tolerances] = unpack_setup(testCase)    
    settings = testCase.TestData.baseline_settings;
    parameters = testCase.TestData.baseline_parameters;
    tolerances = testCase.TestData.tolerances;
end

% Define an absolute tolerance for floating point comparisons

%A minimally modified version of the HACT code for comparison.  The main difference is generality and the boundary value at 0.  See /graveyard
function baseline_HACT_test(testCase)
    [settings, ~, tolerances] = unpack_setup(testCase);  
    
    %Rewriting parameters entirely.
    mu_bar = -0.01; %Drift.  Sign changes the upwind direction.
    sigma_bar = 0.01; %Variance
    S_bar = 10.0; %the value of stopping
    gamma = 0.5; %u(x) = x^gamma

    %Relevant functions for u(x), S(x), mu(x) and sigma(x) for a general diffusion dx_t = mu(x) dt + sigma(x) dW_t, for W_t brownian motion
    parameters.rho = 0.05; %Discount rate
    parameters.x_min = 0.1; %Reflecting barrier at x_min.  i.e. v'(x_min) = 0 as a boundary value
    parameters.x_max = 1.0; %Reflecting barrier at x_max.  i.e. v'(x_max) = 0 as a boundary value

    parameters.u_x = @(x) x.^gamma; %u(x) = x^gamma in this example
    parameters.S_x = @(x) S_bar.*ones(numel(x),1); %S(x) = S_bar in this example
    parameters.mu_x = @(x) mu_bar * ones(numel(x),1); %i.e. mu(x) = mu_bar
    parameters.sigma_2_x = @(x) (sigma_bar*x).^2; %i.e. sigma(x) = sigma_bar x
    
    %Create uniform grid and determine step sizes.
    results = simple_optimal_stopping_diffusion(parameters, settings);
    v = results.v;

    %Check all values
    v_old = dlmread(strcat(mfilename,'_1_v_output.csv')); %Loads old value, asserts identical.  Note that the precision of floating points in the .csv matters, and can't be lower than test_tol.
    verifyTrue(testCase, max(abs(v - v_old)) < tolerances.test_tol, 'Value of solution no longer matches HACT example');
end

function convex_u_x_test(testCase)
    [settings, ~, tolerances] = unpack_setup(testCase);  
    
    %Rewriting parameters entirely.
    mu_bar = -0.01; %Drift.  Sign changes the upwind direction.
    sigma_bar = 0.01; %Variance
    S_bar = 10.0; %the value of stopping
    gamma = 2; %u(x) = x^gamma

    %Relevant functions for u(x), S(x), mu(x) and sigma(x) for a general diffusion dx_t = mu(x) dt + sigma(x) dW_t, for W_t brownian motion
    parameters.rho = 0.05; %Discount rate
    parameters.x_min = 0.1; %Reflecting barrier at x_min.  i.e. v'(x_min) = 0 as a boundary value
    parameters.x_max = 1.0; %Reflecting barrier at x_max.  i.e. v'(x_max) = 0 as a boundary value

    parameters.u_x = @(x) x.^gamma; %u(x) = x^2 in this example
    parameters.S_x = @(x) S_bar.*ones(numel(x),1); %S(x) = S_bar in this example
    parameters.mu_x = @(x) mu_bar * ones(numel(x),1); %i.e. mu(x) = mu_bar
    parameters.sigma_2_x = @(x) (sigma_bar*x).^2; %i.e. sigma(x) = sigma_bar x
    
    %Create uniform grid and determine step sizes.
    results = simple_optimal_stopping_diffusion(parameters, settings);
    v = results.v;
    
    %dlmwrite(strcat(mfilename, '_2_v_output.csv'), results.v, 'precision', tolerances.default_csv_precision); %Uncomment to save again
    %Check all values
    v_old = dlmread(strcat(mfilename,'_2_v_output.csv')); %Loads old value, asserts identical.  Note that the precision of floating points in the .csv matters, and can't be lower than test_tol.
    verifyTrue(testCase, max(abs(v - v_old)) < tolerances.test_tol, 'Value of solution no longer matches negative u(x) for small x example');
end


function u_x_is_negative_for_small_x_test(testCase)
    [settings, ~, tolerances] = unpack_setup(testCase);  
    
    %Rewriting parameters entirely.
    mu_bar = -0.01; %Drift.  Sign changes the upwind direction.
    sigma_bar = 0.01; %Variance
    S_bar = 10.0; %the value of stopping

    %Relevant functions for u(x), S(x), mu(x) and sigma(x) for a general diffusion dx_t = mu(x) dt + sigma(x) dW_t, for W_t brownian motion
    parameters.rho = 0.05; %Discount rate
    parameters.x_min = 0.1; %Reflecting barrier at x_min.  i.e. v'(x_min) = 0 as a boundary value
    parameters.x_max = 1.0; %Reflecting barrier at x_max.  i.e. v'(x_max) = 0 as a boundary value

    parameters.u_x = @(x) x - 0.2; %u(x) = x - 0.2 in this example
    parameters.S_x = @(x) S_bar.*ones(numel(x),1); %S(x) = S_bar in this example
    parameters.mu_x = @(x) mu_bar * ones(numel(x),1); %i.e. mu(x) = mu_bar
    parameters.sigma_2_x = @(x) (sigma_bar*x).^2; %i.e. sigma(x) = sigma_bar x
    
    %Create uniform grid and determine step sizes.
    results = simple_optimal_stopping_diffusion(parameters, settings);
    v = results.v;
    
    %dlmwrite(strcat(mfilename, '_3_v_output.csv'), results.v, 'precision', tolerances.default_csv_precision); %Uncomment to save again
    %Check all values
    v_old = dlmread(strcat(mfilename,'_3_v_output.csv')); %Loads old value, asserts identical.  Note that the precision of floating points in the .csv matters, and can't be lower than test_tol.
    verifyTrue(testCase, max(abs(v - v_old)) < tolerances.test_tol, 'Value of solution no longer matches negative u(x) for small x example');
end


function negaive_S_x_test(testCase)
    [settings, ~, tolerances] = unpack_setup(testCase);  
    
    %Rewriting parameters entirely.
    mu_bar = -0.01; %Drift.  Sign changes the upwind direction.
    sigma_bar = 0.01; %Variance
    S_bar = -10.0; %the value of stopping
    gamma = 0.5; %u(x) = x^gamma

    %Relevant functions for u(x), S(x), mu(x) and sigma(x) for a general diffusion dx_t = mu(x) dt + sigma(x) dW_t, for W_t brownian motion
    parameters.rho = 0.05; %Discount rate
    parameters.x_min = 0.1; %Reflecting barrier at x_min.  i.e. v'(x_min) = 0 as a boundary value
    parameters.x_max = 1.0; %Reflecting barrier at x_max.  i.e. v'(x_max) = 0 as a boundary value

    parameters.u_x = @(x) x.^gamma; %u(x) = x^gamma in this example
    parameters.S_x = @(x) S_bar.*ones(numel(x),1); %S(x) = S_bar in this example
    parameters.mu_x = @(x) mu_bar * ones(numel(x),1); %i.e. mu(x) = mu_bar
    parameters.sigma_2_x = @(x) (sigma_bar*x).^2; %i.e. sigma(x) = sigma_bar x
    
    %Create uniform grid and determine step sizes.
    results = simple_optimal_stopping_diffusion(parameters, settings);
    v = results.v;
    
    %dlmwrite(strcat(mfilename, '_4_v_output.csv'), results.v, 'precision', tolerances.default_csv_precision); %Uncomment to save again
    %Check all values
    v_old = dlmread(strcat(mfilename,'_4_v_output.csv')); %Loads old value, asserts identical.  Note that the precision of floating points in the .csv matters, and can't be lower than test_tol.
    verifyTrue(testCase, max(abs(v - v_old)) < tolerances.test_tol, 'Value of solution no longer matches negative u(x) for small x example');
end



function S_x_increases_in_x_test(testCase)
    [settings, ~, tolerances] = unpack_setup(testCase);  
    
    %Rewriting parameters entirely.
    mu_bar = -0.01; %Drift.  Sign changes the upwind direction.
    sigma_bar = 0.01; %Variance
    S_bar = 10.0; %the value of stopping
    gamma = 0.5; %u(x) = x^gamma

    %Relevant functions for u(x), S(x), mu(x) and sigma(x) for a general diffusion dx_t = mu(x) dt + sigma(x) dW_t, for W_t brownian motion
    parameters.rho = 0.05; %Discount rate
    parameters.x_min = 0.1; %Reflecting barrier at x_min.  i.e. v'(x_min) = 0 as a boundary value
    parameters.x_max = 1.0; %Reflecting barrier at x_max.  i.e. v'(x_max) = 0 as a boundary value

    parameters.u_x = @(x) x.^gamma; %u(x) = x^gamma in this example
    parameters.S_x = @(x) x; %S(x) = x in this example
    parameters.mu_x = @(x) mu_bar * ones(numel(x),1); %i.e. mu(x) = mu_bar
    parameters.sigma_2_x = @(x) (sigma_bar*x).^2; %i.e. sigma(x) = sigma_bar x
    
    %Create uniform grid and determine step sizes.
    results = simple_optimal_stopping_diffusion(parameters, settings);
    v = results.v;
    
    %dlmwrite(strcat(mfilename, '_5_v_output.csv'), results.v, 'precision', tolerances.default_csv_precision); %Uncomment to save again
    %Check all values
    v_old = dlmread(strcat(mfilename,'_5_v_output.csv')); %Loads old value, asserts identical.  Note that the precision of floating points in the .csv matters, and can't be lower than test_tol.
    verifyTrue(testCase, max(abs(v - v_old)) < tolerances.test_tol, 'Value of solution no longer matches negative u(x) for small x example');
end

function S_x_decreases_in_x_test(testCase)
    [settings, ~, tolerances] = unpack_setup(testCase);  
    
    %Rewriting parameters entirely.
    mu_bar = -0.01; %Drift.  Sign changes the upwind direction.
    sigma_bar = 0.01; %Variance
    S_bar = 10.0; %the value of stopping
    gamma = 0.5; %u(x) = x^gamma

    %Relevant functions for u(x), S(x), mu(x) and sigma(x) for a general diffusion dx_t = mu(x) dt + sigma(x) dW_t, for W_t brownian motion
    parameters.rho = 0.05; %Discount rate
    parameters.x_min = 0.1; %Reflecting barrier at x_min.  i.e. v'(x_min) = 0 as a boundary value
    parameters.x_max = 1.0; %Reflecting barrier at x_max.  i.e. v'(x_max) = 0 as a boundary value

    parameters.u_x = @(x) x.^gamma; %u(x) = x^gamma in this example
    parameters.S_x = @(x) S_bar - x; %S(x) = S_bar - x in this example
    parameters.mu_x = @(x) mu_bar * ones(numel(x),1); %i.e. mu(x) = mu_bar
    parameters.sigma_2_x = @(x) (sigma_bar*x).^2; %i.e. sigma(x) = sigma_bar x
    
    %Create uniform grid and determine step sizes.
    results = simple_optimal_stopping_diffusion(parameters, settings);
    v = results.v;
    
    %dlmwrite(strcat(mfilename, '_6_v_output.csv'), results.v, 'precision', tolerances.default_csv_precision); %Uncomment to save again
    %Check all values
    v_old = dlmread(strcat(mfilename,'_6_v_output.csv')); %Loads old value, asserts identical.  Note that the precision of floating points in the .csv matters, and can't be lower than test_tol.
    verifyTrue(testCase, max(abs(v - v_old)) < tolerances.test_tol, 'Value of solution no longer matches negative u(x) for small x example');
end



function negative_mu_test(testCase)
    [settings, ~, tolerances] = unpack_setup(testCase);  
    
    %Rewriting parameters entirely.
    mu_bar = -0.01; %Drift.  Sign changes the upwind direction.
    sigma_bar = 0.01; %Variance
    S_bar = 10.0; %the value of stopping
    gamma = 0.5; %u(x) = x^gamma

    %Relevant functions for u(x), S(x), mu(x) and sigma(x) for a general diffusion dx_t = mu(x) dt + sigma(x) dW_t, for W_t brownian motion
    parameters.rho = 0.05; %Discount rate
    parameters.x_min = 0.1; %Reflecting barrier at x_min.  i.e. v'(x_min) = 0 as a boundary value
    parameters.x_max = 1.0; %Reflecting barrier at x_max.  i.e. v'(x_max) = 0 as a boundary value

    parameters.u_x = @(x) x.^gamma; %u(x) = x^0.5 in this example
    parameters.S_x = @(x) S_bar.*ones(numel(x),1); %S(x) = S_bar in this example
    parameters.mu_x = @(x) mu_bar * ones(numel(x),1); %i.e. mu(x) = mu_bar
    parameters.sigma_2_x = @(x) (sigma_bar*x).^2; %i.e. sigma(x) = sigma_bar x
    
    %Create uniform grid and determine step sizes.
    results = simple_optimal_stopping_diffusion(parameters, settings);
    v = results.v;
    
    %dlmwrite(strcat(mfilename, '_8_v_output.csv'), results.v, 'precision', tolerances.default_csv_precision); %Uncomment to save again
    %Check all values
    v_old = dlmread(strcat(mfilename,'_8_v_output.csv')); %Loads old value, asserts identical.  Note that the precision of floating points in the .csv matters, and can't be lower than test_tol.
    verifyTrue(testCase, max(abs(v - v_old)) < tolerances.test_tol, 'Value of solution no longer matches negative u(x) for small x example');
end

% function positive_mu_test(testCase)
%     [settings, ~, tolerances] = unpack_setup(testCase);  
%     
%     %Rewriting parameters entirely.
%     mu_bar = 0.01; %Drift.  Sign changes the upwind direction.
%     sigma_bar = 0.01; %Variance
%     S_bar = 10.0; %the value of stopping
%     gamma = 0.5; %u(x) = x^gamma
% 
%     %Relevant functions for u(x), S(x), mu(x) and sigma(x) for a general diffusion dx_t = mu(x) dt + sigma(x) dW_t, for W_t brownian motion
%     parameters.rho = 0.05; %Discount rate
%     parameters.x_min = 0.1; %Reflecting barrier at x_min.  i.e. v'(x_min) = 0 as a boundary value
%     parameters.x_max = 1.0; %Reflecting barrier at x_max.  i.e. v'(x_max) = 0 as a boundary value
% 
%     parameters.u_x = @(x) x.^gamma; %u(x) = x^0.5 in this example
%     parameters.S_x = @(x) S_bar.*ones(numel(x),1); %S(x) = S_bar in this example
%     parameters.mu_x = @(x) mu_bar * ones(numel(x),1); %i.e. mu(x) = mu_bar
%     parameters.sigma_2_x = @(x) (sigma_bar*x).^2; %i.e. sigma(x) = sigma_bar x
%     
%     %Create uniform grid and determine step sizes.
%     results = simple_optimal_stopping_diffusion(parameters, settings);
%     v = results.v;
%     
%     dlmwrite(strcat(mfilename, '_9_v_output.csv'), results.v, 'precision', tolerances.default_csv_precision); %Uncomment to save again
%     %Check all values
%     v_old = dlmread(strcat(mfilename,'_9_v_output.csv')); %Loads old value, asserts identical.  Note that the precision of floating points in the .csv matters, and can't be lower than test_tol.
%     verifyTrue(testCase, max(abs(v - v_old)) < tolerances.test_tol, 'Value of solution no longer matches negative u(x) for small x example');
% end

function zero_mu_test(testCase)
    [settings, ~, tolerances] = unpack_setup(testCase);  
    
    %Rewriting parameters entirely.
    mu_bar = 0; %Drift.  Sign changes the upwind direction.
    sigma_bar = 0.01; %Variance
    S_bar = 10.0; %the value of stopping
    gamma = 0.5; %u(x) = x^gamma

    %Relevant functions for u(x), S(x), mu(x) and sigma(x) for a general diffusion dx_t = mu(x) dt + sigma(x) dW_t, for W_t brownian motion
    parameters.rho = 0.05; %Discount rate
    parameters.x_min = 0.1; %Reflecting barrier at x_min.  i.e. v'(x_min) = 0 as a boundary value
    parameters.x_max = 1.0; %Reflecting barrier at x_max.  i.e. v'(x_max) = 0 as a boundary value

    parameters.u_x = @(x) x.^gamma; %u(x) = x^0.5 in this example
    parameters.S_x = @(x) S_bar.*ones(numel(x),1); %S(x) = S_bar in this example
    parameters.mu_x = @(x) mu_bar * ones(numel(x),1); %i.e. mu(x) = mu_bar
    parameters.sigma_2_x = @(x) (sigma_bar*x).^2; %i.e. sigma(x) = sigma_bar x
    
    %Create uniform grid and determine step sizes.
    results = simple_optimal_stopping_diffusion(parameters, settings);
    v = results.v;
    
    %dlmwrite(strcat(mfilename, '_10_v_output.csv'), results.v, 'precision', tolerances.default_csv_precision); %Uncomment to save again
    %Check all values
    v_old = dlmread(strcat(mfilename,'_10_v_output.csv')); %Loads old value, asserts identical.  Note that the precision of floating points in the .csv matters, and can't be lower than test_tol.
    verifyTrue(testCase, max(abs(v - v_old)) < tolerances.test_tol, 'Value of solution no longer matches negative u(x) for small x example');
end

% function negative_mu_min_and_positive_mu_max_test(testCase)
%     [settings, ~, tolerances] = unpack_setup(testCase);  
%     
%     %Rewriting parameters entirely.
%     mu_bar = -0.01; %Drift.  Sign changes the upwind direction.
%     sigma_bar = 0.01; %Variance
%     S_bar = 10.0; %the value of stopping
%     gamma = 0.5; %u(x) = x^gamma
% 
%     %Relevant functions for u(x), S(x), mu(x) and sigma(x) for a general diffusion dx_t = mu(x) dt + sigma(x) dW_t, for W_t brownian motion
%     parameters.rho = 0.05; %Discount rate
%     parameters.x_min = 0.1; %Reflecting barrier at x_min.  i.e. v'(x_min) = 0 as a boundary value
%     parameters.x_max = 1.0; %Reflecting barrier at x_max.  i.e. v'(x_max) = 0 as a boundary value
% 
%     parameters.u_x = @(x) x.^gamma; %u(x) = x^0.5 in this example
%     parameters.S_x = @(x) S_bar.*ones(numel(x),1); %S(x) = S_bar in this example
%     parameters.mu_x = @(x) x - 0.5; %i.e. mu(x) = x - 0.5;
%     parameters.sigma_2_x = @(x) (sigma_bar*x).^2; %i.e. sigma(x) = sigma_bar x
%     
%     %Create uniform grid and determine step sizes.
%     results = simple_optimal_stopping_diffusion(parameters, settings);
%     v = results.v;
%     
%     dlmwrite(strcat(mfilename, '_11_v_output.csv'), results.v, 'precision', tolerances.default_csv_precision); %Uncomment to save again
%     %Check all values
%     v_old = dlmread(strcat(mfilename,'_11_v_output.csv')); %Loads old value, asserts identical.  Note that the precision of floating points in the .csv matters, and can't be lower than test_tol.
%     verifyTrue(testCase, max(abs(v - v_old)) < tolerances.test_tol, 'Value of solution no longer matches negative u(x) for small x example');
% end
% 
% 
% 
% function positive_mu_min_and_negative_mu_max_test(testCase)
%     [settings, ~, tolerances] = unpack_setup(testCase);  
%     
%     %Rewriting parameters entirely.
%     mu_bar = -0.01; %Drift.  Sign changes the upwind direction.
%     sigma_bar = 0.01; %Variance
%     S_bar = 10.0; %the value of stopping
%     gamma = 0.5; %u(x) = x^gamma
% 
%     %Relevant functions for u(x), S(x), mu(x) and sigma(x) for a general diffusion dx_t = mu(x) dt + sigma(x) dW_t, for W_t brownian motion
%     parameters.rho = 0.05; %Discount rate
%     parameters.x_min = 0.1; %Reflecting barrier at x_min.  i.e. v'(x_min) = 0 as a boundary value
%     parameters.x_max = 1.0; %Reflecting barrier at x_max.  i.e. v'(x_max) = 0 as a boundary value
% 
%     parameters.u_x = @(x) x.^gamma; %u(x) = x^0.5 in this example
%     parameters.S_x = @(x) S_bar.*ones(numel(x),1); %S(x) = S_bar in this example
%     parameters.mu_x = @(x) -x + 0.5; %i.e. mu(x) = -x + 0.5;
%     parameters.sigma_2_x = @(x) (sigma_bar*x).^2; %i.e. sigma(x) = sigma_bar x
%     
%     %Create uniform grid and determine step sizes.
%     results = simple_optimal_stopping_diffusion(parameters, settings);
%     v = results.v;
%     
%     dlmwrite(strcat(mfilename, '_12_v_output.csv'), results.v, 'precision', tolerances.default_csv_precision); %Uncomment to save again
%     %Check all values
%     v_old = dlmread(strcat(mfilename,'_12_v_output.csv')); %Loads old value, asserts identical.  Note that the precision of floating points in the .csv matters, and can't be lower than test_tol.
%     verifyTrue(testCase, max(abs(v - v_old)) < tolerances.test_tol, 'Value of solution no longer matches negative u(x) for small x example');
% end


% function negative_mu_and_zero_sigma_test(testCase)
%     [settings, ~, tolerances] = unpack_setup(testCase);  
%     
%     %Rewriting parameters entirely.
%     mu_bar = -0.01; %Drift.  Sign changes the upwind direction.
%     sigma_bar = 0; %Variance
%     S_bar = 10.0; %the value of stopping
%     gamma = 0.5; %u(x) = x^gamma
% 
%     %Relevant functions for u(x), S(x), mu(x) and sigma(x) for a general diffusion dx_t = mu(x) dt + sigma(x) dW_t, for W_t brownian motion
%     parameters.rho = 0.05; %Discount rate
%     parameters.x_min = 0.1; %Reflecting barrier at x_min.  i.e. v'(x_min) = 0 as a boundary value
%     parameters.x_max = 1.0; %Reflecting barrier at x_max.  i.e. v'(x_max) = 0 as a boundary value
% 
%     parameters.u_x = @(x) x.^gamma; %u(x) = x^0.5 in this example
%     parameters.S_x = @(x) S_bar.*ones(numel(x),1); %S(x) = S_bar in this example
%     parameters.mu_x = @(x) mu_bar * ones(numel(x),1); %i.e. mu(x) = mu_bar
%     parameters.sigma_2_x = @(x) sigma_bar * ones(numel(x),1); %i.e. sigma(x) = sigma_bar
%     
%     %Create uniform grid and determine step sizes.
%     results = simple_optimal_stopping_diffusion(parameters, settings);
%     v = results.v;
%     
%     dlmwrite(strcat(mfilename, '_13_v_output.csv'), results.v, 'precision', tolerances.default_csv_precision); %Uncomment to save again
%     %Check all values
%     v_old = dlmread(strcat(mfilename,'_13_v_output.csv')); %Loads old value, asserts identical.  Note that the precision of floating points in the .csv matters, and can't be lower than test_tol.
%     verifyTrue(testCase, max(abs(v - v_old)) < tolerances.test_tol, 'Value of solution no longer matches negative u(x) for small x example');
% end


% function positive_mu_and_zero_sigma_test(testCase)
%     [settings, ~, tolerances] = unpack_setup(testCase);  
%     
%     %Rewriting parameters entirely.
%     mu_bar = 0.01; %Drift.  Sign changes the upwind direction.
%     sigma_bar = 0; %Variance
%     S_bar = 10.0; %the value of stopping
%     gamma = 0.5; %u(x) = x^gamma
% 
%     %Relevant functions for u(x), S(x), mu(x) and sigma(x) for a general diffusion dx_t = mu(x) dt + sigma(x) dW_t, for W_t brownian motion
%     parameters.rho = 0.05; %Discount rate
%     parameters.x_min = 0.1; %Reflecting barrier at x_min.  i.e. v'(x_min) = 0 as a boundary value
%     parameters.x_max = 1.0; %Reflecting barrier at x_max.  i.e. v'(x_max) = 0 as a boundary value
% 
%     parameters.u_x = @(x) x.^gamma; %u(x) = x^0.5 in this example
%     parameters.S_x = @(x) S_bar.*ones(numel(x),1); %S(x) = S_bar in this example
%     parameters.mu_x = @(x) mu_bar * ones(numel(x),1); %i.e. mu(x) = mu_bar
%     parameters.sigma_2_x = @(x) sigma_bar * ones(numel(x),1); %i.e. sigma(x) = sigma_bar
%     
%     %Create uniform grid and determine step sizes.
%     results = simple_optimal_stopping_diffusion(parameters, settings);
%     v = results.v;
%     
%     dlmwrite(strcat(mfilename, '_14_v_output.csv'), results.v, 'precision', tolerances.default_csv_precision); %Uncomment to save again
%     %Check all values
%     v_old = dlmread(strcat(mfilename,'_14_v_output.csv')); %Loads old value, asserts identical.  Note that the precision of floating points in the .csv matters, and can't be lower than test_tol.
%     verifyTrue(testCase, max(abs(v - v_old)) < tolerances.test_tol, 'Value of solution no longer matches negative u(x) for small x example');
% end

%This test runs the test case with only the default parameters in settings.
function default_parameters_test(testCase)
    [~, parameters, tolerances] = unpack_setup(testCase); 
    
    %default parameters, but note that settings is not used.    
    settings.I = 1000; %Only the number of points is provided.
    
    %Create uniform grid and determine step sizes.
    results = simple_optimal_stopping_diffusion(parameters, settings);    
    dlmwrite(strcat(mfilename,'_22_v_output.csv'), results.v, 'precision', tolerances.default_csv_precision); %To save results again
    v_old = dlmread(strcat(mfilename,'_22_v_output.csv')); %Loads old value, asserts identical.  Note that the precision of floating points in the .csv matters, and can't be lower than test_tol.    
    verifyTrue(testCase, max(abs(results.v - v_old)) < tolerances.test_tol, 'Value of solution no longer matches default value');        
end
