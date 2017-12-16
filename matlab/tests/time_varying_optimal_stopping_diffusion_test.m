%Using the unit testing framework in matlab.  See https://www.mathworks.com/help/matlab/matlab_prog/write-simple-test-case-with-functions.html
%To run tests:
% runtests %would run all of them in the current directory
% runtests('my_test') %runs just the my_test.m file
% runtests('my_test/my_function_test') %runs only `my_function_test function in `my_test'.

function tests = time_varying_optimal_stopping_diffusion_test
    tests = functiontests(localfunctions);
end

%This is run at the beginning of the test.  Not required.
function setupOnce(testCase)
    addpath('../lib/');
    
    testCase.TestData.tolerances.test_tol = 1e-9;    
    testCase.TestData.tolerances.test_tol_less = 1e-5;    
    testCase.TestData.tolerances.test_tol_much_less = 1e-3;    
    testCase.TestData.tolerances.default_csv_precision = '%.10f'; %Should be higher precision than test_tol    
end

function baseline_one_period_test(testCase)
    tolerances = testCase.TestData.tolerances;
    parameters.rho = 0.05; %Discount rate
    parameters.x_min = 0.01; %Reflecting barrier at x_min.  i.e. v'(x_min) = 0 as a boundary value
    parameters.x_max = 1.0; %Reflecting barrier at x_max.  i.e. v'(x_max) = 0 as a boundary value

    %Rewriting parameters entirely.
    mu_bar = -0.01; %Drift.  Sign changes the upwind direction.
    sigma_bar = 0.01; %Variance
    S_bar = 10.0; %the value of stopping
    gamma = 0.5; %u(x) = x^gamma
    parameters.rho = 0.05; %Discount ra te
    parameters.u = @(t,x) x.^gamma + 0*t; %u(x) = x^gamma in this example
    parameters.S = @(t,x) S_bar + 0*x + 0*t; %S(x) = S_bar in this example
    parameters.mu = @(t,x) mu_bar + 0*x + 0*t; %i.e. mu(x) = mu_bar
    parameters.sigma_2 = @(t,x) (sigma_bar*x).^2 + 0*t; %i.e. sigma(x) = sigma_bar x
    
    %Grids
    x_min = 0.01; %Reflecting barrier at x_min.  i.e. v'(x_min) = 0 as a boundary value
    x_max = 1.0; %Reflecting barrier at x_max.  i.e. v'(x_max) = 0 as a boundary value
    I = 100;

    parameters.t = 0; %One time period!
    parameters.x = linspace(x_min, x_max, I)';

    %Create uniform grid and determine step sizes.
    settings.method = 'yuval';
    tic;
    disp('yuval method');
    results = optimal_stopping_diffusion(parameters, settings);

    plot(parameters.x, results.v, parameters.x, results.S)
end

function baseline_repeated_period_test(testCase)
    tolerances = testCase.TestData.tolerances;
    parameters.rho = 0.05; %Discount rate
    parameters.x_min = 0.01; %Reflecting barrier at x_min.  i.e. v'(x_min) = 0 as a boundary value
    parameters.x_max = 1.0; %Reflecting barrier at x_max.  i.e. v'(x_max) = 0 as a boundary value

    %Rewriting parameters entirely.
    mu_bar = -0.01; %Drift.  Sign changes the upwind direction.
    sigma_bar = 0.01; %Variance
    S_bar = 10.0; %the value of stopping
    gamma = 0.5; %u(x) = x^gamma
    parameters.rho = 0.05; %Discount ra te
    parameters.u = @(t,x) x.^gamma + 0*t; %u(x) = x^gamma in this example
    parameters.S = @(t,x) S_bar + 0*x + 0*t; %S(x) = S_bar in this example
    parameters.mu = @(t,x) mu_bar + 0*x + 0*t; %i.e. mu(x) = mu_bar
    parameters.sigma_2 = @(t,x) (sigma_bar*x).^2 + 0*t; %i.e. sigma(x) = sigma_bar x
    
    %Grids
    x_min = 0.01; %Reflecting barrier at x_min.  i.e. v'(x_min) = 0 as a boundary value
    x_max = 1.0; %Reflecting barrier at x_max.  i.e. v'(x_max) = 0 as a boundary value
    I = 100;
    N = 3;
    parameters.t = linspace(0, 1, N)'; %One time period!
    parameters.x = linspace(x_min, x_max, I)';

    %Create uniform grid and determine step sizes.
    settings.method = 'yuval';
    tic;
    disp('yuval method');
    results = optimal_stopping_diffusion(parameters, settings);
    v = reshape(results.v, [I N]); %in three dimensions now
    S = reshape(results.S, [I N]); %in three dimensions now
    
%Useful to display the value funciton and the stopping value    
%    surf(parameters.t, parameters.x, v); hold on;
%    surf(parameters.t, parameters.x, S,'FaceAlpha',0.5, 'EdgeColor', 'none'); %SHows the S a litle different.
    
    %%TODO: Make sure that the v is nearly identical for of the time periods, and that it is the same as the optimal_stopping_diffusion we previously calculated.
end

function changing_S_test(testCase)
    tolerances = testCase.TestData.tolerances;
    parameters.rho = 0.05; %Discount rate
    parameters.x_min = 0.01; %Reflecting barrier at x_min.  i.e. v'(x_min) = 0 as a boundary value
    parameters.x_max = 1.0; %Reflecting barrier at x_max.  i.e. v'(x_max) = 0 as a boundary value

    %Rewriting parameters entirely.
    mu_bar = -0.01; %Drift.  Sign changes the upwind direction.
    sigma_bar = 0.01; %Variance
    S_bar = 10.0; %the value of stopping
    gamma = 0.5; %u(x) = x^gamma
    parameters.rho = 0.05; %Discount rate
    parameters.u = @(t,x) x.^gamma + 0*t;
    parameters.S = @(t,x) S_bar + 0*x + .1 * t; %NOTE: value increasing over time! Should have less stopping
    parameters.mu = @(t,x) mu_bar + 0*x + 0*t;
    parameters.sigma_2 = @(t,x) (sigma_bar*x).^2 + 0*t;
    
    %Grids
    x_min = 0.01; %Reflecting barrier at x_min.  i.e. v'(x_min) = 0 as a boundary value
    x_max = 1.0; %Reflecting barrier at x_max.  i.e. v'(x_max) = 0 as a boundary value
    I = 100;
    N = 10;
    parameters.t = linspace(0, 1, N)'; %One time period!
    parameters.x = linspace(x_min, x_max, I)';

    %Create uniform grid and determine step sizes.
    settings.method = 'yuval';
    tic;
    results = optimal_stopping_diffusion(parameters, settings);
    v = reshape(results.v, [I N]); %in three dimensions now
    S = reshape(results.S, [I N]); %in three dimensions now
%    surf(parameters.t, parameters.x, v); hold on;
%    surf(parameters.t, parameters.x, S,'FaceAlpha',0.5, 'EdgeColor', 'none'); %SHows the S a litle different.
    
    %%TODO: Make sure that this makes sense, that the v has less stopping than the one with an identical stopping value, etc.
end

