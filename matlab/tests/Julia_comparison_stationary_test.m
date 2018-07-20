% This is test function that try to replicate Julia problem with same set
% up

%Using the unit testing framework in matlab.  See https://www.mathworks.com/help/matlab/matlab_prog/write-simple-test-case-with-functions.html
%To run tests:
% runtests %would run all of them in the current directory
% runtests('my_test') %runs just the my_test.m file
% runtests('my_test/my_function_test') %runs only `my_function_test function in `my_test'.

function tests = Julia_comparison_stationary_test
    tests = functiontests(localfunctions);
end

%This is run at the beginning of the test.  Not required.
function setupOnce(testCase)
    addpath('../lib/');
    testCase.TestData.tolerances.test_tol = 1e-9;    
    testCase.TestData.tolerances.lower_test_tol = 1e-8;    %For huge matrices, the inf norm can get a little different.
    testCase.TestData.tolerances.default_csv_precision = '%.10f'; %Should be higher precision than tolerances.test_tol
end

function uniform_more_grid_test(testCase)
    tolerances = testCase.TestData.tolerances;
    r = 0.05;
    zeta = 14.5;
    gamma = 0.005;
    g = 0.020758;
    mu_x = @(x) gamma-g;
    sigma_bar = 0.02;
    sigma_2_x = @(x) (sigma_bar).^2+0.0*x;
    u_x = @(x) exp(x);
    rho = r-g;
    x_min = 0.0;
    x_max = 5.0;
    I = 500;
    x = linspace(x_min,x_max,I)';
    A = discretize_univariate_diffusion(x, mu_x(x), sigma_2_x(x));
    u = u_x(x);
    
    %Solve the simple problem: rho v(x) = u(x) + A v(x) for the above process.
    [v, success] = simple_HJBE_discretized_univariate(A, x, u, rho);
    %dlmwrite(strcat(mfilename, '_1_v_output.csv'), v, 'precision', tolerances.default_csv_precision); %Uncomment to save again
    %v_check = dlmread(strcat(mfilename, '_1_v_output.csv'));    
    %verifyTrue(testCase,norm(v - v_check, Inf) < tolerances.test_tol, 'v value no longer matches');
    %verifyTrue(testCase, success==true, 'unsuccesful');
    
    %Solve with a uniform grid and check if similar after interpolation
    x_2 = linspace(x_min, x_max, 701)'; %Twice as many points to be sure.
    A_2 = discretize_univariate_diffusion(x_2, mu_x(x_2), sigma_2_x(x_2));
    %Solve the simple problem: rho v(x) = u(x) + A v(x) for the above process.
    [v_2, success] = simple_HJBE_discretized_univariate(A_2, x_2, u_x(x_2), rho);

    %Make sure within range.  This seems too large.
    verifyTrue(testCase, norm(interp1(x_2, v_2, x) - v,Inf) < 0.02, 'Not within range of interpolation');    
end