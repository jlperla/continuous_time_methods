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
    u_x = @(x) exp(0.5*x);
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
    
    figure()

    v_2int=interp1(x_2, v_2, x);
    dif=v-v_2int;
    plot(x,dif,'LineWidth',2);hold on
    title('difference of v and v_2 along z_grid, z_bar=5')

    %Make sure within range.  This seems too large.
    verifyTrue(testCase, norm(interp1(x_2, v_2, x) - v,Inf) < 0.02, 'Not within range of interpolation');    
end

function uniform_more_grid_test2(testCase)
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
    x_max = 2.0; % lower xbar
    I = 301;
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

    figure()

    v_2int=interp1(x_2, v_2, x);
    dif=v-v_2int;
    plot(x,dif,'LineWidth',2);hold on
    title('difference of v and v_2 along z_grid,z_bar=2')
    %Make sure within range.  This seems too large.
    verifyTrue(testCase, norm(interp1(x_2, v_2, x) - v,Inf) < 0.02, 'Not within range of interpolation');    
end

function uniform_plot_test(testCase)
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
    x_max = 5.0; % lower xbar
    
    %Solve the simple problem: rho v(x) = u(x) + A v(x) for the above process.
    I_c=[101 201 301 401 501];
    for i=1:5
    I = I_c(i);
    x = linspace(x_min,x_max,I)';
    A = discretize_univariate_diffusion(x, mu_x(x), sigma_2_x(x));
    u = u_x(x);
    [v_{i}, success] = simple_HJBE_discretized_univariate(A, x, u, rho);
    x_{i}=x;
    end
    
    figure()
    for i=1:5
    plot(x_{i}(end-5:end),v_{i}(end-5:end)); hold on
    end
    legend('101','201','301','401','500')
    title('value function at last 5 grid point')
  
end

function uniform_to_nonuniform_test(testCase)
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
    x_max = 2.0; % lower xbar
    
    %Solve the simple problem: rho v(x) = u(x) + A v(x) for the above process.
    I = 301;
    x = linspace(x_min,x_max,I)';
    A = discretize_univariate_diffusion(x, mu_x(x), sigma_2_x(x));
    u = u_x(x);
    [v, success] = simple_HJBE_discretized_univariate(A, x, u, rho);
    
    I2=floor(I*1/3); % propotional adding grid points
    x_2 = unique([linspace(x_min, x_max, I)'; linspace(1.7,x_max,I2)']); %Twice as many points to be sure.
    A_2 = discretize_univariate_diffusion(x_2, mu_x(x_2), sigma_2_x(x_2));
    %Solve the simple problem: rho v(x) = u(x) + A v(x) for the above process.
    [v_2, success] = simple_HJBE_discretized_univariate(A_2, x_2, u_x(x_2), rho);
    v_2int=interp1(x_2, v_2, x);
    
    %Make sure within range.  This seems too large.
    verifyTrue(testCase, norm(interp1(x_2, v_2, x) - v,Inf) < 0.02, 'Not within range of interpolation');    
end