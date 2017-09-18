%Using the unit testing framework in matlab.  See https://www.mathworks.com/help/matlab/matlab_prog/write-simple-test-case-with-functions.html
%To run tests:
% runtests %would run all of them in the current directory
% runtests('my_test') %runs just the my_test.m file
% runtests('my_test/my_function_test') %runs only `my_function_test function in `my_test'.

function tests = HJBE_discretized_univariate_test
    tests = functiontests(localfunctions);
end

%This is run at the beginning of the test.  Not required.
function setupOnce(testCase)
    addpath('../lib/');
    testCase.TestData.tolerances.test_tol = 1e-9;    
    testCase.TestData.tolerances.default_csv_precision = '%.10f'; %Should be higher precision than tolerances.test_tol
end

function simple_value_function_test(testCase)
    tolerances = testCase.TestData.tolerances;
    mu_x = @(x) -0.01 * x;
    sigma_bar = 0.1;
    sigma_2_x = @(x) (sigma_bar*x).^2;
    u_x = @(x) exp(x);
    rho = 0.05;
    x_min = 1;
    x_max = 2;
    I = 20;
    x = linspace(x_min, x_max, I)';
    A = discretize_univariate_diffusion(x, mu_x(x), sigma_2_x(x));     
    u = u_x(x);
    
    %Solve the simple problem: rho v(x) = u(x) + A v(x) for the above process.
    [v, success] = simple_HJBE_discretized_univariate(A, x, u, rho);
    %dlmwrite(strcat(mfilename, '_1_v_output.csv'), v, 'precision', tolerances.default_csv_precision); %Uncomment to save again
    v_check = dlmread(strcat(mfilename, '_1_v_output.csv'));    
    verifyTrue(testCase,norm(v - v_check, Inf) < tolerances.test_tol, 'v value no longer matches');
    verifyTrue(testCase, success==true, 'unsuccesful');
end

function bigger_value_function_test(testCase)
    tolerances = testCase.TestData.tolerances;
    mu_x = @(x) -0.01 * x;
    sigma_bar = 0.1;
    sigma_2_x = @(x) (sigma_bar*x).^2;
    u_x = @(x) log(x);
    rho = 0.05;
    x_min = .01;
    x_max = 10;
    I = 10000;
    x = linspace(x_min, x_max, I)';
    A = discretize_univariate_diffusion(x, mu_x(x), sigma_2_x(x));         
    u = u_x(x);
    
    %Solve the simple problem: rho v(x) = u(x) + A v(x) for the above process.
    [v, success] = simple_HJBE_discretized_univariate(A, x, u, rho);
    %dlmwrite(strcat(mfilename, '_2_v_output.csv'), v, 'precision', tolerances.default_csv_precision); %Uncomment to save again
    v_check = dlmread(strcat(mfilename, '_2_v_output.csv'));    
    verifyTrue(testCase,norm(v - v_check, Inf) < tolerances.test_tol, 'v value no longer matches');
    verifyTrue(testCase, success==true, 'unsuccesful');
end
% 
%This will try to jointly solve for the KFE and the value function with linear least squares.
function joint_value_function_KFE_test(testCase)
    tolerances = testCase.TestData.tolerances;
    
    mu_x = @(x) -0.01 * ones(numel(x),1);
    sigma_bar = 0.1;
    sigma_2_x = @(x) (sigma_bar * ones(numel(x),1)).^2;
    u_x = @(x) log(x);
    rho = 0.05;
    x_min = 0.1;
    x_max = 2;
    I = 1000;
    x = linspace(x_min, x_max, I)';
    A = discretize_univariate_diffusion(x, mu_x(x), sigma_2_x(x));     
    u = u_x(x);
    
    %Solve the simple problem: rho v(x) = u(x) + A v(x) for the above process.
    disp('Solving HJBE on its own');
    tic;
    [v, success] = simple_HJBE_discretized_univariate(A, x, u, rho);
    toc;
    %dlmwrite(strcat(mfilename, '_3_v_output.csv'), v, 'precision', tolerances.default_csv_precision); %Uncomment to save again
    v_check = dlmread(strcat(mfilename, '_3_v_output.csv'));    
    verifyTrue(testCase,norm(v - v_check, Inf) < tolerances.test_tol, 'v value no longer matches');
    verifyTrue(testCase, success==true, 'unsuccesful');
    
    %Now find the KFE with linear least squares and the default preconditioner
    settings.display = true;
    %settings.method = 'LLS';
    %settings.method = 'eigenproblem_all';
    disp('Solving stationary distribution on its own');
    tic;
    f = stationary_distribution_discretized_univariate(A, x, settings);    
    toc;
    %dlmwrite(strcat(mfilename, '_3_f_output.csv'), f, 'precision', tolerances.default_csv_precision); %Uncomment to save again
    f_check = dlmread(strcat(mfilename, '_3_f_output.csv'));    
    verifyTrue(testCase,norm(f - f_check, Inf) < tolerances.test_tol, 'f value no longer matches');
        
    %Now, solve as a joint problem with LLS
    %settings.preconditioner = 'none';
    settings.sparse = 'false';
    disp('Solving dense');
    tic;
    [v, f, success] = simple_joint_HJBE_stationary_distribution_univariate(A, x, u, rho, settings);
    verifyTrue(testCase,norm(v - v_check, Inf) < tolerances.test_tol, 'v value no longer matches');        
    verifyTrue(testCase,norm(f - f_check, Inf) < tolerances.test_tol, 'f value no longer matches');    
    toc;

    %Now try as a sparse system.
    disp('Solving sparse with no preconditioner');
    settings.sparse = 'true';   
    tic;
    [v, f, success] = simple_joint_HJBE_stationary_distribution_univariate(A, x, u, rho, settings);
    toc;
    verifyTrue(testCase,norm(v - v_check, Inf) < tolerances.test_tol, 'v value no longer matches');        
    verifyTrue(testCase,norm(f - f_check, Inf) < tolerances.test_tol, 'f value no longer matches');    

    %Now try as a sparse system with jacobi preconditioner
    disp('Solving sparse with jacobi preconditioner');
    settings.sparse = 'jacobi';    
    tic;
    [v, f, success] = simple_joint_HJBE_stationary_distribution_univariate(A, x, u, rho, settings);
    toc;
    verifyTrue(testCase,norm(v - v_check, Inf) < tolerances.test_tol, 'v value no longer matches');        
    verifyTrue(testCase,norm(f - f_check, Inf) < tolerances.test_tol, 'f value no longer matches');        
end

