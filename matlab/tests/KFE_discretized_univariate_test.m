%Using the unit testing framework in matlab.  See https://www.mathworks.com/help/matlab/matlab_prog/write-simple-test-case-with-functions.html
%To run tests:
% runtests %would run all of them in the current directory
% runtests('my_test') %runs just the my_test.m file
% runtests('my_test/my_function_test') %runs only `my_function_test function in `my_test'.

function tests = KFE_discretized_univariate_test
    tests = functiontests(localfunctions);
end

%This is run at the beginning of the test.  Not required.
function setupOnce(testCase)
    addpath('../lib/');
    testCase.TestData.tolerances.test_tol = 1e-9;    
    testCase.TestData.tolerances.lower_test_tol = 1e-6; %Too high precision for some tests with big matrices
    testCase.TestData.tolerances.default_csv_precision = '%.10f'; %Should be higher precision than tolerances.test_tol
end

function [A, x] = baseline_negative_drift_discretization(I, testCase) %Used by other test cases as a baseline.
    mu_x = @(x) -0.01 * x;
    sigma_bar = 0.1;
    sigma_2_x = @(x) (sigma_bar*x).^2;
    x_min = 1;
    x_max = 2;
    x = linspace(x_min, x_max, I)';
    A = discretize_univariate_diffusion(x, mu_x(x), sigma_2_x(x));
end

function small_LLS_vs_eigenvalue_test(testCase) %Simple baseline check.
    tolerances = testCase.TestData.tolerances;
    
    I = 20; %Small matrix.
    [A, x] = baseline_negative_drift_discretization(I, testCase);

    f = stationary_distribution_discretized_univariate(A, x);
%    dlmwrite(strcat(mfilename, '_1_f_output.csv'), f, 'precision', tolerances.default_csv_precision); %Uncomment to save again
    f_check = dlmread(strcat(mfilename, '_1_f_output.csv'));    
    verifyTrue(testCase,norm(f - f_check, Inf) < tolerances.test_tol, 'f value no longer matches');
    
    settings.method = 'LLS'; %Now using the LLS method
    f_lls = stationary_distribution_discretized_univariate(A,x, settings);
    verifyTrue(testCase,norm(f_lls - f_check, Inf) < tolerances.lower_test_tol, 'f value no longer matches');    
    
    %With a perfect initial guess!
    settings.initial_guess = f_check;
    f_lls = stationary_distribution_discretized_univariate(A,x, settings);
    verifyTrue(testCase,norm(f_lls - f_check, Inf) < tolerances.lower_test_tol, 'f value no longer matches');    
    
end

function medium_LLS_vs_eigenvalue_test(testCase) %Larger system using default options, etc.
    tolerances = testCase.TestData.tolerances;
    
    I = 1000; %Larger grid
    settings.print_level = 1;    
    settings.method = 'eigenproblem'; %Will try to only find the appropriate eigenvector, which can be faster.
    settings.num_basis_vectors = 100; %In general, need to tweak this to use the default eigenvector approach for large I
    [A, x] = baseline_negative_drift_discretization(I, testCase);

    tic;
    f = stationary_distribution_discretized_univariate(A, x, settings);
    toc;
    %dlmwrite(strcat(mfilename, '_2_f_output.csv'), f, 'precision', tolerances.default_csv_precision); %Uncomment to save again
    f_check = dlmread(strcat(mfilename, '_2_f_output.csv'));    
    verifyTrue(testCase,norm(f - f_check, Inf) < tolerances.test_tol, 'f value no longer matches');
    
    settings.method = 'LLS'; %Now using the LLS method with the default pre-conditioner
    tic;    
    f_lls = stationary_distribution_discretized_univariate(A, x, settings);
    toc;
    verifyTrue(testCase,norm(f_lls - f_check, Inf) < tolerances.lower_test_tol, 'f value no longer matches');    
    
    settings.method = 'eigenproblem_all'; %Now using the eigenvalues, but calculating all
    tic;    
    f_eigen_all = stationary_distribution_discretized_univariate(A, x, settings);
    toc;
    verifyTrue(testCase,norm(f_eigen_all - f_check, Inf) < tolerances.lower_test_tol, 'f value no longer matches');       
end

function medium_preconditioner_test(testCase) %tests of the various preconditioners.  Some not worth much.
    tolerances = testCase.TestData.tolerances;
    
    I = 1000; %Larger grid
    settings.print_level = 1;    
    [A, x] = baseline_negative_drift_discretization(I, testCase);

    %Use the eigenproblem approach as a comparison
    disp('Eigenvalue based solution');
    tic;
    f = stationary_distribution_discretized_univariate(A, x, settings);
    toc;
    %dlmwrite(strcat(mfilename, '_3_f_output.csv'), f, 'precision', tolerances.default_csv_precision); %Uncomment to save again
    f_check = dlmread(strcat(mfilename, '_3_f_output.csv'));    
    verifyTrue(testCase,norm(f - f_check, Inf) < tolerances.test_tol, 'f value no longer matches');

    %Setup basic LLS
    settings.method = 'LLS'; %Now using the LLS method with the default pre-conditioner    
    settings.max_iterations = 50000; %Only really need this many to test the no preconditioner version.
    settings.tolerance = 1E-8;
    settings.print_level = 1;

    %Use LLS with no preconditioners
    tic;
    disp('LLS no preconditioner');
    settings.preconditioner = 'none';
    f = stationary_distribution_discretized_univariate(A, x, settings);
    toc;
    verifyTrue(testCase,norm(f - f_check, Inf) < tolerances.test_tol, 'f value no longer matches');

    tic;
    disp('LLS with incomplete LU preconditioner and perfect initial guess');
    settings.preconditioner = 'incomplete_LU';
    settings.initial_guess = f_check;
    f = stationary_distribution_discretized_univariate(A, x, settings);
    toc;
    verifyTrue(testCase,norm(f - f_check, Inf) < tolerances.test_tol, 'f value no longer matches');        
    
    tic;
    disp('LLS incomplete LU preconditioner');
    settings.preconditioner = 'incomplete_LU';
    f = stationary_distribution_discretized_univariate(A, x, settings);
    toc;
    verifyTrue(testCase,norm(f - f_check, Inf) < tolerances.test_tol, 'f value no longer matches');    
    
    tic;
    disp('LLS jacobi preconditioner');
    settings.preconditioner = 'jacobi';
    f = stationary_distribution_discretized_univariate(A, x, settings);
    toc;
    verifyTrue(testCase,norm(f - f_check, Inf) < tolerances.test_tol, 'f value no longer matches');      
    
    tic;
    disp('LLS incomplete_cholesky preconditioner');
    settings.preconditioner = 'incomplete_cholesky';
    f = stationary_distribution_discretized_univariate(A, x, settings);
    toc;
    verifyTrue(testCase,norm(f - f_check, Inf) < tolerances.test_tol, 'f value no longer matches');      
    
    tic;
    disp('LLS with incomplete LU preconditioner and mediocre initial guess');
    settings.preconditioner = 'incomplete_LU';
    settings.initial_guess = linspace(.005,0,I)';
    f = stationary_distribution_discretized_univariate(A, x, settings);
    toc;
    verifyTrue(testCase,norm(f - f_check, Inf) < tolerances.test_tol, 'f value no longer matches');        
end

