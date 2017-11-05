%Using the unit testing framework in matlab.  See https://www.mathworks.com/help/matlab/matlab_prog/write-simple-test-case-with-functions.html
%To run tests:
% runtests %would run all of them in the current directory
% runtests('my_test') %runs just the my_test.m file
% runtests('my_test/my_function_test') %runs only `my_function_test function in `my_test'.

function tests = discretize_nonuniform_univariate_diffusion_test
    tests = functiontests(localfunctions);
end

%This is run at the beginning of the test.  Not required.
function setupOnce(testCase)
    addpath('../lib/');
    testCase.TestData.tolerances.test_tol = 1e-9;    
    testCase.TestData.tolerances.default_csv_precision = '%.10f'; %Should be higher precision than tolerances.test_tol
end

function zero_drift_test(testCase)%Simple and small with zero drift with uniform grid
    tolerances = testCase.TestData.tolerances;
    
    mu_x = @(x) zeros(numel(x),1);
    sigma_bar = 0.1;
    sigma_2_x = @(x) (sigma_bar*x).^2;
    x_min = 0.01;
    x_max = 1;
    I = 5;
    x = logspace(log10(x_min),log10(x_max),I)';
    [A, Delta_p, Delta_m] = discretize_univariate_diffusion(x, mu_x(x), sigma_2_x(x));     

    %dlmwrite(strcat(mfilename, '_1_A_output.csv'), full(A), 'precision', tolerances.default_csv_precision); %Uncomment to save again
    A_check = dlmread(strcat(mfilename, '_1_A_output.csv'));    
    
    verifyTrue(testCase,norm(A - A_check, Inf) < tolerances.test_tol, 'A value no longer matches');
    
    %The following are worth testing for almost every matrix in the test suit.
    verifyTrue(testCase,is_stochastic_matrix(testCase, A), 'Intensity matrix rows do not sum to 0');
    verifyTrue(testCase,is_negative_diagonal(testCase, A), 'Intensity Matrix diagonal has positive elements');
    verifyTrue(testCase,isbanded(A,1,1), 'Intensity Matrix is not tridiagonal');
    verifyTrue(testCase,is_negative_definite(testCase, A), 'Intensity Matrix is not positive definite');
end    
function negative_drift_uniform_grid_test(testCase)
    tolerances = testCase.TestData.tolerances;
    
    mu_x = @(x) ones(numel(x),1) * (-1);
    sigma_bar = 0.1;
    sigma_2_x = @(x) (sigma_bar*x).^2;
    x_min = 0.01;
    x_max = 1;
    I = 1001;
    x = logspace(log10(x_min),log10(x_max),I)';
    A = discretize_univariate_diffusion(x, mu_x(x), sigma_2_x(x));     

    
    %To save the file again, can uncomment this.
    %[indices_i, indices_j, values_ij] = find(A); %Uncomment to save again
    %dlmwrite(strcat(mfilename, '_3_A_output.csv'), [indices_i indices_j values_ij], 'precision', tolerances.default_csv_precision); %Uncomment to save again
    
    %Load and check against the sparse matrix file.
    A_check = spconvert(dlmread(strcat(mfilename, '_3_A_output.csv')));
    verifyTrue(testCase,norm(A - A_check, Inf) < tolerances.test_tol, 'A value no longer matches');
    
    %The following are worth testing for almost every matrix in the test suit.
    verifyTrue(testCase,is_stochastic_matrix(testCase, A), 'Intensity matrix rows do not sum to 0');
    verifyTrue(testCase,is_negative_diagonal(testCase, A), 'Intensity Matrix diagonal has positive elements');
    verifyTrue(testCase,isbanded(A,1,1), 'Intensity Matrix is not tridiagonal');
end

%These are utility functions for testing returned matrices.
function result = is_stochastic_matrix(testCase, A)
    result = (max(abs(full(sum(A,2)))) < testCase.TestData.tolerances.test_tol);
end

function result = is_negative_diagonal(testCase, A)
    result = (max(full(diag(A))) < 0);
end

function result = is_negative_definite(testCase, A)
    result =all(eig(full(A)) < testCase.TestData.tolerances.test_tol);
end
