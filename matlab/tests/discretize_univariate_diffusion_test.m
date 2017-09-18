%Using the unit testing framework in matlab.  See https://www.mathworks.com/help/matlab/matlab_prog/write-simple-test-case-with-functions.html
%To run tests:
% runtests %would run all of them in the current directory
% runtests('my_test') %runs just the my_test.m file
% runtests('my_test/my_function_test') %runs only `my_function_test function in `my_test'.

function tests = discretize_univariate_diffusion_test
    tests = functiontests(localfunctions);
end

%This is run at the beginning of the test.  Not required.
function setupOnce(testCase)
    addpath('../lib/');
    testCase.TestData.tolerances.test_tol = 1e-9;    
    testCase.TestData.tolerances.default_csv_precision = '%.10f'; %Should be higher precision than tolerances.test_tol
end

function zero_drift_uniform_grid_test(testCase)%Simple and small with zero drift with uniform grid
    tolerances = testCase.TestData.tolerances;
    
    mu_x = @(x) zeros(numel(x),1);
    sigma_bar = 0.1;
    sigma_2_x = @(x) (sigma_bar*x).^2;
    x_min = 0.01;
    x_max = 1;
    I = 5;
    x = linspace(x_min, x_max, I)';
    A = discretize_univariate_diffusion(x, mu_x(x), sigma_2_x(x));     

    %dlmwrite(strcat(mfilename, '_1_A_output.csv'), full(A), 'precision', tolerances.default_csv_precision); %Uncomment to save again
    A_check = dlmread(strcat(mfilename, '_1_A_output.csv'));    
    
    verifyTrue(testCase,norm(A - A_check, Inf) < tolerances.test_tol, 'A value no longer matches');
    
    %The following are worth testing for almost every matrix in the test suit.
    verifyTrue(testCase,is_stochastic_matrix(testCase, A), 'Intensity matrix rows do not sum to 0');
    verifyTrue(testCase,is_negative_diagonal(testCase, A), 'Intensity Matrix diagonal has positive elements');
    verifyTrue(testCase,isbanded(A,1,1), 'Intensity Matrix is not tridiagonal');
    verifyTrue(testCase,is_negative_definite(testCase, A), 'Intensity Matrix is not positive definite');
end    
    
function large_zero_drift_uniform_grid_test(testCase)
    tolerances = testCase.TestData.tolerances;
    
    mu_x = @(x) zeros(numel(x),1);
    sigma_bar = 0.1;
    sigma_2_x = @(x) (sigma_bar*x).^2;
    x_min = 0.01;
    x_max = 1;
    I = 1000000; %million X million matrix, but sparse.
    x = linspace(x_min, x_max, I)';
    A = discretize_univariate_diffusion(x, mu_x(x), sigma_2_x(x));     
   
    %The following are worth testing for almost every matrix in the test suit.
    verifyTrue(testCase, (nnz(A) == 2999998), 'Number of non-zero values is wrong'); %Should have about 3million non-zeros.  Tridiagonal.
    verifyTrue(testCase,is_stochastic_matrix(testCase, A), 'Intensity matrix rows do not sum to 0');
    verifyTrue(testCase,is_negative_diagonal(testCase, A), 'Intensity Matrix diagonal has positive elements');
    verifyTrue(testCase,isbanded(A,1,1), 'Intensity Matrix is not tridiagonal');
end    

function negative_drift_uniform_grid_test(testCase)
    tolerances = testCase.TestData.tolerances;
    
    mu_x = @(x) zeros(numel(x),1);
    sigma_bar = 0.1;
    sigma_2_x = @(x) (sigma_bar*x).^2;
    x_min = 0.01;
    x_max = 1;
    I = 1001;
    x = linspace(x_min, x_max, I)';
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


% TODO: Should add these checks on equations.
% 	% Variation 1: construct the A assuming that mu < 0 (i.e., the direction of the finite differences is known a-priori)
% 	X_var1 = - mu/Delta + sigma_2/(2*Delta_2);
% 	Y_var1 = mu/Delta - sigma_2/Delta_2;
% 	Z_var1 = sigma_2/(2*Delta_2);
% 	A_var1 = spdiags(Y_var1, 0, I, I) + spdiags(X(2:I), -1, I, I) + spdiags([0;Z(1:I-1)], 1, I, I);
% 	A_var1(I,I)= Y(I) + sigma_2(I)/(2*Delta_2);
% 	% Variation 2: construct the A with a for loop, essentially adding in each row as an equation.  Map to exact formulas in a latex document.
% 	S = zeros(I+2, I+2);
% 	for i = 1: I
% 	  x_i = -mu(i)/Delta + sigma_2(i)/(2*Delta_2);  % equation (8)
% 	  y_i = mu(i)/Delta - sigma_2(i)/Delta_2;       % equation (9)
% 	  z_i = sigma_2(i)/(2*Delta_2);              % equation (10)
% 	  S(i+1, i) = x_i;
% 	  S(i+1, i+1) = y_i;
% 	  S(i+1, i+2) = z_i;
% 	end
% 	S(I+1, I+1) = mu(I)/Delta - sigma_2(I)/(2*Delta_2);  %% equation (11)
% 	A_var2 = sparse(S(2: I+1, 2: I+1));
% 	



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
