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

    %%dlmwrite(strcat(mfilename, '_1_A_output.csv'), full(A), 'precision', tolerances.default_csv_precision); %Uncomment to save again
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
   
    %The following are worth testing for almost every matrix in the test suite.
    verifyTrue(testCase, (nnz(A) == 2999998), 'Number of non-zero values is wrong'); %Should have about 3million non-zeros.  Tridiagonal.
    verifyTrue(testCase,is_stochastic_matrix(testCase, A), 'Intensity matrix rows do not sum to 0');
    verifyTrue(testCase,is_negative_diagonal(testCase, A), 'Intensity Matrix diagonal has positive elements');
    verifyTrue(testCase,isbanded(A,1,1), 'Intensity Matrix is not tridiagonal');
end    

function negative_drift_uniform_grid_test(testCase)
    tolerances = testCase.TestData.tolerances;
    
    mu_x = @(x) ones(numel(x),1) * (-1);
    sigma_bar = 0.1;
    sigma_2_x = @(x) (sigma_bar*x).^2;
    x_min = 0.01;
    x_max = 1;
    I = 1001;
    x = linspace(x_min, x_max, I)';
    A = discretize_univariate_diffusion(x, mu_x(x), sigma_2_x(x));     

    
    %To save the file again, can uncomment this.
    [indices_i, indices_j, values_ij] = find(A); %Uncomment to save again
    %%dlmwrite(strcat(mfilename, '_3_A_output.csv'), [indices_i indices_j values_ij], 'precision', tolerances.default_csv_precision); %Uncomment to save again
    
    %Load and check against the sparse matrix file.
    A_check = spconvert(dlmread(strcat(mfilename, '_3_A_output.csv')));
    verifyTrue(testCase,norm(A - A_check, Inf) < tolerances.test_tol, 'A value no longer matches');
    
    %The following are worth testing for almost every matrix in the test suit.
    verifyTrue(testCase,is_stochastic_matrix(testCase, A), 'Intensity matrix rows do not sum to 0');
    verifyTrue(testCase,is_negative_diagonal(testCase, A), 'Intensity Matrix diagonal has positive elements');
    verifyTrue(testCase,isbanded(A,1,1), 'Intensity Matrix is not tridiagonal');
end


function x_min_is_less_than_zero_test(testCase)    % x_min < 0
    tolerances = testCase.TestData.tolerances;
    
    mu_x = @(x) zeros(numel(x),1);
    sigma_bar = 0.1;
    sigma_2_x = @(x) (sigma_bar*x).^2;
    x_min = -0.49;
    x_max = 0.5;
    I = 5;
    x = linspace(x_min, x_max, I)';
    A = discretize_univariate_diffusion(x, mu_x(x), sigma_2_x(x)); 

    %%dlmwrite(strcat(mfilename, '_4_A_output.csv'), full(A), 'precision', tolerances.default_csv_precision); %Uncomment to save again
    A_check = dlmread(strcat(mfilename, '_4_A_output.csv'));
    
    verifyTrue(testCase,norm(A - A_check, Inf) < tolerances.test_tol, 'A value no longer matches');
    
    %The following are worth testing for almost every matrix in the test suit.
    verifyTrue(testCase,is_stochastic_matrix(testCase, A), 'Intensity matrix rows do not sum to 0');
    verifyTrue(testCase,is_negative_diagonal(testCase, A), 'Intensity Matrix diagonal has positive elements');
    verifyTrue(testCase,isbanded(A,1,1), 'Intensity Matrix is not tridiagonal');
end


%Don't remember what this was for, but doesn't apply now that we have Delta back in denominator.
% function rescale_x_min_and_x_max_test(testCase)   % Change in the scale of x_min and x_max for given I
%     tolerances = testCase.TestData.tolerances;
%     
%     mu_x = @(x) zeros(numel(x),1);
%     sigma_bar = 0.1;
%     sigma_2_x = @(x) (sigma_bar*x).^2;
%     x_min = 0.01;
%     x_max = 1;
%     I = 1001; 
%     scaler = 6;   % Just pick 6 as a random scaler to rescale x.
%     x = linspace(x_min, x_max, I)';
%     x_rescale = scaler * x;  
%     A = discretize_univariate_diffusion(x, mu_x(x), sigma_2_x(x)); 
%     A_rescale = discretize_univariate_diffusion(x_rescale, mu_x(x_rescale), sigma_2_x(x_rescale)) / scaler; 
%    
%     verifyTrue(testCase, norm(A - A_rescale, Inf) < tolerances.test_tol, 'A value no longer matches');
%     
%     %The following are worth testing for almost every matrix in the test suit.
%     verifyTrue(testCase,is_stochastic_matrix(testCase, A), 'Intensity matrix rows do not sum to 0');
%     verifyTrue(testCase,is_negative_diagonal(testCase, A), 'Intensity Matrix diagonal has positive elements');
%     verifyTrue(testCase,isbanded(A,1,1), 'Intensity Matrix is not tridiagonal');    
% end


function construction_test(testCase)  %Use variations of construction with mu<0
    tolerances = testCase.TestData.tolerances;
    
    mu_x = @(x) ones(numel(x),1) * (-2);  % A random mu, which is less than 0
    sigma_bar = 0.1;
    sigma_2_x = @(x) (sigma_bar*x).^2;
    x_min = 0.01;
    x_max = 1;
    I = 1001; 
    x = linspace(x_min, x_max, I)';
    Delta = x(2) - x(1);
    Delta_2 = Delta^2;
    A = discretize_univariate_diffusion(x, mu_x(x), sigma_2_x(x));
    mu = mu_x(x);
    sigma_2 = sigma_2_x(x);
    
    % Variation 1: construct the A assuming that mu < 0 (i.e., the direction of the finite differences is known a-priori)
    X_var1 = - mu/Delta + sigma_2/(2*Delta_2);
 	Y_var1 = mu/Delta - sigma_2/Delta_2;
 	Z_var1 = sigma_2/(2*Delta_2);
 	A_var1 = spdiags(Y_var1, 0, I, I) + spdiags(X_var1(2:I), -1, I, I) + spdiags([0;Z_var1(1:I-1)], 1, I, I);
    A_var1(1, 1) = Y_var1(1) + X_var1(1);
 	A_var1(I,I)= Y_var1(I) + sigma_2(I)/(2*Delta_2);
    
    % Variation 2: construct the A with a for loop, essentially adding in each row as an equation.  Map to exact formulas in a latex document.
 	S = zeros(I, I+2);
 	for i = 1: I
 	  x_i = -mu(i)/Delta + sigma_2(i)/(2*Delta_2);  
 	  y_i = mu(i)/Delta - sigma_2(i)/Delta_2;       
 	  z_i = sigma_2(i)/(2*Delta_2);              
 	  S(i, i) = x_i;
 	  S(i, i+1) = y_i;
 	  S(i, i+2) = z_i;
    end
    S(1, 2) = S(1, 2) + S(1, 1);
	S(I, I+1) = mu(I)/Delta - sigma_2(I)/(2*Delta_2);  
 	A_var2 = sparse(S(:, 2: I+1));
 	
    
    verifyTrue(testCase, norm(A - A_var1, Inf) < tolerances.test_tol, 'A is different from A_var1');
    verifyTrue(testCase, norm(A - A_var2, Inf) < tolerances.test_tol, 'A is different from A_var2');
    
    %The following are worth testing for almost every matrix in the test suit.
    verifyTrue(testCase,is_stochastic_matrix(testCase, A), 'Intensity matrix rows do not sum to 0');
    verifyTrue(testCase,is_negative_diagonal(testCase, A), 'Intensity Matrix diagonal has positive elements');
    verifyTrue(testCase,isbanded(A,1,1), 'Intensity Matrix is not tridiagonal'); 
end


function monotonically_increasing_mu_test(testCase)    % mu is monotonically increasing in x with mu(x_min)<0 and mu(x_max)>0
    tolerances = testCase.TestData.tolerances;
    
    mu_x = @(x) (x - 0.5);
    sigma_bar = 0.1;
    sigma_2_x = @(x) (sigma_bar*x).^2;
    x_min = 0.01;
    x_max = 1;
    I = 5;
    x = linspace(x_min, x_max, I)';
    A = discretize_univariate_diffusion(x, mu_x(x), sigma_2_x(x));     
  
    %dlmwrite(strcat(mfilename, '_7_A_output.csv'), full(A), 'precision', tolerances.default_csv_precision); %Uncomment to save again
    A_check = dlmread(strcat(mfilename, '_7_A_output.csv'));
    
   
    verifyTrue(testCase,norm(A - A_check, Inf) < tolerances.test_tol, 'A value no longer matches');
    
    %The following are worth testing for almost every matrix in the test suit.
    verifyTrue(testCase,is_stochastic_matrix(testCase, A), 'Intensity matrix rows do not sum to 0');
    verifyTrue(testCase,is_negative_diagonal(testCase, A), 'Intensity Matrix diagonal has positive elements');
    verifyTrue(testCase,isbanded(A,1,1), 'Intensity Matrix is not tridiagonal');
    verifyTrue(testCase,is_negative_definite(testCase, A), 'Intensity Matrix is not positive definite');
end    


function monotonically_decreasing_mu_test(testCase)    % mu is monotonically decreasing in x with mu(x_min)>0 and mu(x_max)<0
    tolerances = testCase.TestData.tolerances;
    
    mu_x = @(x) (x - 0.5) * (-1);
    sigma_bar = 0.1;
    sigma_2_x = @(x) (sigma_bar*x).^2;
    x_min = 0.01;
    x_max = 1;
    I = 5;
    x = linspace(x_min, x_max, I)';
    A = discretize_univariate_diffusion(x, mu_x(x), sigma_2_x(x));     

    %dlmwrite(strcat(mfilename, '_8_A_output.csv'), full(A), 'precision', tolerances.default_csv_precision); %Uncomment to save again
    A_check = dlmread(strcat(mfilename, '_8_A_output.csv'));
    
    verifyTrue(testCase,norm(A - A_check, Inf) < tolerances.test_tol, 'A value no longer matches');
    
    %The following are worth testing for almost every matrix in the test suit.
    verifyTrue(testCase,is_stochastic_matrix(testCase, A), 'Intensity matrix rows do not sum to 0');
    verifyTrue(testCase,is_negative_diagonal(testCase, A), 'Intensity Matrix diagonal has positive elements');
    verifyTrue(testCase,isbanded(A,1,1), 'Intensity Matrix is not tridiagonal');
    verifyTrue(testCase,is_negative_definite(testCase, A), 'Intensity Matrix is not positive definite');
end 


function concave_mu_test(testCase)     % mu is concave in x with mu(x_min)<0, mu(x_max)<0 and mu(x)>0 for some x
    tolerances = testCase.TestData.tolerances;
    
    mu_x = @(x) (-(x - 0.5).^2 + 0.1);
    sigma_bar = 0.1;
    sigma_2_x = @(x) (sigma_bar*x).^2;
    x_min = 0.01;
    x_max = 1;
    I = 5;
    x = linspace(x_min, x_max, I)';
    A = discretize_univariate_diffusion(x, mu_x(x), sigma_2_x(x));     

    %dlmwrite(strcat(mfilename, '_9_A_output.csv'), full(A), 'precision', tolerances.default_csv_precision); %Uncomment to save again
    A_check = dlmread(strcat(mfilename, '_9_A_output.csv'));
    
    verifyTrue(testCase,norm(A - A_check, Inf) < tolerances.test_tol, 'A value no longer matches');
    
    %The following are worth testing for almost every matrix in the test suit.
    verifyTrue(testCase,is_stochastic_matrix(testCase, A), 'Intensity matrix rows do not sum to 0');
    verifyTrue(testCase,is_negative_diagonal(testCase, A), 'Intensity Matrix diagonal has positive elements');
    verifyTrue(testCase,isbanded(A,1,1), 'Intensity Matrix is not tridiagonal');
    verifyTrue(testCase,is_negative_definite(testCase, A), 'Intensity Matrix is not positive definite');
end

%Removed.  Should be checking that this throws errors, but not sure how to do it with matlab tests.
% function zero_sigma_everywhere_test(testCase)
%     I = 5;
%     %mu_x = @(x) zeros(numel(x),1);
%     mu_x = @(x) -.01 * ones(numel(x),1);
%     sigma_2_x = @(x) zeros(numel(x),1);
%     x_min = 0.01;
%     x_max = 1;
%     x = linspace(x_min, x_max, I)';
%     %A = discretize_univariate_diffusion(x, mu_x(x), sigma_2_x(x));     
%     %testCase.assertFail(@() discretize_univariate_diffusion(x, mu_x(x), sigma_2_x(x))); ???? Not working, 
% end

function zero_sigma_somewhere_test(testCase)
    tolerances = testCase.TestData.tolerances;
    
    mu_x = @(x) (x - 0.5) * (-1);
    sigma_bar = 0.1;
    sigma_2_x = @(x) (sigma_bar*x).^2;
    x_min = 0.01;
    x_max = 1;
    I = 5;
    x = linspace(x_min, x_max, I)';
    sigma_2 = sigma_2_x(x);
    sigma_2(1, 1) = 0;
    sigma_2(3, 1) = 0;
    sigma_2(5, 1) = 0;
    A = discretize_univariate_diffusion(x, mu_x(x), sigma_2);     

    %dlmwrite(strcat(mfilename, '_11_A_output.csv'), full(A), 'precision', tolerances.default_csv_precision); %Uncomment to save again
    A_check = dlmread(strcat(mfilename, '_11_A_output.csv'));
    
    verifyTrue(testCase,norm(A - A_check, Inf) < tolerances.test_tol, 'A value no longer matches');
    
    %The following are worth testing for almost every matrix in the test suit.
    verifyTrue(testCase,is_stochastic_matrix(testCase, A), 'Intensity matrix rows do not sum to 0');
    verifyTrue(testCase,is_negative_diagonal(testCase, A), 'Intensity Matrix diagonal has positive elements'); 
    verifyTrue(testCase,isbanded(A,1,1), 'Intensity Matrix is not tridiagonal');
    verifyTrue(testCase,is_negative_definite(testCase, A), 'Intensity Matrix is not positive definite');
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
