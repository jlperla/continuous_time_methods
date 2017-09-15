addpath('../lib/');
main_script_tested = 'discretize_univariate_diffusion';
default_csv_precision = '%.10f'; %Should be higher precision than test_tol
% Define an absolute tolerance for floating point comparisons
test_tol = 1e-9;

verify_stochastic_matrix = @(A) max(abs(full(sum(A,2)))) < test_tol; %Ensures that all rows sum to 0, which is essential for an intensity matrix.  Wouldn't hold without reflective barriers in this csae.
verify_negative_diagonal = @(A) max(full(diag(A))) < 0;  %intensity matrices need to have negatives along the diagonal

%% Test 1: Simple and small with zero drift with uniform grid
    mu_x = @(x) zeros(numel(x),1);
    sigma_bar = 0.1;
    sigma_2_x = @(x) (sigma_bar*x).^2;
    x_min = 0.01;
    x_max = 1;
    I = 5;
    [A, x] = discretize_univariate_diffusion_uniform_driver(mu_x, sigma_2_x, I, x_min, x_max);
    %dlmwrite(strcat(main_script_tested, '_1_A_output.csv'), full(A), 'precision', default_csv_precision); %Uncomment to save again
    A_check = dlmread(strcat(main_script_tested, '_1_A_output.csv'));    
    
    assert(norm(A - A_check, Inf) < test_tol, 'A value no longer matches');
    
    %The following are worth testing for almost every matrix in the test suit.
    assert(verify_stochastic_matrix(A), 'Intensity matrix rows do not sum to 0');
    assert(verify_negative_diagonal(A), 'Intensity Matrix diagonal has positive elements');
    assert(isbanded(A,1,1), 'Intensity Matrix is not tridiagonal');
    
    
%% Test 2: Large Matrix, no drift with uniform grid
    mu_x = @(x) zeros(numel(x),1);
    sigma_bar = 0.1;
    sigma_2_x = @(x) (sigma_bar*x).^2;
    x_min = 0.01;
    x_max = 1;
    I = 1000000; %million X million matrix, but sparse.
    [A, x] = discretize_univariate_diffusion_uniform_driver(mu_x, sigma_2_x, I, x_min, x_max);

   
    %The following are worth testing for almost every matrix in the test suit.
    assert(nnz(A) == 2999998, 'Number of non-zero values is wrong'); %Should have about 3million non-zeros.  Tridiagonal.
    assert(verify_stochastic_matrix(A), 'Intensity matrix rows do not sum to 0');
    assert(verify_negative_diagonal(A), 'Intensity Matrix diagonal has positive elements');
    assert(isbanded(A,1,1), 'Intensity Matrix is not tridiagonal');
    %No need to compare to a stored version at this point, more a quetsion of ensuring the speed doesn't collapse
    
%% Test 3: Negative drift with uniform grid
    mu_x = @(x) zeros(numel(x),1);
    sigma_bar = 0.1;
    sigma_2_x = @(x) (sigma_bar*x).^2;
    x_min = 0.01;
    x_max = 1;
    I = 1001;
    [A, x] = discretize_univariate_diffusion_uniform_driver(mu_x, sigma_2_x, I, x_min, x_max);
    
    %To save the file again, can uncomment this.
    %[indices_i, indices_j, values_ij] = find(A); %Uncomment to save again
    %dlmwrite(strcat(main_script_tested, '_3_A_output.csv'), [indices_i indices_j values_ij], 'precision', default_csv_precision); %Uncomment to save again
    
    %Load and check against the sparse matrix file.
    A_check = spconvert(dlmread(strcat(main_script_tested, '_3_A_output.csv')));
    assert(norm(A - A_check, Inf) < test_tol, 'A value no longer matches');
    
    %The following are worth testing for almost every matrix in the test suit.
    assert(verify_stochastic_matrix(A), 'Intensity matrix rows do not sum to 0');
    assert(verify_negative_diagonal(A), 'Intensity Matrix diagonal has positive elements');
    assert(isbanded(A,1,1), 'Intensity Matrix is not tridiagonal');

% TODO: Moving to a separate test file.
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

function [A, x] = discretize_univariate_diffusion_uniform_driver(mu_x, sigma_2_x, I, x_min, x_max) %This is a driver for many tests with a uniform grid.
     x = linspace(x_min, x_max, I)';
     A = discretize_univariate_diffusion(x, mu_x(x), sigma_2_x(x));     
end