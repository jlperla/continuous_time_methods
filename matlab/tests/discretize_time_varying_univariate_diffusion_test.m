%Using the unit testing framework in matlab.  See https://www.mathworks.com/help/matlab/matlab_prog/write-simple-test-case-with-functions.html
%To run tests:
% runtests %would run all of them in the current directory
% runtests('my_test') %runs just the my_test.m file
% runtests('my_test/my_function_test') %runs only `my_function_test function in `my_test'.

function tests = discretize_time_varying_univariate_diffusion_test
    tests = functiontests(localfunctions);
end

%This is run at the beginning of the test.  Not required.
function setupOnce(testCase)
    addpath('../lib/');
    testCase.TestData.tolerances.test_tol = 1e-9;    
    testCase.TestData.tolerances.default_csv_precision = '%.10f'; %Should be higher precision than tolerances.test_tol
end

function basic_test(testCase)
    %This will create a non-time-varying setup, and should show that  
    tolerances = testCase.TestData.tolerances;
    
    mu_tx = @(t, x) -0.1 + t + .1*x;%cludge if constant since bsxfun gets confused otherwise
    sigma_bar = 0.1;
    sigma_2_tx = @(t, x) (sigma_bar*x).^2;
    
    %Grid
    x_min = 0.01;
    x_max = 1;
    I = 5;
	t_min = 0.0;
	t_max = 10.0;
	N = 4;
    x = linspace(x_min, x_max, I)';
	t = linspace(t_min, t_max, N)';
    
   [t_grid, x_grid] = meshgrid(t,x); %Generates permutations (stacked by t first, as we want) could look at: [t_grid(:) x_grid(:)]
   state_permutations = [t_grid(:) x_grid(:)];
  
   %should check that this does the correct permuations in order and calls the underlying function.
   mu = bsxfun(mu_tx, state_permutations(:,1), state_permutations(:,2)); %applies binary function to these, and remains in the correct stacked order.
   sigma_2 = bsxfun(sigma_2_tx, state_permutations(:,1), state_permutations(:,2)); %applies binary function to these, and remains in the correct stacked order.
  
   %Discretize the operator
   [A, Delta_p, Delta_m, h_p, h_m] = discretize_time_varying_univariate_diffusion(t, x, mu, sigma_2);
   
end    

