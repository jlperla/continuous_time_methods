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
    %This will create a time-varying setup, and should show that  
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
   
   %dlmwrite(strcat(mfilename, '_1_A_output.csv'), full(A), 'precision', tolerances.default_csv_precision); %Uncomment to save again
   
   A_check = dlmread(strcat(mfilename, '_1_A_output.csv'));    
    
   verifyTrue(testCase,norm(A - A_check, Inf) < tolerances.test_tol, 'A value no longer matches');
end    

function non_time_varying_test(testCase)
    %This will create a non-time-varying setup, and should show that  
    tolerances = testCase.TestData.tolerances;
    
    % This test use mu(x)=mu*x;sigma(x)^2=sigma^2*x^2 and u(x)=exp(x)
    mu_tx = @(t,x) -0.01 * x+0*t;
    mu_x  = @(x) -0.01 * x;
    sigma_bar = 0.1;
    sigma_2_tx = @(t,x) (sigma_bar*x).^2+0*t;
    sigma_2_x  = @(x) (sigma_bar*x).^2;
    u_tx = @(t,x) exp(x)+t*0;
    u_x  = @(x) exp(x);
    rho = 0.05;
    x_min = 0.1;
    x_max = 3;
    t_min = 0.0;
    t_max = 10.0;
    I = 1500;
    N = 10;
    x = linspace(x_min, x_max, I)';
    t = linspace(t_min,t_max, N)';
    
    % generate permutations of x and t, tacking t first then x
    
    [t_grid, x_grid] = meshgrid(t,x); %Generates permutations (stacked by t first, as we want) could look at: [t_grid(:) x_grid(:)]
    state_permutations = [t_grid(:) x_grid(:)];

    mu = bsxfun(mu_tx, state_permutations(:,1), state_permutations(:,2)); %applies binary function to these, and remains in the correct stacked order.
    sigma_2 = bsxfun(sigma_2_tx, state_permutations(:,1), state_permutations(:,2)); %applies binary function to these, and remains in the correct stacked order.
  
    % Solve for A_n and v_n for non time-varying method, use as check
    %A_n = discretize_univariate_diffusion(x, mu_x(x), sigma_2_x(x));
    %dlmwrite(strcat(mfilename, '_2_An_output.csv'), full(A_n), 'precision', tolerances.default_csv_precision); %Uncomment to save again
    A_check = dlmread(strcat(mfilename, '_2_An_output.csv'));
    u_n=u_x(x);
    
    %Solve the simple problem: rho v(x) = u(x) + A v(x) for the non-time
    %varying sprocess
    %v_n = simple_HJBE_discretized_univariate(A_n, x, u_n, rho);
    %dlmwrite(strcat(mfilename, '_2_vn_output.csv'), v_n, 'precision', tolerances.default_csv_precision); %Uncomment to save again
    v_check = dlmread(strcat(mfilename, '_2_vn_output.csv'));    
    
    % Solve for A and v for time-varying method
    [A, Delta_p, Delta_m, h_p, h_m] = discretize_time_varying_univariate_diffusion(t, x, mu, sigma_2);
    
    u = bsxfun(u_tx, state_permutations(:,1), state_permutations(:,2)); % sp(:,1) is time sp(:,2) is x
    
    [v,success] = simple_HJBE_discretized_univariate(A, state_permutations(:,1), u, rho); % state_perm need to be in size N*I
    % dlmwrite(strcat(mfilename, '_2_v_output.csv'), v, 'precision', tolerances.default_csv_precision); %Uncomment to save again
    v_check_1 = dlmread(strcat(mfilename, '_2_v_output.csv'));
    
    % tests for 1. A is same as A_n; 2. A is same across time; 3. v is same
    % as v_n ; 4. v is same across time
    
    % initial test: make sure the v of time varying result didn't change.
    verifyTrue(testCase,norm(v - v_check_1, Inf) < 1e-9, 'v time-varying result no longer matches');
    
    % this checks the matrix A1 in time varying case same as A in
    % non-time-varying case
    verifyTrue(testCase,norm(A(1:I,1:I) + A(1:I,I+1:2*I) - A_check, Inf) < 1e-5, 'A_1 not match with non-time-varying result'); 
    
    % This checks the matrix A1 in time varying case same as A_n, n
    % randomly drawn
    luck=max(floor(rand*N),2);
    indx=(luck-1)*I+1;
    verifyTrue(testCase,norm(A(1:I,1:I) + A(1:I,I+1:2*I) - A(indx:indx+I-1,indx:indx+I-1) - A(indx:indx+I-1,indx+I:indx+2*I-1) , Inf) < 1e-9, 'A_1 not match with A_n'); 
    %verifyTrue(testCase, success==true, 'unsuccesful');
   
    % This checks the value function of t=0 is same as v in
    % non-time-varying result
    verifyTrue(testCase,norm(v(1:I) - v_check, Inf) < 1e-5, 'v_1 not match with non-time-varying result');
    
    % This checks the v_1 in time varying case same as v_n
    luck=max(floor(rand*N),2);
    indx=(luck-1)*I+1;
    verifyTrue(testCase,norm(v(1:I) - v(indx:indx+I-1), Inf) < 1e-9, 'v_1 not match with v_n');

end    
