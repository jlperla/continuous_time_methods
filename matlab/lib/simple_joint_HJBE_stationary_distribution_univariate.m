%Takes the discretized operator A, the grid x, and finds the stationary distribution f.
function [v, f, success] = simple_joint_HJBE_stationary_distribution_univariate(A, x, u, rho, settings)
    I = length(x);
    Delta = x(2) - x(1);
    if nargin < 5
       settings.default = true; %Just creates as required.
    end

    if(~isfield(settings, 'normalization'))
       settings.normalization = 'sum'; %The only method supported right now is a direct sum
       %Otherwise, could consider better quadrature methods using x such as trapezoidal or simpsons rule.
    end
    if ~isfield(settings, 'display')
       settings.display = false; %Tolerance
    end

     if(isfield(settings, 'max_iterations')) %OTherwise use the default
        max_iterations = settings.max_iterations; %Number of iterations
     else
        max_iterations = 10*I; %The default is too small for our needs
     end
     if(isfield(settings, 'tolerance')) %OTherwise use the default
        tolerance = settings.tolerance; %Number of iterations
     else
        tolerance = []; %Empty tells it to use default
     end 
     if(~isfield(settings, 'preconditioner'))
         settings.preconditioner = 'none'; %Default is no preconditioner.
     end
     if(~isfield(settings, 'sparse'))
         settings.sparse = true;
     end


     if(isfield(settings, 'initial_guess'))         
         initial_guess = settings.initial_guess;
     else
         initial_guess = [];
     end
     
     
    %Create the joint system
    y = sparse([Delta * u; sparse(I,1);Delta * 1]); %(41)
    X = [(Delta * rho * speye(I) - A) sparse(I,I); sparse(I,I) A'; sparse(1,I) Delta*ones(1,I)]; %(42).  Only supporting simple sum.

    if(settings.sparse == true)
         if(strcmp(settings.preconditioner,'jacobi'))
             preconditioner = diag(diag(X)); %Jacobi preconditioner is easy to calculate.  Helps a little
         elseif(strcmp(settings.preconditioner,'incomplete_LU'))
             %Matter if it is negative or positive?
             [L,U] = ilu(X(1:end-1,:))
              preconditioner = L;
         elseif(strcmp(settings.preconditioner,'none'))
             preconditioner = [];
         else
             assert(false, 'unsupported preconditioner');
         end

        [val,flag,relres,iter] = lsqr(X, y, tolerance, max_iterations, preconditioner, [], initial_guess); %Linear least squares.  Note tolerance changes with I
        if(flag==0)
            success = true;
            %Extracts the solution.
            v = val(1:I);
            f = val(I+1:end);
        else
            if(settings.display)
                disp('Failure to converge: flag and residual');
                [flag, relres]
            end
            success = false;
            f = NaN;
            v = NaN;
        end
    else %Otherwise solve as a dense system
       val = full(X) \ full(y);
       v = val(1:I);
       f = val(I+1:end);
       success = true;
    end
end	
