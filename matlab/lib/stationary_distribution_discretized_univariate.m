%Takes the discretized operator A, the grid x, and finds the stationary distribution f.
function [f, success] = stationary_distribution_discretized_univariate(A, x, settings)
   I = length(x);
   if nargin < 3
       settings.default = true; %Just creates as required.
   end
   
   %TODO: Consider adding in a 'dense' option for small matrices.  
   
   if(~isfield(settings, 'method'))
        settings.method = 'eigenproblem_all';
   end
   if(~isfield(settings, 'normalization'))
       settings.normalization = 'sum'; %The only method supported right now is a direct sum
       %Otherwise, could consider better quadrature methods using x such as trapezoidal or simpsons rule.
   end
   if ~isfield(settings, 'display')
       settings.display = false; %Tolerance
   end
   assert(I == size(A,1) && I == size(A,2)); %Make sure sizes match
   
   if(strcmp(settings.method, 'eigenproblem')) %Will use sparsity
        opts.isreal = true;
        if(isfield(settings, 'num_basis_vectors')) %Otherwise use the default
           opts.p = settings.num_basis_vectors; %Number of Lanczos basis vectors.  Need to increase often
        end
        if(isfield(settings, 'max_iterations')) %OTherwise use the default
          opts.maxit = settings.max_iterations; %Number of iterations
        end
        [V, D, flag] = eigs(A',1,'sm', opts);%The eigenvalue with the smallest magnitude should be the zero eigenvalue
        if((flag ~= 0) || (abs(D - 0.0) > 1E-9)) %The 'sm' one is hopefully the zero, but maybe not if there are convergence issues.  Also, the algorithm may simply not converge.
            if(settings.display)
                disp('The eigenvalue is not zero or did not converge.  Try increasing the num_basis_vectors or max_iterations.  Otherwise, consider eigenproblem_all');
            end
            success = false;
            f = NaN;
            return;
        end
        f = V / sum(V); %normalize to sum to 1.  Could add other normalizations using the grid 'x' depending on settings.normalization 
        success = true;
   elseif(strcmp(settings.method, 'eigenproblem_all')) %Will use sparsity but computes all of the eigenvaluse/eigenvectors.  Use if `eigenproblem' didn't work.
        opts.isreal = true;
        if(isfield(settings, 'num_basis_vectors')) %OTherwise use the default
           opts.p = settings.num_basis_vectors; %Number of Lanczos basis vectors.
        end
         if(isfield(settings, 'max_iterations')) %OTherwise use the default
          opts.maxit = settings.max_iterations; %Number of iterations
        end
        [V,D] = eigs(A', I, 'sm',opts); %Gets all of the eigenvalues and eigenvectors.  Might be slow, so try `eigenproblem` first.
        zero_index = find(abs(diag(D) - 0) < 1E-9);       
        
        if(isempty(zero_index))
            if(settings.display)
                disp('Cannot find eigenvalue of 0.');
            end
            success = false;
            f = NaN;
            return;
        end
        f = V(:,zero_index) / sum(V(:,zero_index)); %normalize to sum to 1.  Could add other normalizations using the grid 'x' depending on settings.normalization 
        success = true;
        
     elseif(strcmp(settings.method, 'LLS')) %Solves a linear least squares problem adding in the sum constraint
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
             settings.preconditioner = 'incomplete_LU'; %Default is incomplete_LU
         end
         
         if(strcmp(settings.preconditioner,'jacobi'))
             preconditioner = diag(diag(A)); %Jacobi preconditioner is easy to calculate.  Helps a little
         elseif(strcmp(settings.preconditioner,'incomplete_cholesky'))
             %Matter if it is negative or positive?  Possible this is doing it incorrectly.
             preconditioner =-ichol(-A, struct('type','ict','droptol',1e-3,'diagcomp',1));% ichol(A, struct('diagcomp', 10, 'type','nofill','droptol',1e-1)); %matlab formula exists
         elseif(strcmp(settings.preconditioner,'incomplete_LU'))
             %Matter if it is negative or positive?
             [L,U] = ilu(A);
              preconditioner = L;
         elseif(strcmp(settings.preconditioner,'none'))
             preconditioner = [];
         else
             assert(false, 'unsupported preconditioner');
         end
        
         if(isfield(settings, 'initial_guess'))
             initial_guess = settings.initial_guess / sum(settings.initial_guess); %It normalized to 1 for simplicity.
         else
             initial_guess = [];
         end
        
        Delta = x(2) - x(1);         
        [f,flag,relres,iter] = lsqr([A';ones(1,I)], sparse([zeros(I,1);1]), tolerance, max_iterations, preconditioner, [], initial_guess); %Linear least squares.  Note tolerance changes with I
        if(flag==0)
            success = true;
        else
            if(settings.display)
                disp('Failure to converge: flag and residual');
                [flag, relres]
            end
            success = false;
            f = NaN;
        end
    end
end	
