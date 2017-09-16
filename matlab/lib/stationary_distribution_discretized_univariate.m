%Takes the discretized operator A, the grid x, and finds the stationary distribution f.
function [f, success] = stationary_distribution_discretized_univariate(A, x, settings)
   I = length(x);
   if nargin < 3
       settings.default = true; %Just creates as required.
   end
   if(~isfield(settings, 'method'))
        settings.method = 'eigenproblem';
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
        if(isfield(settings, 'num_basis_vectors')) %OTherwise use the default
           opts.p = settings.num_basis_vectors; %Number of Lanczos basis vectors.
        end
         if(isfield(settings, 'max_iterations')) %OTherwise use the default
          opts.maxit = settings.max_iterations; %Number of iterations
        end
        [V,D] = eigs(A' + speye(I),1,'lr', opts);%The eigenvalue with the largest real part should be the unity one
        if(abs(D - 1.0) > 1E-9) %The 'lr' one is hopefully unity, but maybe not.
            if(settings.display)
                disp('The eigenvalue is not unity.  Try increasing the num_basis_vectors or max_iterations.  Then eigenproblem_all');
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
        [V,D] = eigs(A' + speye(I), I, opts); %Gets all of the eigenvalues and eigenvectors.  Might be slow, so try `eigenproblem` first.
        
        f = V / sum(V); %normalize to sum to 1.  Could add other normalizations using the grid 'x' depending on settings.normalization 
        success = true;      
   elseif(strcmp(settings.method, 'LLS')) %Solves a linear least squares problem adding in the sum constraint
         if(isfield(settings, 'max_iterations')) %OTherwise use the default
            max_iterations = settings.max_iterations; %Number of iterations
         else
            max_iterations = 50; %The default is too small for our needs
         end
         if(isfield(settings, 'tolerance')) %OTherwise use the default
            tolerance = settings.tolerance; %Number of iterations
         else
            tolerance = []; %Empty tells it to use default
         end         
        [f,flag,relres,iter] = lsqr([A';ones(1,I)], sparse([zeros(I,1);1]), tolerance, max_iterations); %Linear least squares.  Note tolerance changes with I
        if(flag==0)
            success = true;
        else
            if(setttings.display)
                disp('Failure to converge: flag and residual');
                [flag, relres]
            end
            success = false;
            f = NaN;
        end
    end
end	
