%Takes the discretized operator A, the grid x, and finds the stationary distribution f.
function [v, success] = simple_HJBE_discretized_univariate(A, x, u, rho, settings)
   I = length(x);
   assert(I == size(A,1) && I == size(A,2)); %Make sure sizes match
   Delta = x(2) - x(1);
   
   if nargin < 5
       settings.default = true; %Just creates as required.
   end
   if(~isfield(settings, 'method'))
        settings.method = 'sparse_system';
   end
    if ~isfield(settings, 'display')
       settings.display = false; %Tolerance
    end
   
    if(strcmp(settings.method, 'sparse_system'))
        %Solve as a simple sparse system of equations.
        %More advanced solvers could use preconditioners, etc.
        v = (Delta * rho * speye(I) - A)\(Delta * u);
        success = true;
    end
 
end	
