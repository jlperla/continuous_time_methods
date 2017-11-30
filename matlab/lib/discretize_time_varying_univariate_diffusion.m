%For algebra and equation numbers, see the 'operator_discretization_finite_differences.pdf'

%This function takes a grid on [x_min, x_max], [t_min, t_max] and discretizing a general diffusion defined by the following SDE
%d x_t = mu(t, x_t)dt + sigma(t, x_t)^2 dW_t
%Subject to reflecting barrier at x_min and x_max and a stationary requirement at t_max

%Pass in the vector of the grid x, and the vectors of mu and sigma_2 at the nodes, and returns a sparse discretized operator.
%The mu and sigma_2 are assumed to be already stacked correctly (i.e., keeping all time together).


%This so far is the explicit time procedure (Nov.26)
function [A, Delta_p, Delta_m, h_p, h_m] = discretize_time_varying_univariate_diffusion(t, x, mu, sigma_2, check_absorbing_states)
    if nargin < 4
        check_absorbing_states = true;
    end
	I = length(x); %number of grid variables for x
	N = length(t);

	%Could check if the grid is uniform
	Delta_p = [diff(x)' (x(I)-x(I-1))]'; %(35) Find distances between grid points.
    Delta_m = [x(2)-x(1) diff(x)']'; % %(34) \Delta_{i, -}
	h_p = [diff(t)' (t(N)-t(N-1))]'; % (67) h_{+}
    h_m = [t(2)-t(1) diff(t)']'; % %(68) h{i, -}
    
    % stack delta's into R^NI
    Delta_stack_p=repmat(Delta_p,N,1);
    Delta_stack_m=repmat(Delta_m,N,1); 
    D_h_stack_p = kron(1./h_p(1:N), ones(I,1)); %Stacks up 1/h_+ for each spatial dimension.
    
   %% Construct sparse A matrix with non-uniform grid (uniform case is just a generalization of non-uniform)
    %For non-uniform grid, \Delta_{i, +}=x_{i+1} - x_{i} and \Delta_{i, -}=x_{i} - x_{i-1}
    mu_m = min(mu,0); %General notation of plus/minus.
    mu_p = max(mu,0); 		
    X = - mu_m./Delta_stack_m + sigma_2 ./(Delta_stack_m.*(Delta_stack_p + Delta_stack_m)); %(74)
    Y = - mu_p./Delta_stack_p + mu_m./Delta_stack_m - sigma_2./(Delta_stack_p .* Delta_stack_m); % (75)
    Z =  mu_p./Delta_stack_p + sigma_2 ./ (Delta_stack_p.*(Delta_stack_p + Delta_stack_m)); % (76)
    
    %Creates A^n matrices that's are n-specific(time grid specific)
    
    % Construct the A matrix in pieces
    A = spdiags([nan(I,1); D_h_stack_p], [I], N*I, N*I); %Start with the off-diagonal of 1/h_p
    for n=1:N
        indx = I*(n-1)+1;
        Xn = X(indx:indx+I-1);
        Yn = Y(indx:indx+I-1);
        Zn = Z(indx:indx+I-1);
        A_n = spdiags([[Xn(2:I); NaN] Yn [NaN; Zn(1:I - 1)]], [-1 0 1], I,I);% (77) for each time node indexed by n, A_n is different as mu_n changes. The procedure similar to that in time-invariant case   
        A_n(1,1) = Yn(1) + Xn(1);%Reflecting barrier, top corner of (77)
        A_n(I,I) = Yn(I) + Zn(I);%Reflecting barrier, bottom corner of (77)
        A(indx:indx+I-1,indx:indx+I-1)=A_n;
    end
    A = A + spdiags(-[D_h_stack_p(1:I*(N-1)); zeros(I,1)], [0], N*I, N*I); %Putting 0's at the end

end	