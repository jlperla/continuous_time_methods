%For algebra and equation numbers, see the 'operator_discretization_finite_differences.pdf'

%This function takes a grid on [x_min, x_max], [t_min, t_max] and discretizing a general diffusion defined by the following SDE
%d x_t = mu(t, x_t)dt + sigma(t, x_t)^2 dW_t
%Subject to reflecting barrier at x_min and x_max and a stationary requirement at t_max

%Pass in the vector of the grid x, and the vectors of mu and sigma_2 at the nodes, and returns a sparse discretized operator.
%The mu and sigma_2 are assumed to be already stacked correctly (i.e., keeping all time together).
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
    
    %Use mu, sigma_2, etc. to calculate A    
	A = NaN; %TODO: Calculate correctly

end	