%Use this to test the creation of the A matrix after it has been separated from the LCP solver as a separate function

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