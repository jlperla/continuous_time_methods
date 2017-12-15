% this is a function take input of Delta_p and Delta_m,v and create
% v1-omega*v
function residual = ValueMatch(v,z,Delta_p,Delta_m,h_p,h_m)
    I = length(Delta_p);
    N = length(v)/I; % v is N*I
    alpha = 2.1;
    F_p = @(z) alpha*exp(-alpha*z);
    eta = 1;
   
    %Trapezoidal weights, adjusted for non-uniform trapezoidal weighting
    omega_bar = (Delta_p + Delta_m)/2; %(52) though the corners are wrong
    omega_bar(1) = Delta_p(1)/2; %(51)
    omega_bar(end) = Delta_m(end)/2; %(51)  
    omega = omega_bar .* F_p(z);  %(19) Adjusted for the PDF to make (20) easy.
    
    %Stacking for the Omega
    %This is the eta = 0 case
    %omega_tilde = ([1;zeros(I-1,1)] - omega)';
    %Omega = kron(speye(N),omega_tilde); %omega_tsilde as block diagonal 
    
    %Alternative Omega with the alpha weighting
    Omega = sparse(N, N*I); %prealloate
    for n=1:N-1
        Omega(n, (n-1)*I+1:(n+1)*I) = [([1;zeros(I-1,1)] - (1-eta) * omega)'   -eta * omega'];
    end
    Omega(N, end - I + 1:end) = ([1;zeros(I-1,1)] - omega)'; %Adds in corner without the eta weighting.
    
    residual = Omega * v;
return
    
    