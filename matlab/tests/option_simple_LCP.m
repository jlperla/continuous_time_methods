%Minor modifications to http://www.princeton.edu/~moll/HACTproject/option_simple_LCP.m
%clearvars -except  MADMaxDenseN MADMaxSparseFracForFull MADMinSparseN;
%format compact;clc; close all;
addpath('../lib/');
rho = 0.05;

I= 1000;
xmin = 0; xmax=1;
x = linspace(xmin,xmax,I)';
dx = x(2)-x(1);
dx2 = dx^2;

u = x.^(0.5);
mu_bar = -0.01;
mu = ones(I,1).*mu_bar;

sig_bar = 0.01;
sig = x.*sig_bar;
sig2 = sig.^2;

S = 10*ones(I,1);
%S = 10*exp(0.05*x);


%CONSTRUCT MATRIX
X = - min(mu,0)/dx + sig2/(2*dx2);
Y = - max(mu,0)/dx + min(mu,0)/dx - sig2/dx2;
Z =  max(mu,0)/dx + sig2/(2*dx2);

A =spdiags(Y,0,I,I)+spdiags(X(2:I),-1,I,I)+spdiags([0;Z(1:I-1)],1,I,I);
A(I,I)= Y(I) + sig2(I)/(2*dx2); A(I,I-1) = X(I);

%Solve option problem as LCP
B = rho*speye(I) - A;
q = -u + B*S; 

%TRY DIFFERENT ALGORITHMS FOR SOLVING LCP
tic
%disp('Solving LCP')

%Test 1: Yuval Tassa's Newton-based LCP solver, download from http://www.mathworks.com/matlabcentral/fileexchange/20952
z0 = zeros(I,1); l = zeros(I,1); u = Inf*ones(I,1);
z = LCP(B,q,l,u,z0,false);

%Test 2: Andreas Almqvist pivoting (Lemke) LCP solver, download from http://www.mathworks.com/matlabcentral/fileexchange/41485
%[w,z,retcode] = LCPSolve(B,q);

LCP_error = max(abs(z.*(B*z + q)));
if LCP_error > 10^(-6)
    disp('LCP not solved')
    exit;
end
    
V_LCP = z+S; %calculate value function
toc;

% toc
% 
% error = z.*(B*z + q);
% plot(x,error)

v = V_LCP;
test_output.x = x;
test_output.v = v;
test_output.S = S;
test_output.error = LCP_error;
save('output_option_simple_LCP_HACT_raw.mat', 'test_output');

% plot(x,V_LCP,x,S,'--','LineWidth',2)
% set(gca,'FontSize',16)
% legend('v(x)','S(x)','Location','NorthWest')
% xlabel('x')
% %print -depsc option_simple.eps
