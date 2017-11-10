%% Modified from http://www.princeton.edu/~moll/HACTproject/HJB_NGM_implicit.m by Benjamin Moll
clc;
close all;
clearvars -except  MADMaxDenseN MADMaxSparseFracForFull MADMinSparseN;
format compact

tic;

%parameters for neoclassical growth model 
gamma = 2;
alpha = 0.3;
delta = 0.05;
rho = 0.05;
Aprod = 1;

%calculate steady state level of capital 
kss = (alpha*Aprod/(rho+delta))^(1/(1-alpha)); %(9) 

%set up uniform grid with I discrete points
I=10000;
kmin = 0.001*kss;
kmax = 2*kss;
k = linspace(kmin,kmax,I)';
dk = (kmax-kmin)/(I-1);

crit = 10^(-6);
Delta = 1000;

dVf = zeros(I,1);
dVb = zeros(I,1);
c = zeros(I,1);

%INITIAL GUESS
v0 = (Aprod.*k.^alpha).^(1-gamma)/(1-gamma)/rho;
v = v0;

maxit=10;
for n=1:maxit
    V = v;
    % forward difference
    dVf(1:I-1) = (V(2:I)-V(1:I-1))/dk; %(11)
    dVf(I) = (Aprod.*kmax.^alpha - delta.*kmax)^(-gamma); %(13) state constraint, for stability
    % backward difference
    dVb(2:I) = (V(2:I)-V(1:I-1))/dk; %(12)
    dVb(1) = (Aprod.*kmin.^alpha - delta.*kmin)^(-gamma); %(14) state constraint, for stability
        
    %consumption and savings with forward difference
    cf = dVf.^(-1/gamma);
    muf = Aprod.*k.^alpha - delta.*k - cf; %(17)
    %consumption and savings with backward difference
    cb = dVb.^(-1/gamma);
    mub = Aprod.*k.^alpha - delta.*k - cb; %(18)
    %consumption and derivative of value function at steady state
    c0 = Aprod.*k.^alpha - delta.*k;
    dV0 = c0.^(-gamma);
    
    % dV_upwind makes a choice of forward or backward differences based on
    % the sign of the drift    
    If = muf > 0; %below steady state
    Ib = mub < 0; %above steady state
    I0 = (1-If-Ib); %at steady state

    dV_Upwind = dVf.*If + dVb.*Ib + dV0.*I0; %(19) important to include third term
    c = dV_Upwind.^(-1/gamma);
    u = c.^(1-gamma)/(1-gamma); %(24)
    
    %CONSTRUCT MATRIX
    X = -min(mub,0)/dk;
    Y = -max(muf,0)/dk + min(mub,0)/dk;
    Z = max(muf,0)/dk;
    
    %full matrix: slower
    %     for i=2:I-1
    %         A(i,i-1) = x(i);
    %         A(i,i) = y(i);
    %         A(i,i+1) = z(i);
    %     end
    %     A(1,1)=y(1); A(1,2) = z(1);
    %     A(I,I)=y(I); A(I,I-1) = x(I);
   
    %sparse matrix: faster
    A =spdiags(Y,0,I,I)+spdiags(X(2:I),-1,I,I)+spdiags([0;Z(1:I-1)],1,I,I);
    B = (rho + 1/Delta)*speye(I) - A;
    
    b = u + V/Delta;
    V = B\b; %SOLVE SYSTEM OF EQUATIONS
    Vchange = V - v;
    v = V;   

    dist(n) = max(abs(Vchange));
    if dist(n)<crit
        disp('Value Function Converged, Iteration = ')
        disp(n)
        break
    end
end
toc;

% Graphs
set(gca,'FontSize',14)
plot(dist,'LineWidth',2)
grid
xlabel('Iteration')
ylabel('||V^{n+1} - V^n||')

kdot = Aprod.*k.^alpha - delta.*k - c;
Verr = c.^(1-gamma)/(1-gamma) + dV_Upwind.*kdot - rho.*V;

set(gca,'FontSize',14)
plot(k,Verr,'LineWidth',2)
grid
xlabel('k')
ylabel('Error in HJB Equation')
xlim([kmin kmax])

set(gca,'FontSize',12)
plot(k,V,'LineWidth',2)
grid
xlabel('k')
ylabel('V(k)')
xlim([kmin kmax])

set(gca,'FontSize',14)
plot(k,c,'LineWidth',2)
grid
xlabel('k')
ylabel('c(k)')
xlim([kmin kmax])

set(gca,'FontSize',14)
plot(k,kdot,k,zeros(1,I),'--','LineWidth',2)
grid
xlabel('$k$','FontSize',16,'interpreter','latex')
ylabel('$s(k)$','FontSize',16,'interpreter','latex')
xlim([kmin kmax])
%print -depsc HJB_NGM.eps