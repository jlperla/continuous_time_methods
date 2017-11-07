%% Modified from http://www.princeton.edu/~moll/HACTproject/HJB_NGM_implicit.m by Benjamin Moll
clc;
close all;
clearvars -except  MADMaxDenseN MADMaxSparseFracForFull MADMinSparseN;
format compact

tic;

s = 2;
a = 0.3;
d = 0.05;
r = 0.05;
Aprod = 1;

kss = (a*Aprod/(r+d))^(1/(1-a));

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
v0 = (Aprod.*k.^a).^(1-s)/(1-s)/r;
v = v0;

maxit=10;
for n=1:maxit
    V = v;
    % forward difference
    dVf(1:I-1) = (V(2:I)-V(1:I-1))/dk;
    dVf(I) = (Aprod.*kmax.^a - d.*kmax)^(-s); %state constraint, for stability
    % backward difference
    dVb(2:I) = (V(2:I)-V(1:I-1))/dk;
    dVb(1) = (Aprod.*kmin.^a - d.*kmin)^(-s); %state constraint, for stability
        
    %consumption and savings with forward difference
    cf = dVf.^(-1/s);
    muf = Aprod.*k.^a - d.*k - cf;
    %consumption and savings with backward difference
    cb = dVb.^(-1/s);
    mub = Aprod.*k.^a - d.*k - cb;
    %consumption and derivative of value function at steady state
    c0 = Aprod.*k.^a - d.*k;
    dV0 = c0.^(-s);
    
    % dV_upwind makes a choice of forward or backward differences based on
    % the sign of the drift    
    If = muf > 0; %below steady state
    Ib = mub < 0; %above steady state
    I0 = (1-If-Ib); %at steady state

    dV_Upwind = dVf.*If + dVb.*Ib + dV0.*I0; %(19) important to include third term
    c = dV_Upwind.^(-1/s);
    u = c.^(1-s)/(1-s);
    
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
    B = (r + 1/Delta)*speye(I) - A;
    
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

kdot = Aprod.*k.^a - d.*k - c;
Verr = c.^(1-s)/(1-s) + dV_Upwind.*kdot - r.*V;

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