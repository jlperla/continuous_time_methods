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
Aprod = 1; %factor of production A 

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

dV_F = zeros(I,1);
dV_B = zeros(I,1);
c = zeros(I,1);

%INITIAL GUESS
v0 = (Aprod.*k.^alpha).^(1-gamma)/(1-gamma)/rho;
v = v0;

maxit=10;
for n=1:maxit
    V = v;
    % forward difference
    dV_F(1:I-1) = (V(2:I)-V(1:I-1))/dk; %(11)
    dV_F(I) = (Aprod.*kmax.^alpha - delta.*kmax)^(-gamma); %(13) state constraint, for stability
    % backward difference
    dV_B(2:I) = (V(2:I)-V(1:I-1))/dk; %(12)
    dV_B(1) = (Aprod.*kmin.^alpha - delta.*kmin)^(-gamma); %(14) state constraint, for stability
        
    %consumption and savings with forward difference
    c_F = dV_F.^(-1/gamma);
    mu_F = Aprod.*k.^alpha - delta.*k - c_F; %(17)
    %consumption and savings with backward difference
    c_B = dV_B.^(-1/gamma);
    mu_B = Aprod.*k.^alpha - delta.*k - c_B; %(18)
    %consumption and derivative of value function at steady state
    c_0 = Aprod.*k.^alpha - delta.*k;
    dV_0 = c_0.^(-gamma);
    
    %average of forward and backward differences
    dV_avg = (dV_B + dV_F)/2;
    
    % dV_upwind makes a choice of forward or backward differences based on
    % the sign of the drift
    
    % original code 
    % strict or weak inequality makes no difference in this case
    
    If = mu_F > 0; %below steady state
    Ib = mu_B < 0; %above steady state
    I0 = (1-If-Ib); %at steady state
    dV_Upwind = dV_F.*If + dV_B.*Ib + dV_0.*I0; %(19) important to include third term
    %dV_Upwind = dV_F.*If + dV_B.*Ib + dV_avg.*I0; %substituting
                                                    %average
                                                            
    %variation 1: 
    %If = mu_F >= 0;
    %Ib = (mu_B <= 0) & (mu_F < 0); 
    %I0 = (mu_F < 0) & (mu_B > 0);
    %dV_Upwind = dV_F.*If + dV_B.*Ib + dV_0.*I0;
    
    %variation 2: using only mu_F
    %If = mu_F >=0;
    %Ib = mu_F < 0;
    %dV_Upwind = dV_F.*If + dV_B.*Ib;
    
    %variation 3: taking the average  
    %If = mu_F >= 0; 
    %Ib = (mu_B <= 0) & (mu_F < 0); 
    %I0 = (mu_F < 0) & (mu_B > 0); 
    %dV_Upwind = dV_F.*If + dV_B.*Ib + dV_avg.*I0;
            
    c = dV_Upwind.^(-1/gamma);
    u = c.^(1-gamma)/(1-gamma); %(24)
    
    %CONSTRUCT MATRIX
    mu_B_m = min(mu_B,0); %General notation of plus/minus.
	mu_F_p = max(mu_F,0); 
    X = -mu_B_m/dk;    
    Y = -mu_F_p/dk + mu_B_m/dk;
    Z = mu_F_p/dk;
    
    %full matrix: slower
    %     for i=2:I-1
    %         A(i,i-1) = x(i);
    %         A(i,i) = y(i);
    %         A(i,i+1) = z(i);
    %     end
    %     A(1,1)=y(1); A(1,2) = z(1);
    %     A(I,I)=y(I); A(I,I-1) = x(I);
   
    % Construct sparse A matrix  
    A =spdiags(Y,0,I,I)+spdiags(X(2:I),-1,I,I)+spdiags([0;Z(1:I-1)],1,I,I);
    B = (rho + 1/Delta)*speye(I) - A; %(25)
    
    b = u + V/Delta; %(26) 
    V = B\b; %(27) SOLVE SYSTEM OF EQUATIONS
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