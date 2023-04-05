%% Data-Driven Perron-Forbenius Optimal Control for Scalar System
clear; close all; clc;
set(0,'DefaultLineLineWidth',1.8) %linewidh on plots
set(0,'defaultfigurecolor',[1 1 1])
rng(2141444)
%% Parameters Simulaiton
% Dynamics
a  =0.5; b =1;
fL_u = @(t, x, u) a*x.^3 +b*u; % dynamics
n=1; % dimension
deltaT = 0.01; % sample time
R_L = 0.1;  % local weight for input of LQR
R_NL =0.1; %globall weight forData=Driven PF Control
Nrbf =75;   % numer of bases totally
Dom = [-20 20]; % Domain
dis = 55;    %magnitud of domain
%
% Data generated from one-step and Euler discretization
nb_IC =1e4; % number the initial conditions
X0 = rand(1, nb_IC)*Dom(2)*2 - Dom(2);
[X , Y ] = Data(@(t, x, u)fL_u(t, x, 0), X0, deltaT, 0);
[X1 , Y1] = Data(@(t, x, u)fL_u(t, x, 1), X0, deltaT, 1 );
data = [X Y X1 Y1];

%% ************** Phase 1. PF approximation  ******************
%
%  ********* Basis functions **********
sig = (0.49)*(dis/(Nrbf-1)); % sigma in  function of distance
% RBF centers
cent = linspace(Dom(1), Dom(2), Nrbf);
% Remove origen
eq_0 = 0;
[~,idx] = min(sqrt(sum((cent-eq_0).^2,1))); % Find the center closest to the equilibrium point
cent(:,idx) = eq_0;
[idx] = find(sqrt(sum((cent-eq_0).^2,1))<=(1*sig)); % Find the centers within 1 sigma of origin
%liftFun = @(xx)(rbf(xx,cent,'gauss'));
Psi = @(X)(GaussRBF(X,cent, sig));

%  ********* Lifted the data  *********
PsiX =  Psi(X);
PsiY =  Psi(Y);
PsiX1 =  Psi(X1);
PsiY1 =  Psi(Y1);
%  ********* Compute the operators using NSDMD  *********
Lam = eye(Nrbf); %
%Lam = (pi*sig^2/2)^(n/2)*exp(-squareform(pdist(cent').^2)./(2*sig^2));
gamma = 0; % Factor of regularization
P_f = PF_Operator(PsiX, PsiY, Lam, gamma);
P_fg = PF_Operator(PsiX1,PsiY1, Lam, gamma);

%% ************** Phase II. Control Design  ******************
%
% ******  Cost Function Approximaiton  ******
%  First term of cost
d = zeros(Nrbf,1);
qx = @(x)x.^2; % cost function associated to the states
for i = 1:Nrbf
    dfunc =@(x) qx(x).*exp(-(1/sig)^2*sum( (x - (cent(i))).^2 ,1));
    d(i) = integral(dfunc,Dom(1),Dom(2));
end
%
%  Second term of cost
%%%%% Numerically:
D = (pi*sig^2/2)^(n/2)*exp(-squareform(pdist(cent').^2)./(2*sig^2));
% ******  Remove rows and colums the origen  ******
P_f(idx,:) = [];P_f(:,idx) = [];
P_fg(idx,:) = [];P_fg(:,idx) = [];
Lam(idx,:) = []; Lam(:,idx) = [];
D(idx,:) = []; D(:,idx) = [];
d(idx) = [];
lg_idx = length(idx);

% ****** Convex  Optimization Problem  ******
cvx_begin quiet
%cvx_solver sedumi
%cvx_solver mosek
variable v(Nrbf-lg_idx) nonnegative
variable w(Nrbf-lg_idx)
I = eye(Nrbf-lg_idx);
ctr_cost = 0;
for i = 1:(Nrbf-lg_idx)
    for j = 1:(Nrbf-lg_idx)
        if (i==j)  ~=  (i==j-1) ~=  (i-1==j) ~=  (i-2==j)  ~=  (i ==j-2) %omit for final result
        %if i==j
            ctr_cost = ctr_cost + D(i,j)*quad_over_lin(w(i),v(j));
            %toc
        end
    end
end
%
minimize(d'*v + R_NL*ctr_cost)
%minimize(0)
%
subject to
mm= 1e-3;
alpha = 1;
h = mm*ones(Nrbf-lg_idx,1);
%
(I-alpha*P_f)*v - alpha*(P_fg-P_f)*w >= h;
%
cvx_end
cvx_optval

% **** Global Controller ****
K = w./v;
K0 = zeros(Nrbf,1);
K0(setdiff(1:end,idx)) = K;
K = K0;

%**** Local controller ****
x = sym('x',[1;1]); u = sym('u',[1;1]);
Al = double(subs(jacobian(fL_u(0,x,u),x),[x;u],[0;0]));
Bl = double(subs(jacobian(fL_u(0,x,u),u),[x;u],[0;0]));
Klqr = lqr(Al, Bl, eye(size(a, 2)),R_L);
Fcl = (a - b*Klqr);
P_Lyap = lyap(Fcl, eye(n));
% Loical density function
r =200000;
rho_L=  @(x) max(((x'*P_Lyap*x)^-3  - r), 0);

%%  ************************  Close-Loop Simualation   ******************************
close all
% Euler discretization:
f_ud = @(t,x,u) (x + deltaT*fL_u(t,x,u));
% build vN and Wn acording dimensions of Lifted space
vN = zeros(Nrbf,1); wN = zeros(Nrbf,1);
vN(setdiff(1:Nrbf, idx)) =v;
wN(setdiff(1:Nrbf, idx)) =w;
%
Tmax = 300;
Nsim = Tmax/1;
Ninit = 10;
rng(2156488);   lw =1;
% Initial condition for test
x0 = rand(n,Ninit)*6-3;

% Analytical feeback control
u_ANL = @(x) -(a*x^3+x*sqrt(a^2*x^4+(1/R_NL)));
% Tests
for j = 1:Ninit
    
    x_CL = x0(:,j);
    x_open = x0(:,j);
    x_LQR = x0(:,j);
    x_NLu = x0(:,j);
    u_PF = [];
    u_LQR = [];
    u_NL = [];
    
    for i = 0:Nsim-1
        % Data-driven control
        u_PF(i+1) =(rho_L(x_CL(:,end))/(rho_L(x_CL(:,end))+Psi(x_CL(:,end))'*vN))*(-Klqr*x_CL(:,end))...
            +  ((Psi(x_CL(:,end))'*vN)/(rho_L(x_CL(:,end))+Psi(x_CL(:,end))'*vN))*(Psi(x_CL(:,end))'*K);
        % Local control
        u_LQR(i+1) = -Klqr*x_LQR(:,end);
        % Analyticcal control
        u_NL(i+1) = u_ANL(x_NLu(:,end));
        % Evaluaiton
        x_open= [x_open, f_ud(0,x_open(:,end), 0)];
        x_CL  = [x_CL,   f_ud(0,x_CL(:,end), u_PF(end))];
        %x_LQR   = [x_LQR,   f_ud(0,x_LQR(:,end), u_LQR(end))];
        x_NLu  = [x_NLu,   f_ud(0,x_NLu(:,end),   u_NL(end)) ];
    end
    u_PF(i+2) =  0;
    u_LQR(i+2) =  0;
    u_NL(i+2) =  0;
    
    figure(1)
    hold on
    plot([0:Nsim]*deltaT,x_CL(1,:),'-.r');
    %plot([0:Nsim]*deltaT,x_LQR(1,:),'-.r');
    plot([0:Nsim]*deltaT,x_NLu(1,:),'--k');
%     %
%     figure(2)
%     %JLQR = cumsum(x_LQR(1,:).^2 + ( u_LQR).^2);
%     JPF = cumsum(x_CL(1,:).^2 + (u_PF).^2);
%     J_NL = cumsum(x_NLu(1,:).^2 + (u_NL).^2);
%     plot([0:Nsim]*deltaT, JPF,'-.r','linewidth',lw); hold on
%     %plot([0:Nsim]*deltaT, JLQR,'-.r','linewidth',lw);
%     plot([0:Nsim]*deltaT, J_NL,'--k','linewidth',lw);
%     %
%     figure(3)
%     plot([0:Nsim]*deltaT, u_PF(1,:),'-.r','linewidth',1); hold on
%     plot([0:Nsim]*deltaT,u_NL(1,:) ,'--k','linewidth',1);
    
end
figure(1)
grid on;
xlabel('Time [s]','interpreter','latex');
ylabel('$x$','interpreter','latex');
set(gca,'fontsize',20)
ellip=sqrt((1/r)^(1/3)/P_Lyap);
yline(ellip,'m','linewidth',lw);
yline(-ellip,'m','linewidth',lw);
LEG = legend('$x$(P-F)','$x$(NL)');%,'control u');
set(LEG,'interpreter','latex')
title(['RBFs = ' num2str( Nrbf ) ', \sigma = ' num2str( sig ) ',   RL=' num2str(R_L) ',  RNL=' num2str(R_NL)   ])
%


%%
%x_span = Dom(1):10/Nrbf:Dom(2);
x_span = -11:0.01:11;
J_NL = 0;
U_x=[];
U_xNL=[];
Vfunc=[];
V =@(x) 0.025*(x^4 + 200*asinh(0.0707107*x^2)+ sqrt(x^4 +200*x^2));
for i=1:length(x_span)
    XX = x_span(i);
    U_x(i) = (rho_L(XX)/(rho_L(XX)+Psi(XX)'*vN))*(-Klqr*XX)...
        +  ((Psi(XX)'*vN)/(rho_L(XX)+Psi(XX)'*vN))*(Psi(XX)'*K);
    U_xNL(i) = u_ANL(XX);
    U_LQR(i) = -Klqr*XX;
    Vfunc(i) = V(XX);
end
set(0,'DefaultLineLineWidth',3.2) 
figure
plot(x_span, U_x,'k')
hold on
plot(x_span, U_xNL, '--r')
%plot(x_span, U_LQR, '--b')
spanL = -ellip:0.01:ellip;
for i=1:length(spanL)
    XX = spanL(i);
    U_Lo(i) = -Klqr*XX;
end
plot(spanL, U_Lo, 'm')
grid on;
ylabel('$u(x)$','interpreter','latex');
xlabel('$x$','interpreter','latex');
set(gca,'fontsize', 15)
LEG = legend('$u(x)$(P-F)','$u(x)$(NL)');%,'control u');
set(LEG,'interpreter','latex')
xlim([-3 3])


