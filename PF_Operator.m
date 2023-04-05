function P = PF_Operator(PsiX, PsiY, Lambda, gamma  )
% This function compute the approximaiton of Perron-Frobenuis Operator
% PsiX:= Lited input data
% PsiY:= Lifted output data
% Lambda:= int(Psi(x)*Psi(x))dx
% gamma:= Factor of regularization
%
yalmip('clear');
if  nargin == 2
    Lambda = eye(size(PsiX,1));
    gamma = 0;
elseif nargin == 3
    gamma = 0;
elseif nargin == 4
else
    disp('Error')
end
%
%NSDMD
M = size(PsiX,2); K = size(PsiX,1);
G = 0;  A = 0;
%
for i = 1:M
    G = G + PsiX(:,i)*PsiX(:,i).';
    A = A + PsiX(:,i)*PsiY(:,i).';
end
G = G/M;
A = A/M;
% \hat{G} and \hat{G}

G1 = G/Lambda;
A1 = A/Lambda;
%
%*********  Constriant Least-Square Problem, YALMIP. *****
Pt = sdpvar(K, K,'full');
%Objective =  norm(G1*Pt-A1,'fro'); % trace(Q); %
Objective =  norm(G1*Pt-A1, 'fro')+ gamma*(norm(Pt, 'fro'));
Constraints = [];
%Constraints = [Constraints, Pt >= 0];
Constraints = [Constraints, Pt*ones(K,1) ==  ones(K,1)];
opt = sdpsettings('solver','gurobi','verbose',0,'cachesolvers',1);
optimize(Constraints, Objective, opt)
%*********
PFT= value(Pt);
P = PFT.';
end