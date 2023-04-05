function K = K_Operator(PsiX, PsiY,  gamma  )
% This function compute the approximaiton of Koopman Operator
% PsiX:= Lited input data
% PsiY:= Lifted output data
% gamma:= Factor of regularization
%
%EDMD
M = size(PsiX,2); K = size(PsiX,1);
G = 0;  A = 0;
%
for i = 1:M
    G = G + PsiX(:,i)*PsiX(:,i).';
    A = A + PsiX(:,i)*PsiY(:,i).';
end
G = G/M;
A = A/M;

if  nargin == 2
    if rcond(G) > eps
        K = G\A;
    else
        K = pinv(G)*A;
    end

elseif nargin == 3
    yalmip('clear');
    %*********  Constriant Least-Square Problem, YALMIP. *****
    Kt = sdpvar(K, K,'full');
    Objective =  norm(G*Kt-A, 'fro')+ gamma*(norm(Kt, 'fro'));
    Constraints = [];
    opt = sdpsettings('solver','gurobi','verbose',0,'cachesolvers',1);
    optimize(Constraints, Objective, opt)
    %*********
    K = value(Kt);
else
    disp('Error')
end
end