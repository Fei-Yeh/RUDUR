function [Aest,Fest,fObj,fWLS,fPrior,fRegF,finalIter] = rudur(Y,M,alpha,beta,gamma,mu)
% Main Function of RUDUR algorithm

% @input :
%   Y : input Data
%   M : Matrix of masks (ROIs)
%   alpha : Weight of the fPrior criterion
%   beta : Weight of the fRegF criterion
%   gamma : Weighting data fidelity constant
%   mu : Parameter of the smoothed-l1 norm (sigma==0 => l1 norm)

% @return : 
%   Aest : Estimation of source images A
%   Fest : Estimation of TACs F
%   fObj : Value of the objective function fRUDUR
%   fWLS : Value of fWLS criterion (relaxed data fidelity)
%   fROI : Value of fROI criterion
%   fTik : Value of fTik criterion

nbIm=size(Y,2);
[nbRow,nbCol,K]=size(M);

M_1D=reshape(M,nbRow*nbCol,K);

%% Stopping criterion
stopIter=0;
maxIter=1000; % Maximum number of iterations
finalIter=maxIter; % Iteration reached for convergence
epsilon=10^-7; % Precision of the relative difference between two iterations to obtain to satisfy stopping criterion
stopIterMax=5; % Number of iteration satisfying stopping criterion necessary to converge

%% Records of Criterion values
fObj=zeros(1,maxIter);
fWLS=zeros(1,maxIter);
fPrior=zeros(1,maxIter);
fRegF=zeros(1,maxIter);

%% Build matrix of distance D with the masks M
[D,~]=buildPrior(M);
D2=D.*D;

%% Build diagonal matrix of weight W
W=buildWeight(D,gamma);
W2=W.*W;

%% Initialize algorithm        
Ainit=M_1D;
Finit=zeros(K,nbIm);
for k=1:K
    Finit(k,:)=mean(Y(M_1D(:,k)==1,:));
end

Yinit=Ainit*Finit;
for t=1:nbIm
    if norm(Yinit(:,t))>0
        Finit(:,t)=Finit(:,t)*norm(Y(:,t))/norm(Yinit(:,t));
    end
end

[Aest,Fest]=CorrectAndNormalize(Ainit,Finit);

varMin=1; % Start minimizing according to F
iter=1;

% Save criterion values with initial values
[fObj(iter),fWLS(iter),fPrior(iter),fRegF(iter)] = compute_fObj(Y,Aest,Fest,W2,D2,alpha,beta,mu);


%% Start optimization
for iter=2:maxIter
    if varMin==0    % Fix F and minimize according to A
        [Aest,Fest]=optimA_ConjGrad(Y,Aest,Fest,W2,D2,alpha,mu);
    else            % Fix A and minimize according to F
        [Aest,Fest]=optimF_ConjGrad(Y,Aest,Fest,W2,D2,beta);
    end
    varMin=1-varMin; % Minimize according to the other matrix on next iteration
    %% Record criterion values
    [fObj(iter),fWLS(iter),fPrior(iter),fRegF(iter)] = compute_fObj(Y,Aest,Fest,W2,D2,alpha,beta,mu);
    
    %% Stopping criterion
    if (iter>2 && (fObj(iter-2)-fObj(iter))/fObj(iter-2)<epsilon)
        stopIter=stopIter+1;
    else
        stopIter=max(0,stopIter-2);
    end
    if stopIter==stopIterMax
        finalIter=iter;
        break;
    end
    
end

end

