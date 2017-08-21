function [A,F]=optimA_ConjGrad(Y,A,F,W2,D2,alpha,mu)
% Optimize according to A

% @input :
%   Y : input Data
%   A : Matrix of factor images (vectorized)
%   F : Matrix of factors
%   W2 : Diagonal matrix with weight for each pixel
%   D2 : Square of the  distance to masks matrix
%   alpha : Weight of the fPrior criterion
%   mu : Parameter of the smoothed-l1 norm (sigma==0 => l1 norm)

% @return :
%   A : Matrix of updated factor images (vectorized)
%   F : Matrix of updates factor curves


maxIterS=5; % Maximum number of iteration of gradient descent
methodCoeff=1; % Method to obtain conjugate gradient

nbPix=size(Y,1);
K=size(A,2);


for iterS=1:maxIterS  
    %% Compute gradient
    dfObj=grad_fObj(Y,A,F,W2,D2,0,alpha,0,mu);
    
    %% Compute conjugate gradient
    if iterS==1
        conjGrad=dfObj;
    else
        conjGrad=conjugateGrad(dfObj,dfObjPrev,conjGradPrev,methodCoeff);
    end

    %% Descent of a "step" in the opposite direction of the conjugate gradient
    A = descentA(A,conjGrad(1:K*nbPix),F,W2);
    
    %% Correction and normalization of A and F
    [A,F] = CorrectAndNormalize(A,F);
    conjGradPrev=conjGrad;
    dfObjPrev=dfObj;
    
end

end

