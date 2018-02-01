function [A,F]=optimF_ConjGrad(Y,A,F,W2,D2,beta,Q,S)
% Optimize according to F

% @input :
%   Y : input Data
%   A : Matrix of factor images (vectorized)
%   F : Matrix of factor curves
%   W2 : Diagonal matrix with weight for each pixel
%   D2 : Square of the  distance to masks matrix
%   beta : Weight of the fRegF criterion
%   Q : Vector with weight for each image
%   S : Vector with weight for Tychonof regularization

% @return :
%   A : Matrix of updated factor images (vectorized)
%   F : Matrix of updates factor curves

maxIterS=5; % Maximum number of iteration of gradient descent
methodCoeff=1; % Method to obtain conjugate gradient

nbPix=size(Y,1);
nbIm=size(Y,2);
K=size(A,2);

for iterS=1:maxIterS  
    %% Compute gradient
    dfObj=grad_fObj(Y,A,F,W2,D2,1,0,beta,0,Q,S);
    
    %% Compute conjugate gradient
    if iterS==1
        conjGrad=dfObj;
    else
        conjGrad=conjugateGrad(dfObj,dfObjPrev,conjGradPrev,methodCoeff);
    end

    %% Descent of a "step" in the opposite direction of the conjugate gradient
    F = descentF(F,conjGrad(K*nbPix+1:K*nbPix+K*nbIm),A,W2,beta,Q,S);
    
    %% Correction and normalization of A and F
    [A,F] = CorrectAndNormalize(A,F);
    
    conjGradPrev=conjGrad;
    dfObjPrev=dfObj;
    
end

end

