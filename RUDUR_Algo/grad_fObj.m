function [dfObj] = grad_fObj(Y,A,F,W2,D2,varMin,alpha,beta,mu,Q,S)
% Compute the gradient of the objective function
% according to A if varMin==0 or else to F

% @input :
%   Y : input Data
%   A : Matrix of factor images (vectorized)
%   F : Matrix of factors
%   W2 : Diagonal matrix with weight for each pixel
%   D2 : Square of the  distance to masks matrix
%   varMin :    0-> Minimize according to A
%               1-> Minimize according to F
%   alpha : Weight of the fPrior criterion
%   beta : Weight of the fRegF criterion
%   mu : Parameter of the smoothed-l1 norm (sigma==0 => l1 norm)
%   Q : Vector with weight for each image
%   S : Vector with weight for Tychonof regularization

% @return : 
%   dfObj : (Vectorized) Gradient of the objective function (ac. to A or F)


nbPix=size(Y,1);
nbIm=size(Y,2);

if varMin==0
    g_fQWLS=grad_fQWLS(Y,A,F,W2,0,Q);
    if alpha==0
        g_fROI=0;
    else
        g_fROI=grad_fROI(A,D2,nbIm,mu);
    end 
    dfObj=g_fQWLS+alpha*g_fROI;
else
    g_fQWLS=grad_fQWLS(Y,A,F,W2,1,Q);
    if beta==0
        g_fTik=0;
    else
        g_fTik=grad_fTik(F,nbPix,S);
    end
    dfObj=g_fQWLS+beta*g_fTik;
end

end

