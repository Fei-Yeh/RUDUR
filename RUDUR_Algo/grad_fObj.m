function [dfObj] = grad_fObj(Y,A,F,W2,D2,varMin,alpha,beta,mu)
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

% @return : 
%   dfObj : (Vectorized) Gradient of the objective function (ac. to A or F)


nbPix=size(Y,1);
nbIm=size(Y,2);

if varMin==0
    g_fWLS=grad_fWLS(Y,A,F,W2,0);
    if alpha==0
        g_fROI=0;
    else
        g_fROI=grad_fROI(A,D2,nbIm,mu);
    end 
    dfObj=g_fWLS+alpha*g_fROI;
else
    g_fWLS=grad_fWLS(Y,A,F,W2,1);
    if beta==0
        g_fTik=0;
    else
        g_fTik=grad_fTik(F,nbPix);
    end
    dfObj=g_fWLS+beta*g_fTik;
end

end

