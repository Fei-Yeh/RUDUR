function [fObj,fQWLS,fROI,fTik] = compute_fObj(Y,A,F,W2,D2,alpha,beta,mu,Q,S)
% Compute the value of the objective function for given factors (F) and
% factor images (A)

% @input :
%   Y : input Data
%   A : Matrix of factor images (vectorized)
%   F : Matrix of factors
%   W2 : Diagonal matrix with weight for each pixel
%   D2 : Square of the  distance to masks matrix
%   alpha : Weight of the fPrior criterion
%   beta : Weight of the fRegF criterion
%   mu : Parameter of the smoothed-l1 norm (sigma==0 => l1 norm)
%   Q : Vector with weight for each image
%   S : Vector with weight for Tychonof regularization

% @return : 
%   fObj : Value of the objective function
%   fQWLS : Value of fWLS criterion (relaxed data fidelity)
%   fROI : Value of fROI criterion
%   fTik : Value of fTik criterion


fQWLS=compute_fQWLS(Y,A,F,W2,Q);

if alpha==0
    fROI=0;
else
    fROI=compute_fROI(A,D2,mu);
end

if beta==0
    fTik=0;
else
    fTik=compute_fTik(F,S);
end

fObj=fQWLS+alpha*fROI+beta*fTik;

end

