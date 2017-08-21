function [fObj] = compute_fObj_bis(step,direction,Y,A,F,W2,D2,alpha,beta,mu)
% Do the same than compute_fObj function, but only return fObj value for
% A=A-stepA*directionA
% and F=F-stepF*directionF
    
% @input :
%   step : Step used for the update
%   direction : Opposite direction of the update (vectorized)
%   Y : input Data
%   A : Matrix of factor images (vectorized)
%   F : Matrix of factors
%   W2 : Diagonal matrix with weight for each pixel
%   D2 : Square of the  distance to masks matrix
%   alpha : Weight of the fPrior criterion
%   beta : Weight of the fRegF criterion
%   mu : Parameter of the smoothed-l1 norm (sigma==0 => l1 norm)


% @return : 
%   fObj : Value of the objective function


    nbPix=size(Y,1);
    nbIm=size(Y,2);
    K=size(A,2);
    
    A=descentA(A,direction(1:K*nbPix),step,0,F,W2);
    F=descentF(F,direction(K*nbPix+1:K*nbPix+K*nbIm),step,A,W2);
    [fObj,~,~,~] = compute_fObj(Y,A,F,W2,D2,alpha,beta,mu);
end

