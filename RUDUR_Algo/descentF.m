function [newF] = descentF(F,direction,A,W2)
% Update F, in the opposite of "direction", with a "step"

% @input :
%   F : Matrix of factors curves
%   direction : Opposite direction of the update (vectorized)
%   A : Matrix of factor images (vectorized)
%   W2 : Diagonal matrix with weight for each pixel

% @return : 
%   newF : Updated matrix of factors

nbIm=size(F,2);
K=size(F,1);


mu=norm(2*(A'.*repmat(W2,K,1)*A));
newF=F'-0.9/mu*reshape(direction,nbIm,K);
newF=newF';

end