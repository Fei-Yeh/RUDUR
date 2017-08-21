function [newA] = descentA(A,direction,F,W2)
% Update A, in the opposite of "direction", with a "step"

% @input :
%   A : Matrix of factor images (vectorized)
%   direction : Opposite direction of the update (vectorized)
%   F : Matrix of factor curves
%   W2 : Diagonal matrix with weight for each pixel

% @return : 
%   newA : Updated matrix of factor images

nbPix=size(A,1);
K=size(A,2);

nFFt=norm(F*F');
mu=2*repmat(W2',1,K)*nFFt;
dir=reshape(direction,nbPix,K);
newA=A-0.9*dir./mu;
end

