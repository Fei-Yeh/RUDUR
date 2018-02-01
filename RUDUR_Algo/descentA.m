function [newA] = descentA(A,direction,F,W2,Q)
% Update A, in the opposite of "direction", with a "step"

% @input :
%   A : Matrix of factor images (vectorized)
%   direction : Opposite direction of the update (vectorized)
%   F : Matrix of factor curves
%   W2 : Diagonal matrix with weight for each pixel
%   Q : Vector with weight for each image

% @return : 
%   newA : Updated matrix of factor images

nbPix=size(A,1);
K=size(A,2);

nFFt=norm(F*diag(Q.*Q)*F','fro');
nu=2*repmat(W2',1,K)*nFFt; % Lipschitz constant
dir=reshape(direction,nbPix,K);
newA=A-0.9*dir./nu;
end

