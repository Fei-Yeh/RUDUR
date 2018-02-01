function [dfROI] = grad_fROI(A,D2,nbIm,mu)
% Compute the gradient of fPrior criterion

% @input :
%   A : Matrix of factor images (vectorized)
%   D2 : Square of the  distance to masks matrix
%   nbIm : Number of images in the data sequence
%   mu : Parameter of the smoothed-l1 norm (sigma==0 => l1 norm)

% @return : 
%   dfROI: (Vectorized) Gradient of fROI


nbPix=size(A,1); % Number of pixel
K=size(A,2); % Number of factors

dfROI=zeros(K*nbPix+K*nbIm,1); % Vectorized

if mu==0 % Compute gradient of l1-norm
    grad=D2.*((A>0)-(A<0));
else % Compute gradient of smoothed l1-norm
    D4=D2.*D2;
    grad=(D4.*A)./sqrt(A.*A.*D4+mu*mu);
end
dfROI(1:K*nbPix)=grad(:);
end

