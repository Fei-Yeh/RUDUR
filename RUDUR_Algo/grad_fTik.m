function [dfTik] = grad_fTik(F,nbPix)
% Compute the gradient of fTik criterion

% @input :
%   F : Matrix of factors
%   nbPix : Number of Pixels in each image

% @return : 
%   dfTik : (Vectorized) Gradient of fRegF


nbIm=size(F,2); % Number of images
K=size(F,1); % Number of factors

dfTik=zeros(K*nbPix+K*nbIm,1); % Vectorized gradient
dfF=zeros(K,nbIm); % gradient according to F

diffF=F(:,2:end)-F(:,1:end-1); % difference between two consecutive images

dfF(:,1)=-2*diffF(:,1);
dfF(:,nbIm)=2*diffF(:,nbIm-1);
dfF(:,2:nbIm-1)=2*(diffF(:,1:nbIm-2)-diffF(:,2:nbIm-1));


for k=1:K
    dfTik(K*nbPix+(k-1)*nbIm+1:K*nbPix+k*nbIm)=dfF(k,:);
end

end

