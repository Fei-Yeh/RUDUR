function [dfWLS] = grad_fWLS(Y,A,F,W2,varMin)
% Compute the gradient of fWLS criterion
% according to A if  varMin==0 or else to F

% @input :
%   Y : input Data
%   A : Matrix of factor images (vectorized)
%   F : Matrix of factors
%   W2 : Vector containing weight for each pixel
%   varMin :    0-> Minimize according to A
%               1-> Minimize according to F

% @return : 
%   dfWLS : (Vectorized) Gradient of fWLS according to A (or F)

nbPix=size(Y,1);
nbIm=size(Y,2);
K=size(A,2);

dfWLS=zeros(K*nbPix+K*nbIm,1); % Vectorized

%% Compute weighted error for each pixel at each time
Y0=A*F;

errorY=(Y0-Y);

if varMin==0 %% Gradient selon A
    gradA=2*errorY*F';
    gradA=gradA.*(repmat(W2',1,K));
    dfWLS(1:K*nbPix)=gradA(:);
else %% Gradient selon F
    gradF=2*errorY'*(A.*repmat(W2',1,K));
    dfWLS(K*nbPix+1:K*nbPix+K*nbIm)=gradF(:);
end

end