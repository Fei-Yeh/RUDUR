function [dfQWLS] = grad_fQWLS(Y,A,F,W2,varMin,Q)
% Compute the gradient of fQWLS criterion
% according to A if  varMin==0 or else to F

% @input :
%   Y : input Data
%   A : Matrix of factor images (vectorized)
%   F : Matrix of factors
%   W2 : Vector containing weight for each pixel
%   varMin :    0-> Minimize according to A
%               1-> Minimize according to F
%   Q : Vector with weight for each image

% @return : 
%   dfQWLS : (Vectorized) Gradient of fWLS according to A (or F)

nbPix=size(Y,1);
nbIm=size(Y,2);
K=size(A,2);

dfQWLS=zeros(K*nbPix+K*nbIm,1); % Vectorized

%% Compute weighted error for each pixel at each time
Y0=A*F;

errorY=(Y0-Y);
for t=1:nbIm
    errorY(:,t)=Q(t)*errorY(:,t);
end

if varMin==0 %% Gradient selon A
    gradA=2*errorY*F';
    gradA=gradA.*(repmat(W2',1,K));
    dfQWLS(1:K*nbPix)=gradA(:);
else %% Gradient selon F
    gradF=2*errorY'*(A.*repmat(W2',1,K));
    dfQWLS(K*nbPix+1:K*nbPix+K*nbIm)=gradF(:);
end

end