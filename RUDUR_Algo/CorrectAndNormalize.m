function [newA,newF] = CorrectAndNormalize(A,F)
% Modify negative coefficient of A and F and normalize A and Fx
% So that mean(Ak)=1/K  for each factor "k" (same for Aest)
% And Fk is counter-normalized accordingly (same for Festk)


% @input :
%   A : Matrix of factor images (vectorized) (Ground Truth)
%   F : Matrix of factors (Ground Truth)

% @return : 
%   newA : Matrix A updated
%   newF : Matrix F updated

epsilon=0.0001;

nbPix=size(A,1);
K=size(A,2);

% Zero on each negative coefficient
newA=A;
newA(newA<0)=0;
newF=F;
newF(newF<0)=0;

%% Normalization of each factor image
% and counter normalization of each factor accordingly
for k=1:K
    normk=sum(newA(:,k))/nbPix*K;
    newA(:,k)=newA(:,k)/normk;
    newF(k,:)=newF(k,:)*normk;
end
newA(newA<epsilon)=epsilon;
newF(newF<epsilon)=epsilon;
end

