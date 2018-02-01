function [newF] = descentF(F,direction,A,W2,beta,Q,S)
% Update F, in the opposite of "direction", with a "step"

% @input :
%   F : Matrix of factors curves
%   direction : Opposite direction of the update (vectorized)
%   A : Matrix of factor images (vectorized)
%   W2 : Diagonal matrix with weight for each pixel
%   Q : Vector with weight for each image
%   S : Vector with weight for Tychonof regularization

% @return : 
%   newF : Updated matrix of factors

nbIm=size(F,2);
K=size(F,1);
if beta~=0
    Gam=eye(nbIm,nbIm);
    Gam(1,1)=0;
    for t=2:nbIm
        Gam(t,t-1)=-1;
    end
end
if sum(Q==ones(1,nbIm))==nbIm
    nQ=1;
else
    nQ=norm(Q.*Q);
end
if beta~=0
    nu=norm(2*(A'.*repmat(W2,K,1)*A),'fro')*nQ+beta*norm(Gam'*diag(S)*Gam,'fro'); % Lipschitz constant
else
    nu=norm(2*(A'.*repmat(W2,K,1)*A),'fro')*nQ;
end
newF=F'-0.9/nu*reshape(direction,nbIm,K);
newF=newF';

end