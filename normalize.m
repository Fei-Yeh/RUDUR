function [A_N,F_N] = normalize(A,F)
% Normalize A and counter normalize F

[nbPix,K]=size(A);
nbIm=size(F,2);
F_N=zeros(K,nbIm);
A_N=zeros(nbPix,K);
for k=1:K
    coeff=nbPix/sum(abs(A(:,k)))/K;
    F_N(k,:)=F(k,:)/coeff;
    A_N(:,k)=A(:,k)*coeff;
end

end

