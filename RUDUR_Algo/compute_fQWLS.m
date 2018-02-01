function [ fQWLS ] = compute_fQWLS(Y,A,F,W2,Q)
% Compute the value of fWLS criterion

% @input :
%   Y : input Data
%   A : Matrix of factor images (vectorized)
%   F : Matrix of factors
%   W2 : Vector containing weight for each pixel
%   Q : Vector with weight for each image

% @return : 
%   fWLS : Value of fWLS

Y0=A*F;
errorYw=(Y0-Y);
for t=1:length(Q)
    errorYw(:,t)=Q(t)*errorYw(:,t);
end
errorY2=(errorYw.*errorYw);
errorY2sum=sum(errorY2,2);
fQWLS=W2*errorY2sum;
end

