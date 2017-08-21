function [ fWLS ] = compute_fWLS(Y,A,F,W2)
% Compute the value of fWLS criterion

% @input :
%   Y : input Data
%   A : Matrix of factor images (vectorized)
%   F : Matrix of factors
%   W2 : Vector containing weight for each pixel

% @return : 
%   fWLS : Value of fWLS

Y0=A*F;
errorYw=(Y0-Y);
errorY2=(errorYw.*errorYw);
errorY2sum=sum(errorY2,2);
fWLS=W2*errorY2sum;
end

