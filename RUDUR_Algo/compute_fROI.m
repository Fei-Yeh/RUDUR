function [ fROI ] = compute_fROI(A,D2,mu)
% Compute the value of fROI criterion

% @input :
%   A : Matrix of factor images (vectorized)
%   D2 : Square of the  distance to masks matrix
%   mu : Parameter of the smoothed-l1 norm (sigma==0 => l1 norm)

% @return : 
%   fROI: Value of fROI criterion


if mu==0 % Compute l1-norm
    fROI=sum(sum(D2.*abs(A)));
else % Compute smoothed l1-norm
    D4=D2.*D2;
    fROI=sum(sum(sqrt((D4.*A.*A)+mu*mu)-mu));
end

end

