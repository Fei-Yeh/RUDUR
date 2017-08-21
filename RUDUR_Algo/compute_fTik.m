function [ fTik ] = compute_fTik(F)
% Compute the value of fTik criterion

% @input :
%   F : Matrix of factors

% @return : 
%   fTik : Value of fTik criterion

diffF=F(:,2:end)-F(:,1:end-1); % difference between two consecutive images
diffF2=diffF.*diffF;

fTik=sum(sum(diffF2));

end

