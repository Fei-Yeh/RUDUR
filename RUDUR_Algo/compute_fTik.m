function [ fTik ] = compute_fTik(F,S)
% Compute the value of fTik criterion

% @input :
%   F : Matrix of factors

% @return : 
%   fTik : Value of fTik criterion

diffF=F(:,2:end)-F(:,1:end-1); % difference between two consecutive images
diffF2=diffF.*diffF;
diffF2s=sum(diffF2,1);

fTik=sum(diffF2s.*S(2:end));

end

