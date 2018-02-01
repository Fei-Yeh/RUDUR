function [W] = buildWeight(D,gamma)
% Build the matrix of weight W relaxing data fidelity


% @input :
%   D : Matrix of distance for each pixel to each mask
%   gamma : Weighting data fidelity constant

% @return : 
%   W : matrix of weight W relaxing data fidelity

epsilonW=0.001;
nbPix=size(D,1);
W=zeros(1,nbPix);
for i=1:nbPix
    W(i)=max(1/(1+gamma*min(D(i,:))),epsilonW);
end

end

