function [priorA,priorA2D] = buildPrior(M_2D)
% Build prior for matrix of factor images thanks to matrix of masks M


% @input :
%   M_2D : Matrix of masks (ROIs)

% @return : 
%   priorA : Matrix of euclidean distance of each pixel to each mask
%   (vectorized)
%   priorA2D : Matrix of euclidean distance of each pixel to each mask
%   (not-vectorized)

nbRow=size(M_2D,1);
nbCol=size(M_2D,2);
K=size(M_2D,3);

%% For each pixel, compute euclidean distance to each mask
priorA2D=zeros(nbRow,nbCol,K);
for k=1:K
    priorA2D(:,:,k)=bwdist(squeeze(M_2D(:,:,k)));
end

%% Vectorize
priorA=reshape(priorA2D,nbRow*nbCol,K);

end

