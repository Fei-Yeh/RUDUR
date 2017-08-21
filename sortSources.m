function [AestNnew,FestNnew]=sortSources(AestN,FestN,A_GTN,F_GTN)
    % Sort estimate sources in the same order than ground truth
    
    nbPix=size(AestN,1);
    nbIm=size(FestN,2);
    K=size(AestN,2);
    
    AestNnew=zeros(nbPix,K);
    FestNnew=zeros(K,nbIm);
    error=zeros(K,K);
    for k1=1:K
        for k2=1:K
            error(k1,k2)=sum(abs(F_GTN(k1,:)-FestN(k2,:)))/sum(F_GTN(k1,:));
        end
    end
    for k1=1:K
        [~,kmin]=min(error(k1,:));
        AestNnew(:,k1)=AestN(:,kmin);
        FestNnew(k1,:)=FestN(kmin,:);
    end
end

