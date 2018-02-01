%% Demo File. Apply RUDUR to a 6-DIG clinical sequence

%% Select a sequence (from 1 to 3)
numSeq=1;

load(strcat('../Dataset/6DIG/Seq',num2str(numSeq),'.mat'));

disp('Display sequence');
implay(reshape(Y,nbLig,nbCol,nbIm)/10);

%% Param RUDUR
alpha=1;
beta=10;
gamma=3;
mu=0.1;
Q=ones(1,nbIm);
S=ones(1,nbIm);
K=3;



%% Apply RUDUR
[Aest,Fest,fObj,fWLS,fPrior,fRegF,finalIter] = rudur(Y,M,alpha,beta,gamma,mu,Q,S);


%% Normalize results
FestN=zeros(K,nbIm);
AestN=zeros(nbPix,K);
coeff=zeros(K,1);
for k=1:K
    coeff(k)=nbPix/sum(Aest(:,k))/K;
    FestN(k,:)=Fest(k,:)/coeff(k);
    AestN(:,k)=Aest(:,k)*coeff(k);
end

%% Plot factor estimation 
figure();
plot(FestN','linewidth',3);%
legend('VD','VG','Myocarde');
set(gca, 'FontSize', 30);
title('Factor estimation (RUDUR)');

%% Plot ROI and factorimage estimation
AestN2D=reshape(AestN,nbLig,nbCol,K);
figure();
for k=1:K   
    subplot(2,K,k); imshow(squeeze(M(:,:,k))); 
    title(strcat('ROI',num2str(k)));
    
    subplot(2,K,K+k); imshow(squeeze(AestN2D(:,:,k))/2.3);
    title(strcat('F',num2str(k),'-RUDUR'));
end
