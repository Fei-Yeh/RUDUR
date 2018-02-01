%% Apply RUDUR to a randomly generated 6-DIG synthetic sequence 
% User can select the input cinetic of the myocardium (k_6,4)

addpath(genpath('../'));

%% Select input cinetic of the myocardium (k_6,4)
cinetic=0.175; %=0.025, 0.075, 0.0125, 0.0175

%% Create synthetic data
[Ytot,A_GTtot,F_GTtot,nbLigtot,nbColtot,nbImtot,Ktot]=createSequence6DIG(cinetic);
Y2Dtot=reshape(Ytot,nbLigtot,nbColtot,nbImtot);

%% Zoom on the heart
rowC=34:59; % row
colC=42:77; % columns
imC=1:30; % images

fC=[2 4 6]; % sources
K=length(fC);

Y2D=Y2Dtot(rowC,colC,imC);
[nbLig,nbCol,nbIm]=size(Y2D);
nbPix=nbLig*nbCol;
Y=reshape(Y2D,nbPix,nbIm);

F_GT=F_GTtot(fC,imC);
A_GTtot2D=reshape(A_GTtot,nbLigtot,nbColtot,Ktot);
A2D_GT=A_GTtot2D(rowC,colC,fC);
A_GT=reshape(A2D_GT,nbPix,K);

%% RUDUR parameters
alpha=1;
beta=10;
gamma=3;
mu=0.1;
Q=ones(1,nbIm);
S=ones(1,nbIm);

%% Select ROI
load('../Dataset/6DIG_Synthetic/ROI_synth.mat');

disp('Display sequence');
implay(Y2D/10);

%% RUDUR
[Aest,Fest,fObj,fWLS,fPrior,fRegF,finalIter] = rudur(Y,M,alpha,beta,gamma,mu,Q,S);

%% Normalize results and ground truth
F_GTN=zeros(K,nbIm);
A_GTN=zeros(nbPix,K);
coeffGT=zeros(K,1);
for k=1:K
    coeffGT(k)=nbPix/sum(A_GT(:,k))/K;
    F_GTN(k,:)=F_GT(k,:)/coeffGT(k);
    A_GTN(:,k)=A_GT(:,k)*coeffGT(k);
end

FestN=zeros(K,nbIm);
AestN=zeros(nbPix,K);
coeff=zeros(K,1);
for k=1:K
    coeff(k)=nbPix/sum(Aest(:,k))/K;
    FestN(k,:)=Fest(k,:)/coeff(k);
    AestN(:,k)=Aest(:,k)*coeff(k);
end

%% Plot factors estimation
F_GTN(3,:)=3*F_GTN(3,:);
FestN(3,:)=3*FestN(3,:);
figure();
plot(F_GTN','--','linewidth',2);ylim([0 13]);hold on;plot(FestN','linewidth',3);%
legend('VD','VG','Myocardium','VD-RUDUR','VG-RUDUR','Myocardium-RUDUR');
set(gca, 'FontSize', 15);
title(strcat('Cinetic=',num2str(cinetic)))

%% Plot factor images estimation
AestN2D=reshape(AestN,nbLig,nbCol,K);
A_GTN2D=reshape(A_GTN,nbLig,nbCol,K);

figure();
for k=1:K
    subplot(3,K,k); imshow(squeeze(A_GTN2D(:,:,k))/2.3); 
    title(strcat('F',num2str(k),'-GT'));
    
    subplot(3,K,K+k); imshow(squeeze(M(:,:,k))); 
    title(strcat('ROI',num2str(k)));
    
    subplot(3,K,2*K+k); imshow(squeeze(AestN2D(:,:,k))/2.3);
    title(strcat('F',num2str(k),'-RUDUR'));
end





