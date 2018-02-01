%% Demo file. Compute estimation errors on the simple dataset
% with different values for parameter alpha (Factor image penalization)

addpath(genpath('../'));


%% Choose dataset parameters
background=0; % 0:without background --- 1:with background
sigma=0.5; % Noise
nbSequences=100;

%% Data Info
nbRow=50;
nbCol=50;
nbPix=nbRow*nbCol;
nbIm=50;

%% Set up Data
load('../Dataset/SimpleDataset/Simple_AF.mat');
if background==1
    K=4;
else
    K=3;
    A_GT=A_GT(:,1:K);
    F_GT=F_GT(1:K,:);
end

Y_GT=A_GT*F_GT;
Yr2=sqrt(Y_GT);

%% RUDUR Parameters
tab_alpha=[0 0.2 0.5 1 2 5 10 20 50];
beta=10;
mu=0.1;
gamma=3;
Q=ones(1,nbIm);
S=ones(1,nbIm);

%% Estimation errors storage
NMAE_F=zeros(length(tab_alpha),nbSequences,K);
NMSE_F=zeros(length(tab_alpha),nbSequences,K);
NMAE_A=zeros(length(tab_alpha),nbSequences,K);
NMSE_A=zeros(length(tab_alpha),nbSequences,K);     

%% ROI Selection
priorA2D=zeros(nbRow,nbCol,K);
priorA2D(2:27,2:34,2)=1;
priorA2D(15:49,2:37,1)=1;
priorA2D(9:49,16:49,3)=1;
priorA=reshape(priorA2D,nbRow*nbCol,K);
if background==1
    %priorA(:,4)=min(sum(priorA(:,1:3),2),1); % Background
    priorA(:,4)=1;
    priorA2D=reshape(priorA,nbRow,nbCol,K);
end
M=priorA2D;

FestNGlobal=zeros(length(tab_alpha),nbSequences,K,nbIm);
AestN2DGlobal=zeros(length(tab_alpha),nbSequences,nbRow,nbCol,K);

disp(strcat('Compute mean estimation errors on',32,num2str(nbSequences),32,'sequences, with',32,num2str(length(tab_alpha)),32,'valeurs différentes pour alpha'));

for i=1:nbSequences
    for l=1:length(tab_alpha)
        alpha=tab_alpha(l);
        disp(strcat('Sequence',32,num2str(i),'/',num2str(nbSequences),32,'---',32,'alpha',32,num2str(l),'/',num2str(length(tab_alpha))));

        % Generate new data
        Y=Y_GT+normrnd(0,sigma,[nbPix nbIm]).*Yr2;

        % Apply RUDUR
        [Aest,Fest] = rudur(Y,M,alpha,beta,gamma,mu,Q,S);

        % Normalize results and ground truth to avoid scale indeterminacy
        F_GTN=zeros(K,nbIm);
        A_GTN=zeros(nbPix,K);
        coeffGT=zeros(K,1);
        for k=1:K
            coeffGT(k)=nbPix/sum(abs(A_GT(:,k)))/K;
            F_GTN(k,:)=F_GT(k,:)/coeffGT(k);
            A_GTN(:,k)=A_GT(:,k)*coeffGT(k);
        end

        FestN=zeros(K,nbIm);
        AestN=zeros(nbPix,K);
        coeff=zeros(K,1);
        for k=1:K
            coeff(k)=nbPix/sum(abs(Aest(:,k)))/K;
            FestN(k,:)=Fest(k,:)/coeff(k);
            AestN(:,k)=Aest(:,k)*coeff(k);
        end
        % De-vectorization of factor images        
        AestN2D=reshape(AestN,50,50,K);
        A_GTN2D=reshape(A_GTN,50,50,K);

        FestNGlobal(l,i,:,:)=FestN;
        AestN2DGlobal(l,i,:,:,:)=AestN2D;

        %% Compute and store estimation errors
        diffF=abs(FestN-F_GTN); diffF2=diffF.*diffF;
        diffA=abs(AestN-A_GTN); diffA2=diffA.*diffA;

        for k=1:K
            NMSE_F(l,i,k)=sum(diffF2(k,:))/sum(F_GTN(k,:).*F_GTN(k,:));
            NMAE_F(l,i,k)=sum(diffF(k,:))/sum(F_GTN(k,:));
            NMSE_A(l,i,k)=sum(diffA2(:,k))/sum(A_GTN(:,k).*A_GTN(:,k));
            NMAE_A(l,i,k)=sum(diffA(:,k))/sum(A_GTN(:,k));
        end
    end
end

%% Diagram - Alpha
NMAE_Falpha=squeeze(permute(NMAE_F,[2 3 4 1]));
NMAE_Falphamean=squeeze(mean(NMAE_Falpha,1));
bar(reshape(NMAE_Falphamean,K,length(tab_alpha)));
set(gca,'xticklabel',{'F1','F2','F3'})
set(gca,'fontsize',18);
legend('\alpha=0','\alpha=0.2','\alpha=0.5','\alpha=1','\alpha=2','\alpha=5','\alpha=10','\alpha=20','\alpha=50','Location','NorthWest');
title(strcat('Mean of NMAE on',32,num2str(nbSequences),32,'sequences with (sigma=',num2str(sigma),')'));


numEtude=5;left=1;
for l=[1 4]
figure();plot(squeeze(FestNGlobal(l,1,:,:))','LineWidth',3);
hold on; plot(F_GTN','--','LineWidth',3);set(gca,'fontsize',15);
ylim=[0 3.5];
legend('Cortex-RUDUR','Medulla-RUDUR','Pelvis-RUDUR','Cortex-GT','Medulla-GT','Pelvis-GT','Location','NorthWest');
title(strcat('alpha=',num2str(tab_alpha(l))));
end

figure();
for k=1:K
    subplot(K,4,1+4*(k-1));imshow(squeeze(M(:,:,k)));title(strcat('A',num2str(k),32,'ROI'));
    subplot(K,4,2+4*(k-1));imshow(squeeze(A_GTN2D(:,:,k)/2.3));title(strcat('A',num2str(k),32,'GT'));
    subplot(K,4,3+4*(k-1));imshow(squeeze(AestN2DGlobal(1,1,:,:,k)/2.3));title(strcat('A',num2str(k),32,'alpha=',num2str(tab_alpha(1))));
    subplot(K,4,4+4*(k-1));imshow(squeeze(AestN2DGlobal(4,1,:,:,k)/2.3));title(strcat('A',num2str(k),32,'alpha=',num2str(tab_alpha(4))));
end
