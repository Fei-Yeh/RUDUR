%% Demo file. Compute estimation errors on the realistic renographies
% with different values for parameter beta (Tykhonov regularization)
% Two different matrix S are tested.

addpath(genpath('../'));

K=3;
nbIm=60;

%% CHOOSE diagonal of matrix S 
S=ones(1,nbIm);
S(1:20)=0; % Comment this line if you don't want to use S (S=Identity)

%% RUDUR Parameters
alpha=1;
tab_beta=[0,0.5,1,2,5,10,20,50]; % 
mu=0.1;
gamma=3;
Q=ones(1,nbIm);

%% Estimation errors storage
NMAE_F=zeros(length(tab_beta),6,2,K);
NMSE_F=zeros(length(tab_beta),6,2,K);
NMAE_A=zeros(length(tab_beta),6,2,K);
NMSE_A=zeros(length(tab_beta),6,2,K);

            
FestNGlobal=zeros(length(tab_beta),6,2,3,nbIm);
            
disp(strcat('Compute mean estimation errors on',32,12,32,'sequences'));
for i=1:6
    for left=1:2
        for l=1:length(tab_beta)
            beta=tab_beta(l);
            disp(strcat('Sequence',32,num2str(2*(i-1)+left),'/',num2str(12),32,'--',32,'beta',32,num2str(l),'/',num2str(length(tab_beta))));
            %% Load Data and ROI
            load(strcat('../Dataset/SyntheticRenography/Seq',num2str(i),'_',num2str(left),'.mat'));
            nbPix=nbRow*nbCol;
            
            %% Apply RUDUR
            [Aest,Fest] = rudur(Y,M,alpha,beta,gamma,mu,Q,S);

            %% Normalize results and ground truth to avoid scale indeterminacy
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
            for k=1:K
                F_GTN(k,:)=F_GTN(k,:)/sum(F_GTN(k,:))*nbIm;
                FestN(k,:)=FestN(k,:)/sum(FestN(k,:))*nbIm;
            end

            % De-vectorization of factor images        
            AestN2D=reshape(AestN,nbRow,nbCol,K);
            A_GTN2D=reshape(A_GTN,nbRow,nbCol,K);

            FestNGlobal(l,i,left,:,:)=FestN;
            
            %% Compute and store estimation errors
            diffF=abs(FestN-F_GTN); diffF2=diffF.*diffF;
            diffA=abs(AestN-A_GTN); diffA2=diffA.*diffA;

            for k=1:K
                NMSE_F(l,i,left,k)=sum(diffF2(k,:))/sum(F_GTN(k,:).*F_GTN(k,:));
                NMAE_F(l,i,left,k)=sum(diffF(k,:))/sum(F_GTN(k,:));
                NMSE_A(l,i,left,k)=sum(diffA2(:,k))/sum(A_GTN(:,k).*A_GTN(:,k));
                NMAE_A(l,i,left,k)=sum(diffA(:,k))/sum(A_GTN(:,k));
            end
        end
    end
end

%% Tychonof diagramme baton
NMAE_Fbeta=squeeze(permute(NMAE_F,[2 3 4 1]));
NMAE_Fbetamean=squeeze(mean(NMAE_Fbeta,1));
bar(reshape(NMAE_Fbetamean,6,length(tab_beta)));
set(gca,'xticklabel',{'Cortex-R','Cortex-L','Medulla-R','Medulla-L','Pelvis-R','Pelvis-L'})
set(gca,'fontsize',18);
legend('\beta=0','\beta=0.5','\beta=1','\beta=2','\beta=5','\beta=10','\beta=20','\beta=50','Location','NorthWest');

numEtude=5;left=1;
for l=[1 6 8]
figure();plot(squeeze(FestNGlobal(l,numEtude,left,:,:))','LineWidth',3);
hold on; plot(F_GTN','--','LineWidth',3);set(gca,'fontsize',15);
ylim=[0 3.5];
legend('Cortex-RUDUR','Medulla-RUDUR','Pelvis-RUDUR','Cortex-GT','Medulla-GT','Pelvis-GT');
title(strcat('Beta=',num2str(tab_beta(l))));
end
