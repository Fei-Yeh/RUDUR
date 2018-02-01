% Demo file. Compute estimation errors on the realistic renographies
% dataset

addpath(genpath('../'));

%% Data Info
K=3;
nbIm=60;

%% Estimation errors storage
NMAE_F=zeros(6,2,K);
NMSE_F=zeros(6,2,K);
NMAE_A=zeros(6,2,K);
NMSE_A=zeros(6,2,K);

%% RUDUR Parameters
alpha=1;
beta=10;
mu=0.1;
gamma=3;
Q=ones(1,nbIm);
S=ones(1,nbIm);

disp(strcat('Compute mean estimation errors on',32,'12',32,'sequences'));
for i=1:6
    for left=1:2
        disp(strcat('Sequence',32,num2str(2*(i-1)+left),'/',num2str(12)));
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


        %% Compute and store estimation errors
        diffF=abs(FestN-F_GTN); diffF2=diffF.*diffF;
        diffA=abs(AestN-A_GTN); diffA2=diffA.*diffA;

        for k=1:K
            NMSE_F(i,left,k)=sum(diffF2(k,:))/sum(F_GTN(k,:).*F_GTN(k,:));
            NMAE_F(i,left,k)=sum(diffF(k,:))/sum(F_GTN(k,:));
            NMSE_A(i,left,k)=sum(diffA2(:,k))/sum(A_GTN(:,k).*A_GTN(:,k));
            NMAE_A(i,left,k)=sum(diffA(:,k))/sum(A_GTN(:,k));
        end
    end
end


%% Display estimation errors
disp('Display estimation errors');

disp('Right Kidney');
disp(strcat('NMSE : F1=',num2str(mean(NMSE_F(:,1,1)),'%.3f'),' --- F2=',num2str(mean(NMSE_F(:,1,2)),'%.3f'),' --- F3=',num2str(mean(NMSE_F(:,1,3)),'%.3f')));
disp(strcat('NMAE : F1=',num2str(mean(NMAE_F(:,1,1)),'%.3f'),' --- F2=',num2str(mean(NMAE_F(:,1,2)),'%.3f'),' --- F3=',num2str(mean(NMAE_F(:,1,3)),'%.3f')));
disp('---');
disp(strcat('NMSE : A1=',num2str(mean(NMSE_A(:,1,1)),'%.3f'),' --- A2=',num2str(mean(NMSE_A(:,1,2)),'%.3f'),' --- A3=',num2str(mean(NMSE_A(:,1,3)),'%.3f')));
disp(strcat('NMAE : A1=',num2str(mean(NMAE_A(:,1,1)),'%.3f'),' --- A2=',num2str(mean(NMAE_A(:,1,2)),'%.3f'),' --- A3=',num2str(mean(NMAE_A(:,1,3)),'%.3f')));

disp('Left Kidney');
disp(strcat('NMSE : F1=',num2str(mean(NMSE_F(:,2,1)),'%.3f'),' --- F2=',num2str(mean(NMSE_F(:,2,2)),'%.3f'),' --- F3=',num2str(mean(NMSE_F(:,2,3)),'%.3f')));
disp(strcat('NMAE : F1=',num2str(mean(NMAE_F(:,2,1)),'%.3f'),' --- F2=',num2str(mean(NMAE_F(:,2,2)),'%.3f'),' --- F3=',num2str(mean(NMAE_F(:,2,3)),'%.3f')));
disp('---');
disp(strcat('NMSE : A1=',num2str(mean(NMSE_A(:,2,1)),'%.3f'),' --- A2=',num2str(mean(NMSE_A(:,2,2)),'%.3f'),' --- A3=',num2str(mean(NMSE_A(:,2,3)),'%.3f')));
disp(strcat('NMAE : A1=',num2str(mean(NMAE_A(:,2,1)),'%.3f'),' --- A2=',num2str(mean(NMAE_A(:,2,2)),'%.3f'),' --- A3=',num2str(mean(NMAE_A(:,2,3)),'%.3f')));
