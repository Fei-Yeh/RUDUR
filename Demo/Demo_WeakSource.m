% Demo file. Compute estimation errors on several sequences of the simple
% dataset
% User can change the value of Q13 to Q19 

addpath(genpath('../'));

%% Choose dataset parameters
sigma=0.5; % Noise
nbSequences=100;

%% Data Info
nbRow=50;
nbCol=50;
nbPix=nbRow*nbCol;
nbIm=50;

%% RUDUR Parameters
alpha=1;
beta=10;
mu=0.1;
gamma=3;
Q=ones(1,nbIm);
S=ones(1,nbIm);

%% Change Q :
Q(13:19)=3;%1,sqrt(2),sqrt(3),2 or 3
alpha=alpha*sum(Q.*Q)/nbIm; % Normalization of alpha
beta=beta*sum(Q.*Q)/nbIm; % Normalization of beta


%% Set up Data
load('../Dataset/SimpleDataset/Weak_AF.mat');
K=3;

Y_GT=A_GT*F_GT;
Yr2=sqrt(Y_GT);


%% Estimation errors storage
NMAE_F=zeros(nbSequences,K);
NMSE_F=zeros(nbSequences,K);
NMAE_A=zeros(nbSequences,K);
NMSE_A=zeros(nbSequences,K);

%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%



%% ROI Selection
priorA2D=zeros(nbRow,nbCol,K);
priorA2D(18:42,13:37,2)=1;
priorA2D(15:49,2:37,1)=1;
priorA2D(9:49,16:49,3)=1;
priorA=reshape(priorA2D,nbRow*nbCol,K);

M=priorA2D;


disp(strcat('Compute mean estimation errors on',32,num2str(nbSequences),32,'sequences'));
for i=1:nbSequences
    disp(strcat('Sequence',32,num2str(i),'/',num2str(nbSequences)));
    %% Generate new data
    Y=Y_GT+normrnd(0,sigma,[nbPix nbIm]).*Yr2;

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
    % De-vectorization of factor images        
    AestN2D=reshape(AestN,50,50,K);
    A_GTN2D=reshape(A_GTN,50,50,K);

    %% Compute and store estimation errors
    diffF=abs(FestN-F_GTN); diffF2=diffF.*diffF;
    diffA=abs(AestN-A_GTN); diffA2=diffA.*diffA;

    for k=1:K
        NMSE_F(i,k)=sum(diffF2(k,:))/sum(F_GTN(k,:).*F_GTN(k,:));
        NMAE_F(i,k)=sum(diffF(k,:))/sum(F_GTN(k,:));
        NMSE_A(i,k)=sum(diffA2(:,k))/sum(A_GTN(:,k).*A_GTN(:,k));
        NMAE_A(i,k)=sum(diffA(:,k))/sum(A_GTN(:,k));
    end 
end

disp('Display last sequence');
implay(reshape(Y,nbRow,nbCol,nbIm)/5);

%% Display estimation errors
disp('Display estimation errors');
disp(strcat('NMSE : F1=',num2str(mean(NMSE_F(:,1)),'%.3f'),' --- F2=',num2str(mean(NMSE_F(:,2)),'%.3f'),' --- F3=',num2str(mean(NMSE_F(:,3)),'%.3f')));
disp(strcat('NMAE : F1=',num2str(mean(NMAE_F(:,1)),'%.3f'),' --- F2=',num2str(mean(NMAE_F(:,2)),'%.3f'),' --- F3=',num2str(mean(NMAE_F(:,3)),'%.3f')));
disp('---');
disp(strcat('NMSE : A1=',num2str(mean(NMSE_A(:,1)),'%.3f'),' --- A2=',num2str(mean(NMSE_A(:,2)),'%.3f'),' --- A3=',num2str(mean(NMSE_A(:,3)),'%.3f')));
disp(strcat('NMAE : A1=',num2str(mean(NMAE_A(:,1)),'%.3f'),' --- A2=',num2str(mean(NMAE_A(:,2)),'%.3f'),' --- A3=',num2str(mean(NMAE_A(:,3)),'%.3f')));

%% Plot last results
% Factors
disp('Plot estimation results');
figure();plot(FestN(:,:)','LineWidth',3);hold on;plot(F_GTN','--','LineWidth',3);
set(gca,'fontsize',15);
if K==3
    legend('F1-RUDUR','F2-RUDUR','F3-RUDUR','F1-GT','F2-GT','F3-GT','Location','NorthWest');
else
    legend('F1-RUDUR','F2-RUDUR','F3-RUDUR','F4-RUDUR','F1-GT','F2-GT','F3-GT','F4-GT','Location','NorthWest');
end
axis([0 50 0 3.5]);
title(strcat('Factor estimation (Fx-RUDUR) vs ground truth (Fx-GT), sigma=',num2str(sigma)));
xlabel('time');
ylabel('Activity');

% Factor images
figure();
for k=1:K
    subplot(3,K,k); imshow(squeeze(A_GTN2D(:,:,k))/2.3); 
    title(strcat('F',num2str(k),'-GT'));
    
    subplot(3,K,K+k); imshow(squeeze(M(:,:,k))); 
    title(strcat('ROI',num2str(k)));
    
    subplot(3,K,2*K+k); imshow(squeeze(AestN2D(:,:,k))/2.3);
    title(strcat('F',num2str(k),'-RUDUR'));
end
