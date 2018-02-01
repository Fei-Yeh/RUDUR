% Demo file. Apply RUDUR on a simple dataset and display/plot results
% User can modify noise intensity and can also add a background

addpath(genpath('../'));

%% Choose dataset parameters
background=0; % 0:without background --- 1:with background
sigma=1; % Noise

%% Data Info
nbRow=50;
nbCol=50;
nbPix=nbRow*nbCol;
nbIm=50;

%% Generate Data
load('../Dataset/SimpleDataset/Simple_AF.mat');

disp('Randomly Generate Data');
if background==1
    K=4;
else
    K=3;
    A_GT=A_GT(:,1:K);
    F_GT=F_GT(1:K,:);
end

Y_GT=A_GT*F_GT;
Yr2=sqrt(Y_GT);
Y=Y_GT+normrnd(0,sigma,[nbPix nbIm]).*Yr2;

%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%

%% RUDUR Parameters
alpha=1;
beta=10;
mu=0.1;
gamma=3;
Q=ones(1,nbIm);
S=ones(1,nbIm);

%% ROI Selection
priorA2D=zeros(nbRow,nbCol,K);
priorA2D(2:27,2:34,2)=1;
priorA2D(15:49,2:37,1)=1;
priorA2D(9:49,16:49,3)=1;
priorA=reshape(priorA2D,nbRow*nbCol,K);
if background==1
    priorA(:,4)=1;
    priorA2D=reshape(priorA,nbRow,nbCol,K);
end
M=priorA2D;

disp('Display sequence');
implay(reshape(Y,nbRow,nbCol,nbIm)/5);

%% Apply RUDUR
disp('Apply RUDUR');
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
        


%% Plot results
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

%% Compute and display estimation errors
disp('Display mean estimation errors');
diffF=abs(FestN-F_GTN); diffF2=diffF.*diffF;
diffA=abs(AestN-A_GTN); diffA2=diffA.*diffA;

NMAE_F=zeros(1,K);
NMSE_F=zeros(1,K);
NMAE_A=zeros(1,K);
NMSE_A=zeros(1,K);

for k=1:K
    NMSE_F(k)=sum(diffF2(k,:))/sum(F_GTN(k,:).*F_GTN(k,:));
    NMAE_F(k)=sum(diffF(k,:))/sum(F_GTN(k,:));
    NMSE_A(k)=sum(diffA2(:,k))/sum(A_GTN(:,k).*A_GTN(:,k));
    NMAE_A(k)=sum(diffA(:,k))/sum(A_GTN(:,k));
end

if K==3
    disp(strcat('NMSE : F1=',num2str(NMSE_F(1),'%.3f'),' --- F2=',num2str(NMSE_F(2),'%.3f'),' --- F3=',num2str(NMSE_F(3),'%.3f')));
    disp(strcat('NMAE : F1=',num2str(NMAE_F(1),'%.3f'),' --- F2=',num2str(NMAE_F(2),'%.3f'),' --- F3=',num2str(NMAE_F(3),'%.3f')));
    disp('---');
    disp(strcat('NMSE : A1=',num2str(NMSE_A(1),'%.3f'),' --- A2=',num2str(NMSE_A(2),'%.3f'),' --- A3=',num2str(NMSE_A(3),'%.3f')));
    disp(strcat('NMAE : A1=',num2str(NMAE_A(1),'%.3f'),' --- A2=',num2str(NMAE_A(2),'%.3f'),' --- A3=',num2str(NMAE_A(3),'%.3f')));
else
    disp(strcat('NMSE : F1=',num2str(NMSE_F(1),'%.3f'),' --- F2=',num2str(NMSE_F(2),'%.3f'),' --- F3=',num2str(NMSE_F(3),'%.3f'),' --- F4=',num2str(NMSE_F(4),'%.3f')));
    disp(strcat('NMAE : F1=',num2str(NMAE_F(1),'%.3f'),' --- F2=',num2str(NMAE_F(2),'%.3f'),' --- F3=',num2str(NMAE_F(3),'%.3f'),' --- F4=',num2str(NMAE_F(4),'%.3f')));
    disp('---');
    disp(strcat('NMSE : A1=',num2str(NMSE_A(1),'%.3f'),' --- A2=',num2str(NMSE_A(2),'%.3f'),' --- A3=',num2str(NMSE_A(3),'%.3f'),' --- A4=',num2str(NMSE_A(4),'%.3f')));
    disp(strcat('NMAE : A1=',num2str(NMAE_A(1),'%.3f'),' --- A2=',num2str(NMAE_A(2),'%.3f'),' --- A3=',num2str(NMAE_A(3),'%.3f'),' --- A4=',num2str(NMAE_A(4),'%.3f')));
end
    


