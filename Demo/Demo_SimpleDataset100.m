% Demo file. Compute estimation errors on several random sequences of the
% simple dataset
% User can add a background

addpath(genpath('../'));

%% Choose dataset parameters
background=1; % 0:without background --- 1:with background
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


%% Estimation errors storage
NMAE_F=zeros(nbSequences,K);
NMSE_F=zeros(nbSequences,K);
NMAE_A=zeros(nbSequences,K);
NMSE_A=zeros(nbSequences,K);

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


%% Display estimation errors
disp(strcat('Display mean estimation errors (on',32,num2str(nbSequences),32,'sequences)'));

if K==3
    disp(strcat('NMSE : F1=',num2str(mean(NMSE_F(:,1)),'%.3f'),' --- F2=',num2str(mean(NMSE_F(:,2)),'%.3f'),' --- F3=',num2str(mean(NMSE_F(:,3)),'%.3f')));
    disp(strcat('NMAE : F1=',num2str(mean(NMAE_F(:,1)),'%.3f'),' --- F2=',num2str(mean(NMAE_F(:,2)),'%.3f'),' --- F3=',num2str(mean(NMAE_F(:,3)),'%.3f')));
    disp('---');
    disp(strcat('NMSE : A1=',num2str(mean(NMSE_A(:,1)),'%.3f'),' --- A2=',num2str(mean(NMSE_A(:,2)),'%.3f'),' --- A3=',num2str(mean(NMSE_A(:,3)),'%.3f')));
    disp(strcat('NMAE : A1=',num2str(mean(NMAE_A(:,1)),'%.3f'),' --- A2=',num2str(mean(NMAE_A(:,2)),'%.3f'),' --- A3=',num2str(mean(NMAE_A(:,3)),'%.3f')));
else
    disp(strcat('NMSE : F1=',num2str(mean(NMSE_F(:,1)),'%.3f'),' --- F2=',num2str(mean(NMSE_F(:,2)),'%.3f'),' --- F3=',num2str(mean(NMSE_F(:,3)),'%.3f'),' --- F4=',num2str(mean(NMSE_F(:,4)),'%.3f')));
    disp(strcat('NMAE : F1=',num2str(mean(NMAE_F(:,1)),'%.3f'),' --- F2=',num2str(mean(NMAE_F(:,2)),'%.3f'),' --- F3=',num2str(mean(NMAE_F(:,3)),'%.3f'),' --- F4=',num2str(mean(NMAE_F(:,4)),'%.3f')));
    disp('---');
    disp(strcat('NMSE : A1=',num2str(mean(NMSE_A(:,1)),'%.3f'),' --- A2=',num2str(mean(NMSE_A(:,2)),'%.3f'),' --- A3=',num2str(mean(NMSE_A(:,3)),'%.3f'),' --- A4=',num2str(mean(NMSE_A(:,4)),'%.3f')));
    disp(strcat('NMAE : A1=',num2str(mean(NMAE_A(:,1)),'%.3f'),' --- A2=',num2str(mean(NMAE_A(:,2)),'%.3f'),' --- A3=',num2str(mean(NMAE_A(:,3)),'%.3f'),' --- A4=',num2str(mean(NMAE_A(:,4)),'%.3f')));
end
