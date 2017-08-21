addpath RUDUR_Algo/
addpath .

%% Parameters of the dataset
backg=0; % add a fourth background factor if backg=1;
sigma=1; % Noise on the dataset

%% Parameters of RUDUR algorithm
alpha=1; % Wegiht of the fPrior criterion
beta=10; % Weight of the fRegF criterion
gamma=3; % Weighting data fidelity constant
mu=0.001;  % Parameter of the smoothed-l1 norm (sigma==0 => l1 norm)

%% Load ground truth
load('dataAF.mat'); % load ground truth matrix A and F

%% Build dataset
if backg==0 % Delete background
    A=A(:,1:3);
    F=F(1:3,:);
end

K=size(A,2);
nbPix=size(A,1);
nbIm=size(F,2);
nbRow=50;nbCol=50;

Y=A*F; % Build matrix Y with A and F (linear mixing model)
Yr2=sqrt(Y);
Y=Y+normrnd(0,sigma,[nbPix nbIm]).*Yr2; %% Add Noise

%% Choice of ROIs
M_2D=zeros(nbRow,nbCol,K);
M_2D(2:27,2:34,2)=1;
M_2D(15:49,2:37,1)=1;
M_2D(9:49,16:49,3)=1;
M_1D=reshape(M_2D,nbRow*nbCol,K);
if backg==1
    M_1D(:,4)=1;
    M_2D=reshape(M_1D,nbRow,nbCol,K);
end

%% Apply RUDUR
[Aest,Fest,fObj,fWLS,fROI,fTik,finalIter] = rudur(Y,M_2D,alpha,beta,gamma,mu);

%% Normalize A and counter normalize F, and sort sources
[A_N,F_N]=normalize(A,F);
[Aest_N,Fest_N]=normalize(Aest,Fest);
[Aest_N,Fest_N]=sortSources(Aest_N,Fest_N,A_N,F_N);

%% Un-vectorized factor images
Aest_N2D=reshape(Aest_N,nbRow,nbCol,K);
A_N2D=reshape(A_N,nbRow,nbCol,K);

%% Compare results obtained with ground truth
% Display data 
%implay(reshape(Y/4,50,50,50));

% Factor curves / TACs
figure();plot(Fest_N');hold on; plot(F_N','--'); 
if K==4
    legend('Fest1','Fest2','Fest3','Fest4','F1','F2','F3','F4');
else
    legend('Fest1','Fest2','Fest3','F1','F2','F3');
end
title('Ground Truth and TACs obtained with RUDUR ');

% Factor images
figure();
for k=1:K
    subplot(K,3,3*(k-1)+1);imshow(M_2D(:,:,k)); title(strcat('ROI',32,num2str(k)));
    subplot(K,3,3*(k-1)+2);imshow(A_N2D(:,:,k)/2.3);title(strcat('A',num2str(k)));
    subplot(K,3,3*(k-1)+3);imshow(Aest_N2D(:,:,k)/2.3);title(strcat('Aest',num2str(k),32,'(RUDUR)'));
end

