%% Test of the robustness of RUDUR to the parameters selection
% Different values are set to alpha, beta, and gamma

addpath(genpath('../'));

nbIm=60;

%% Select array of parameters
tab_alpha=[1,2,5,10,20,50];
tab_beta=[1,2,5,10];
tab_gamma=[0.5 1 3 5];
tab_mu=[0.1];

Q=ones(1,nbIm);
S=ones(1,nbIm);

z1=length(tab_alpha);
z2=length(tab_beta);
z3=length(tab_gamma);
z4=length(tab_mu);
ztot=z1*z2*z3*z4;

%% For estimation errors storage
NMAE_F=zeros(z1,z2,z3,z4,6,2,3);
NMSE_F=zeros(z1,z2,z3,z4,6,2,3);
NMAE_A=zeros(z1,z2,z3,z4,6,2,3);
NMSE_A=zeros(z1,z2,z3,z4,6,2,3);

for numEtude=1:6
    for left=1:2
        for i1=1:z1
            alpha=tab_alpha(i1);
            for i2=1:z2
                beta=tab_beta(i2);
                for i3=1:z3 % Gamma
                    gamma=tab_gamma(i3);
                    for i4=1:z4 %mu
                        mu=tab_mu(i4);
                        
                        disp(strcat('Sequence',32,num2str(2*(numEtude-1)+left),'/','12',32,'---',32,'param',32,num2str(i4+z4*(i3-1)+z4*z3*(i2-1)+z4*z3*z2*(i1-1)),'/',num2str(ztot)));

                        %% Load data
                        load(strcat('../Dataset/SyntheticRenography/Seq',num2str(numEtude),'_',num2str(left),'.mat'));
                        nbPix=nbRow*nbCol;
                        se = strel('disk',1,0);
                        M(:,:,2)=imdilate(M(:,:,2),se); 
                        M(:,1:end-1,1)=M(:,2:end,1); 


                        %% Apply RUDUR
                        [Aest,Fest,fObj,fWLS,fPrior,fRegF,finalIter] = rudur(Y,M,alpha,beta,gamma,mu,Q,S);

                        
                        %% Compute estimation errors
                       
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

                        for k=1:K
                            F_GTN(k,:)=F_GTN(k,:)/sum(F_GTN(k,:))*nbIm;
                            FestN(k,:)=FestN(k,:)/sum(FestN(k,:))*nbIm;
                        end

               

                        diffF=abs(FestN-F_GTN); diffF2=diffF.*diffF;
                        diffA=abs(AestN-A_GTN); diffA2=diffA.*diffA;


                        for k=1:K
                            NMSE_A(i1,i2,i3,i4,numEtude,left,k)=sum(diffA2(:,k))/sum(A_GTN(:,k).*A_GTN(:,k));
                            NMAE_A(i1,i2,i3,i4,numEtude,left,k)=sum(diffA(:,k))/sum(A_GTN(:,k));
                            NMSE_F(i1,i2,i3,i4,numEtude,left,k)=sum(diffF2(k,:))/sum(F_GTN(k,:).*F_GTN(k,:));
                            NMAE_F(i1,i2,i3,i4,numEtude,left,k)=sum(diffF(k,:))/sum(F_GTN(k,:));
                        end
                    end
                end
            end
        end
    end
end

NMSE_A1D=reshape(NMSE_A,z1*z2*z3*z4,6,2,3);
NMSE_F1D=reshape(NMSE_F,z1*z2*z3*z4,6,2,3);
NMAE_A1D=reshape(NMAE_A,z1*z2*z3*z4,6,2,3);
NMAE_F1D=reshape(NMAE_F,z1*z2*z3*z4,6,2,3);

Origin=char(zeros(ztot,6,2,3,6));

for i=1:3
    for j=1:2
        for k=1:6         
            for l=1:ztot
                Origin(l,k,j,i,:)=strcat('P',num2str(k),'O',num2str(j),'F',num2str(i));
            end
        end
    end
end

Origin1DF1=reshape(Origin(:,:,:,1,:),ztot*6*2,6);
Origin1DF2=reshape(Origin(:,:,:,2,:),ztot*6*2,6);
Origin1DF3=reshape(Origin(:,:,:,3,:),ztot*6*2,6);
NMAE_F1DstackF1=reshape(NMAE_F1D(:,:,:,1),ztot*6*2,1);
NMAE_F1DstackF2=reshape(NMAE_F1D(:,:,:,2),ztot*6*2,1);
NMAE_F1DstackF3=reshape(NMAE_F1D(:,:,:,3),ztot*6*2,1);

figure();
boxplot(NMAE_F1DstackF1,Origin1DF1);
ylim([0 1]);set(gca,'fontsize',24);
title('Distribution of NMAE of 1st factor (96 set of parameters are used)');

figure();
boxplot(NMAE_F1DstackF2,Origin1DF2);
ylim([0 1]);set(gca,'fontsize',24);
title('Distribution of NMAE of 2nd factor (96 set of parameters are used)');

figure();
boxplot(NMAE_F1DstackF3,Origin1DF3);
ylim([0 1]);set(gca,'fontsize',24);
title('Distribution of NMAE of 3rd factor (96 set of parameters are used)');



