%% Test of the robustness of RUDUR to the ROI selection
% ROI are translated horizontaly/vertical and/or dilated leading to 
% 5832 set of ROI

addpath(genpath('../'));

%% Select sequence
numSequence=4;
left=1;

%% Load data
load(strcat('../Dataset/SyntheticRenography/Seq',num2str(numSequence),'_',num2str(left),'.mat'));
priorA2D=M;
nbPix=nbRow*nbCol;

NMAE_Ftot=zeros(3,3,3,3,3,3,3,3,3,K);
NMSE_Ftot=zeros(3,3,3,3,3,3,3,3,3,K);
NMAE_Atot=zeros(3,3,3,3,3,3,3,3,3,K);
NMSE_Atot=zeros(3,3,3,3,3,3,3,3,3,K);
se = strel('disk',1,0);


disp('WARNING : These results take a while to be computed');
disp('RUDUR has to be launch 5832 times on the same sequence');

for dil1=2:3
for dil2=2:3
for dil3=2:3
    for dV1=1:3
    for dV2=1:3
    for dV3=1:3
        for dH1=1:3
        for dH2=1:3
        for dH3=1:3
           disp(strcat(num2str(dH3+3*(dH2-1)+9*(dH1-1) + 27*(dV3-1) +81*(dV2-1) + 243*(dV1-1) + 729*(dil3-2) + 1458*(dil2-2)+ 2916*(dil1-2)),'/',num2str(8*3^6)));
           M=priorA2D;     
           M(:,:,1)=imdilate(M(:,:,1),se);   
           M(:,:,2)=imerode(M(:,:,2),se);
            if dil1==1
                M(:,:,1)=imerode(M(:,:,1),se);
            elseif dil1==3
                M(:,:,1)=imdilate(M(:,:,1),se);
            end

            if dil2==1
                M(:,:,2)=imerode(M(:,:,2),se);
            elseif dil2==3
                M(:,:,2)=imdilate(M(:,:,2),se);
            end

            if dil3==1
                M(:,:,3)=imerode(M(:,:,3),se);
            elseif dil3==3
                M(:,:,3)=imdilate(M(:,:,3),se);
            end

            
            if dV1==1
                M(:,1:end-1,1)=M(:,2:end,1);
            elseif dV1==3
                M(:,2:end,1)=M(:,1:end-1,1);
            end

            if dV2==1
                M(:,1:end-1,2)=M(:,2:end,2);
            elseif dV2==3
                M(:,2:end,2)=M(:,1:end-1,2);
             end  

            if dV3==1
                M(:,1:end-1,3)=M(:,2:end,3);
            elseif dV3==3
                M(:,2:end,3)=M(:,1:end-1,3);
            end            


            if dH1==1
                M(1:end-1,:,1)=M(2:end,:,1);
            elseif dH1==3
                M(2:end,:,1)=M(1:end-1,:,1);
            end

            if dH2==1
                M(1:end-1,:,2)=M(2:end,:,2);
            elseif dH2==3
                M(2:end,:,2)=M(1:end-1,:,2);
            end    
            if dH3==1
                M(1:end-1,:,3)=M(2:end,:,3);
            elseif dH3==3
                M(2:end,:,3)=M(1:end-1,:,3);
            end           
            M(:,:,2)=imdilate(M(:,:,2),se); 
            M(:,1:end-1,1)=M(:,2:end,1); 
            
            %% RUDUR parameters
            alpha=1;
            beta=10;
            gamma=3;
            mu=0.1;
            Q=ones(1,nbIm);
            S=ones(1,nbIm);

            %% Methode
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

            Fdiff=abs(F_GTN-FestN);
            Adiff=abs(A_GTN-AestN);

            NMAE_Ftot(dil1,dil2,dil3,dV1,dV2,dV3,dH1,dH2,dH3,:)=mean(Fdiff,2);
            NMAE_Atot(dil1,dil2,dil3,dV1,dV2,dV3,dH1,dH2,dH3,:)=mean(Adiff,1);
            NMSE_Ftot(dil1,dil2,dil3,dV1,dV2,dV3,dH1,dH2,dH3,:)=mean(Fdiff.*Fdiff,2);
            NMSE_Atot(dil1,dil2,dil3,dV1,dV2,dV3,dH1,dH2,dH3,:)=mean(Adiff.*Adiff,1);

            diffF=abs(FestN-F_GTN); diffF2=diffF.*diffF;
            diffA=abs(AestN-A_GTN); diffA2=diffA.*diffA;

            for k=1:K
                NMSE_Atot(dil1,dil2,dil3,dV1,dV2,dV3,dH1,dH2,dH3,k)=sum(diffA2(:,k))/sum(A_GTN(:,k).*A_GTN(:,k));
                NMAE_Atot(dil1,dil2,dil3,dV1,dV2,dV3,dH1,dH2,dH3,k)=sum(diffA(:,k))/sum(A_GTN(:,k));
                NMSE_Ftot(dil1,dil2,dil3,dV1,dV2,dV3,dH1,dH2,dH3,k)=sum(diffF2(k,:))/sum(F_GTN(k,:).*F_GTN(k,:));
                NMAE_Ftot(dil1,dil2,dil3,dV1,dV2,dV3,dH1,dH2,dH3,k)=sum(diffF(k,:))/sum(F_GTN(k,:));
            end

        end
        end
        end
    end
    end
    end
end
end
end




a1=NMAE_Ftot(2:3,2:3,2:3,:,:,:,:,:,:,1);
a2=NMAE_Ftot(2:3,2:3,2:3,:,:,:,:,:,:,2);
a3=NMAE_Ftot(2:3,2:3,2:3,:,:,:,:,:,:,3);

s1=NMAE_Atot(2:3,2:3,2:3,:,:,:,:,:,:,1);
s2=NMAE_Atot(2:3,2:3,2:3,:,:,:,:,:,:,2);
s3=NMAE_Atot(2:3,2:3,2:3,:,:,:,:,:,:,3);

disp('Mean estimation errors (Factor):');
disp(strcat('NMAE F1=',num2str(mean(a1(:)),'%.3f'),'+-',num2str(sqrt(var(a1(:))),'%.3f')));
disp(strcat('NMAE F2=',num2str(mean(a2(:)),'%.3f'),'+-',num2str(sqrt(var(a2(:))),'%.3f')));
disp(strcat('NMAE F3=',num2str(mean(a3(:)),'%.3f'),'+-',num2str(sqrt(var(a3(:))),'%.3f')));

disp(strcat('NMAE A1=',num2str(mean(s1(:)),'%.3f'),'+-',num2str(sqrt(var(s1(:))),'%.3f')));
disp(strcat('NMAE A2=',num2str(mean(s2(:)),'%.3f'),'+-',num2str(sqrt(var(s2(:))),'%.3f')));
disp(strcat('NMAE A3=',num2str(mean(s3(:)),'%.3f'),'+-',num2str(sqrt(var(s3(:))),'%.3f')));


NMAE_Ftot1DK=reshape(NMAE_Ftot(2:3,2:3,2:3,:,:,:,:,:,:,:),8*3^6*3,1);
NMAE_Ftot1D=reshape(NMAE_Ftot(2:3,2:3,2:3,:,:,:,:,:,:,:),8*3^6,3);

NMAE_Atot1DK=reshape(NMAE_Atot(2:3,2:3,2:3,:,:,:,:,:,:,:),8*3^6*3,1);
NMAE_Atot1D=reshape(NMAE_Atot(2:3,2:3,2:3,:,:,:,:,:,:,:),8*3^6,3);

Origin=char(zeros(8*3^6*3,7));
for k=1:8*3^6
    Origin(k,:)='Cortex ';
end
for k=8*3^6+1:8*2*3^6
    Origin(k,:)='Medulla';
end
for k=8*2*3^6+1:8*3*3^6
    Origin(k,:)='Pelvis ';
end

% F
figure();
boxplot(NMAE_Ftot1DK,Origin);
ylim([0 2]);set(gca,'fontsize',24);
title('Distribution of NMAE of factors (5832 differents set of ROI are used)');

%A
figure();
boxplot(NMAE_Atot1DK,Origin);
ylim([0 2]);set(gca,'fontsize',24);
title('Distribution of NMAE of factor images (5832 differents set of ROI are used)');
