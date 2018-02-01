function [A] = createFactorImages()

K=11;
% 1: Veine
% 2: VD
% 3: Poumon
% 4: VG
% 5: Aorte
% 6: Myocarde
% 7: Sang
% 8: Foie
% 9: Muscle
% 10: Autres
% 11: Vessie


nbLig=90;
nbCol=120;
nbPix=nbLig*nbCol;
s1=nbLig;
s2=nbCol;


fen1=23;
fen2=25;
sigma=3;

Poids_Veine=min(PolyConcave(s1,s2,[5,35;41,28;41,30;5,37]) + PolyConcave(s1,s2,[41,28;62,46;62,48;41,30]),1);
Poids_VD=PolyConcave(s1,s2,[51,60;62,61;73,73;57,71;51,60]);
Poids_Poumon=min(PolyConcave(90,120,[33,53;47,35;55,60;46,70;33,66;33,53])+PolyConcave(90,120,[69,51;76,38;91,59;89,70;69,51]),1);
Poids_VG=PolyConcave(s1,s2,[59,55;70,56;79,72;61,64;59,55]);
Poids_Aorte=PolyConcave(s1,s2,[64,33;68,33;68,108;64,108;64,33]);
Poids_Sang=ones(s1,s2);Poids_Sang(:,[1:15 105:120])=0.7;
Poids_Foie=min(PolyConcave(s1,s2,[31,66;50,66;51,71;30,98;31,66])+PolyConcave(90,120,[51,71;70,73;30,98;51,71]),1);
Poids_Muscle=ones(s1,s2);Poids_Muscle(:,[1:15 105:120])=0.7;
Poids_Muscle=Poids_Muscle-0.5*PolyConcave(s1,s2,[31,43;90,43;90,70;30,98;31,43]);
%Poids_Myocarde=PolyConcave(s1,s2,[67,64;75,64;75,73;67,73;67,64]);
Poids_Autres=ones(s1,s2);Poids_Autres(:,[1:15 105:120])=0.7;
load('/localdata/filippma/Documents/MATLAB/Modele6DIG/roiMyo.mat');
Poids_Myocarde=double(roiMyo);

G = fspecial('gaussian',[fen1 fen2],sigma);
G2=fspecial('gaussian',[fen1 fen2],2);
Poids_Veine_G=imfilter(Poids_Veine,G,'same');Poids_Veine_G=Poids_Veine_G/sum(Poids_Veine_G(:));
Poids_VD_G = imfilter(Poids_VD,G,'same');Poids_VD_G=Poids_VD_G/sum(Poids_VD_G(:));
Poids_Poumon_G=imfilter(Poids_Poumon,G,'same');Poids_Poumon_G=Poids_Poumon_G/sum(Poids_Poumon_G(:));
Poids_VG_G=imfilter(Poids_VG,G,'same');Poids_VG_G=Poids_VG_G/sum(Poids_VG_G(:));
Poids_Aorte_G=imfilter(Poids_Aorte,G,'same');Poids_Aorte_G=Poids_Aorte_G/sum(Poids_Aorte_G(:));
Poids_Sang_G=imfilter(Poids_Sang,G,'same');Poids_Sang_G=Poids_Sang_G/sum(Poids_Sang_G(:))*0.4;
Poids_Foie_G=imfilter(Poids_Foie,G,'same');Poids_Foie_G=Poids_Foie_G/sum(Poids_Foie_G(:));
Poids_Muscle_G=imfilter(Poids_Muscle,G,'same');Poids_Muscle_G=Poids_Muscle_G/sum(Poids_Muscle_G(:))*0.4;
Poids_Myocarde_G=imfilter(Poids_Myocarde,G2,'same');Poids_Myocarde_G=Poids_Myocarde_G/sum(Poids_Myocarde_G(:));
Poids_Autres_G=imfilter(Poids_Autres,G,'same');Poids_Autres_G=Poids_Autres_G/sum(Poids_Autres_G(:))*0.4;

% Images Facteurs

A=zeros(nbPix,K);
A(:,1)=Poids_Veine_G(:);
A(:,2)=Poids_VD_G(:);
A(:,3)=Poids_Poumon_G(:);
A(:,4)=Poids_VG_G(:);
A(:,5)=Poids_Aorte_G(:);
A(:,6)=Poids_Myocarde_G(:);
A(:,7)=Poids_Sang_G(:);
A(:,8)=Poids_Foie_G(:);
A(:,9)=Poids_Muscle_G(:);
A(:,10)=Poids_Autres_G(:); % Autres
A(:,11)=0; % Vessie
end

