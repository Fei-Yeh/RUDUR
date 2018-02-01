function [Y,A_GT,F_GT,nbLig,nbCol,nbIm,K]=createSequence6DIG(cinetic)

       
kMyoIn=cinetic;
kMyoOut=0.5*0.5;

%% Meta-Parameters
SNR=1;
dtD=0.1; % En secondes
dt=2; % En secondes - Echantillonage image finale
nbImD=450*20*2;
t12=13.3*60*60;
N0=10^8*t12/log(2)*SNR;
lbd=0.5; % Input function decay

r=0.005; % Rayon du capteur/collim : 2mm
h=0.3; % Hauteur patient/collim : 40cm

nbLig=90;
nbCol=120;

T=(nbImD-1)*dtD;
T2=(nbImD/2-1)*dtD;

%% Build Input function
% Basal
Qin(1:nbImD/2)=N0*exp(-lbd*(0:dtD:T2))/sum(exp(-lbd*(0:dtD:T2)));
% Insuline
Qin(nbImD/2+1:nbImD)=N0*exp(-lbd*(0:dtD:T2))/sum(exp(-lbd*(0:dtD:T2)));

%% Build quantities in each compartiment from model and input function
[Fc] = buildFactorfromModelComp(Qin,dtD,nbImD,kMyoIn,kMyoOut);
K=size(Fc,1);

%% Resample
[FcR,nbIm]=resampleF(Fc,dtD,dt);

%% Radioactivity + Collimation
F_GT=capture(FcR,dt,r,h);

A_GT=createFactorImages();

Y_GT=A_GT*F_GT;
Y=poissrnd(Y_GT);

end
