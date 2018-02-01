function [Fs] = buildFactorfromModelComp(Qin,dt,nbIm,kMyoIn,kMyoOut)

T=(nbIm-1)*dt;

K=11;
% 1: Veine
% 2: VD
% 3: Poumon
% 4: VG
% 5: Aorte
% 6: Sang
% 7: Vessie
% 8: Foie
% 9: Muscle
% 10: Myocarde
% 11: Ailleurs

Fs=zeros(K,nbIm);
echange=zeros(11,11); % echange(i,j) : de i vers j

% 1 - Veine -> VD
echange(1,2)=0.3;

% 2 - VD -> Poumon
echange(2,3)=0.4;

% 3 - Poumon -> VG
echange(3,4)=0.4;

% 4 - VG -> Aorte
echange(4,5)=0.5*0.95;

% 4 - VG -> Myocarde
echange(4,6)=kMyoIn;

% 5 - Aorte -> Sang
echange(5,7)=1;

% 6 - Myocarde -> Sang
echange(6,7)=kMyoOut;

% 7 - Sang -> Vessie
echange(7,11)=0.001;

% 7 - Sang -> Foie
echange(7,8)=0.003;

% 7 - Sang -> Muscle
echange(7,9)=0.003;

% 7 - Sang -> Veine
echange(7,1)=0.005;

% 7 - Sang -> VD
echange(7,2)=0.03;

% 7 - Sang-> Ailleurs
echange(7,10)=0.5;

% 8 - Foie -> Sang
echange(8,7)=0.01;

% 9 - Muscle -> Sang
echange(9,7)=0.01;

% 11 -Ailleurs -> Sang
echange(10,7)=0.3;

%t=1
Fs(1,1)=Qin(1);
for t=2:nbIm
    if t==nbIm/2+1
        echange(4,10)=echange(4,10)*2;
        echange(10,5)=echange(10,5)*2;
    end
% Pour chaque compartiment regarder ce qu'il se passe à t+1...
    for k=1:K
        Fs(k,t)=0;
        for l=1:K
            Fs(k,t)=Fs(k,t-1)+sum(echange(:,k).*Fs(:,t-1))*dt-sum(echange(k,:))*Fs(k,t-1)*dt;
        end
    end
    Fs(1,t)=Fs(1,t)+Qin(t);
end



% A la toute fin on fait une décroissance radioactive !
t12=13.3*3600; % Demi-vie de 13,3 heures
lbd=log(2)/t12;
decroi=exp(-lbd*(0:nbIm-1)*dt);
Fs=Fs.*repmat(decroi,K,1);

end

