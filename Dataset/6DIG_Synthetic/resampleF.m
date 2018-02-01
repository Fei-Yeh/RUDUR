function [F,nbIm]= resampleF(F1,dtD,dt)

[K,nbImD]=size(F1);
nbIm=floor(nbImD*dtD/dt);
pas=dt/dtD;
F=zeros(K,nbIm);
for t=1:nbIm
    F(:,t)=mean(F1(:,1+(t-1)*pas:t*pas),2);
end

end

