function [F] = capture(Fc,dt,r,h)
    
    % Proba de desintegration pendant dt
    t12=13.3*60*60;
    pRadio=log(2)/t12*dt;
    
    % Proba de passage d'une desintegration dans la collim
    phi=r/h;
    pCollim=(1-cos(phi))/2;
    
    F=Fc*pRadio*pCollim;
end

