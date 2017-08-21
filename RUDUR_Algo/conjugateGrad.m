function [conjGrad]=conjugateGrad(grad,gradPrev,conjGradPrev,methodCoeff)
% Compute conjugate gradient

% @input :
%   grad : Current gradient
%   gradPrec : Gradient on the last iteration
%   conjGradPrec : Conjugate gradient on the last iteration
%   methodCoeff : Determine the conjugate gradient method used
%           1 -> Flecher-Reeves
%           2 -> Polak-Ribiere
%           3 -> Hestenes-Stiefel
%           4 -> Dai-Yuan
%           5 -> modified Polak-Ribiere


% @return : 
%   conjGrad : Current conjugate gradient

switch methodCoeff
    case 1 % Fletcher-Reeves (FR)
        coeff=(grad'*grad)/(gradPrev'*gradPrev);
    case 2 % Polak-Ribiere (PR)
        coeff=(grad'*(grad-gradPrev))/(gradPrev'*gradPrev);
    case 3 % Hestenes-Stiefel (HS)
        coeff=-(grad'*(grad-gradPrev))/(conjGradPrev'*(grad-gradPrev));
    case 4 % Dai-Yuan (DY)
        coeff=-(grad'*grad)/(conjGradPrev'*(grad-gradPrev));
    case 5 % modified Polak-Ribiere
        coeff=max(0,(grad'*(grad-gradPrev))/(gradPrev'*gradPrev));
    otherwise
        disp('Choice of conjugate gradient method not valid');
end
conjGrad=grad+coeff*conjGradPrev;

end

