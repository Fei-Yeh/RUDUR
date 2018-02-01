function [Poids] = PolyConcave(s1,s2,Coord)

Poids = zeros (s1,s2);
Coord(:,1)=Coord(:,1)-4;
Coord(:,2)=Coord(:,2)-19;
Coord(:,[1,2])=Coord(:,[2,1]);
% Boucle sur les colonnes
for j=min(Coord(:,1)):max(Coord(:,1))
    % Etape 1 : On recherche le ou les points ou il y a intersection avec
    % une droite. On note leur numero de colonnes x1 et x2
    prec=0;
    l=[];
    aff=0;
    for k=2:size(Coord,1)
        if Coord(k,1)>Coord(k-1,1)
            c1=Coord(k-1,1);c2=Coord(k,1);l1=Coord(k-1,2);l2=Coord(k,2);
        else
            c1=Coord(k,1);c2=Coord(k-1,1);l1=Coord(k,2);l2=Coord(k-1,2);
        end
        if (j>=c1 && j<=c2)
            %Il y a donc une arete intersecté, on cherche donc le numero de
            %ligne 
            if c1==c2
                Poids(j,l1:l2)=1;
                aff=1;
                break;
            else
                prec=prec+1;
                lCur=l1+round(((j-c1)/(c2-c1))*(l2-l1));
                l=[l,lCur];
            end
            % Attention au cas ou Coord(k,1)=Coord(k-1,1)
            % Attention au calcul d'entiers/flottants
        end
    end
    if aff==0
        l=unique(l);
        if size(l,2)==2
            Poids(j,min(l):max(l))=1;
        elseif size(l,2)==1
            Poids(j,l(1))=1;
        elseif size(l,2)>2
            disp('Erreur dans la recherche d abscisse');
        end
    end
end