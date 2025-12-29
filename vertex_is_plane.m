function isPlan = vertex_is_plane(v, facesParSommet, matrice_angles, anglePlane)
% VERTEX_IS_PLANE
% Détermine si un sommet peut être considéré comme appartenant
% à une zone plane du maillage.
%
% Principe :
% - On considère toutes les faces incidentes au sommet v.
% - Le sommet est dit "plan" si l’angle entre les normales de
%   toutes les paires de faces voisines autour de ce sommet
%   est inférieur ou égal à un seuil anglePlane.
%
% Entrées :
%   v              : indice du sommet testé
%   facesParSommet : cell(nbSommets,1), facesParSommet{v} = faces incidentes à v
%   matrice_angles : (nbFaces x nbFaces), angle (en degrés) entre normales
%                    des faces voisines (NaN si non voisines)
%   anglePlane     : seuil d’angle (degrés) pour considérer une zone plane
%
% Sortie :
%   isPlan : booléen
%            true  -> sommet considéré comme plan
%            false -> sommet non plan

    % Faces incidentes au sommet v
    Fv = facesParSommet{v};

    % Si moins de 3 faces, la planéité n'est pas fiable
    if numel(Fv) < 3
        isPlan = false;
        return;
    end

    % Vérification de toutes les paires de faces autour du sommet
    for a = 1:numel(Fv)
        for b = a+1:numel(Fv)

            % Indices des deux faces
            f1 = Fv(a);
            f2 = Fv(b);

            % Angle entre les normales des deux faces
            ang = matrice_angles(f1, f2);

            % Si l'angle est défini et dépasse le seuil -> pas plan
            if ~isnan(ang) && ang > anglePlane
                isPlan = false;
                return;
            end
        end
    end

    % Si toutes les paires respectent le seuil -> sommet plan
    isPlan = true;
end
