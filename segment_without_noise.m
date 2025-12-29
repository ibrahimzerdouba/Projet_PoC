function idGroupe = segment_without_noise(faces, normales, voisins, matrice_angles, params)
% SEGMENT_WITHOUT_NOISE
% Segmentation d’un maillage triangulaire "sans bruit"
% Inspiré de la Section 2 de l’article : segmentation par région (region growing)
% basée sur l’angle entre normales des faces voisines.

% Entrées :
%   faces          : (nbFaces x 3) indices des sommets de chaque face
%   normales       : (nbFaces x 3) normales unitaires
%   voisins        : cell(nbFaces,1), voisins{f} = faces voisines de f
%   matrice_angles : (nbFaces x nbFaces) angles (degrés) entre normales de faces voisines
%   params.angleMax         : seuil d’angle (degrés) pour la croissance du segment
%   params.tailleMinSegment : nb min de faces pour garder un segment tel quel
%   params.angleMaxFusion   : seuil d’angle (degrés) autorisant la fusion
%
% Sortie :
%   idGroupe : (nbFaces x 1) label du segment pour chaque face
%              0 n’est normalement pas utilisé ici (toutes les faces sont segmentées)

    % Nombre total de faces
    nbFaces = size(faces,1);

    % Paramètres de segmentation
    angleMax         = params.angleMax;
    tailleMinSegment = params.tailleMinSegment;
    angleMaxFusion   = params.angleMaxFusion;

    % ============================================================
    % 2.1 Selection of a Group of Triangles (Region Growing)
    % ============================================================

    idGroupe = zeros(nbFaces,1);

    % Compteur de segments
    groupeActuel = 0;

    % Parcours de toutes les faces : chaque face peut être une seed
    for i = 1:nbFaces
        if idGroupe(i) ~= 0
            continue;
        end

        % Nouveau segment
        groupeActuel = groupeActuel + 1;
        idGroupe(i) = groupeActuel;

        % Pile pour le DFS (region growing sur les faces)
        stack = i;

        % Croissance du segment
        while ~isempty(stack)

            % Pop d'une face
            f = stack(end);
            stack(end) = [];

            % Parcours des faces voisines de f
            for v = voisins{f}

                % Si la voisine n'est pas encore affectée
                if idGroupe(v) == 0

                    % Calcul de l’angle entre normales (robuste)
                    dp = dot(normales(f,:), normales(v,:));
                    dp = max(min(dp,1), -1);      
                    angle_fv = acosd(dp);         

                    % Critère d'appartenance au même segment
                    if angle_fv < angleMax
                        idGroupe(v) = groupeActuel;
                        stack(end+1) = v;
                    end
                end
            end
        end
    end

    % ============================================================
    % 2.2 Combination of Small Segments (Fusion des petits segments)
    % ============================================================

    % Liste des segments existants
    segments = unique(idGroupe);
    segments(segments == 0) = [];

    % Taille de chaque segment (nombre de faces)
    tailleSegment = zeros(max(segments),1);
    for s = segments'
        tailleSegment(s) = sum(idGroupe == s);
    end

    % Parcours des segments pour fusionner ceux trop petits
    for s = segments'

        % Si segment assez grand, on ne touche pas
        if tailleSegment(s) >= tailleMinSegment
            continue;
        end

        % Faces du petit segment
        facesPetit = find(idGroupe == s);

        % --------------------------------------------------------
        % 1) Trouver les segments voisins (adjacents via voisinage)
        % --------------------------------------------------------
        groupesVoisins = [];
        for f = facesPetit'
            for v = voisins{f}
                g = idGroupe(v);
                if g ~= s && g ~= 0
                    groupesVoisins(end+1) = g; 
                end
            end
        end
        groupesVoisins = unique(groupesVoisins);

        % --------------------------------------------------------
        % 2) Choisir le "meilleur" segment voisin pour fusionner :
        %    celui qui minimise l’angle de raccord entre faces voisines
        % --------------------------------------------------------
        bestGroup = -1;
        bestAngle = inf;

        for g2 = groupesVoisins'
            for f = facesPetit'
                for v = voisins{f}

                    % Ne comparer que les voisins appartenant au groupe g2
                    if idGroupe(v) == g2
                        ang = matrice_angles(f, v);

                        % Garder le plus petit angle non-NaN
                        if ~isnan(ang) && ang < bestAngle
                            bestAngle = ang;
                            bestGroup = g2;
                        end
                    end
                end
            end
        end

        % --------------------------------------------------------
        % 3) Fusion si le meilleur angle est acceptable
        % --------------------------------------------------------
        if bestGroup ~= -1 && bestAngle <= angleMaxFusion
            idGroupe(facesPetit) = bestGroup;
        end
    end
end
