function idPlane = segment_plane_regions(faces, facesParSommet, matrice_angles, params)
% SEGMENT_PLANE_REGIONS
% Segmentation des régions planes dans un maillage triangulaire.
%
% Idée :
% - On part de sommets "seed" considérés comme plans (vertex_is_plane).
% - On propage (region growing) la région en ajoutant les sommets voisins
%   qui sont eux aussi "plans".
% - Toutes les faces incidentes aux sommets validés sont affectées au segment.
% - On rejette les segments trop petits (minPlaneSize).
%
% Entrées :
%   faces           : (nbFaces x 3) indices des sommets de chaque triangle
%   facesParSommet  : cell(nbSommets,1), facesParSommet{v} = faces incidentes à v
%   matrice_angles  : (nbFaces x nbFaces), angle entre normales des faces voisines (NaN sinon)
%   params.anglePlane   : seuil (en degrés) pour décider si un sommet est "plan"
%   params.minPlaneSize : taille minimale d’un segment plan (en nb de faces)
%
% Sortie :
%   idPlane : (nbFaces x 1)
%            idPlane(f) = identifiant du segment plan auquel appartient la face f
%            0 si la face n’est pas dans un segment plan

    % Nombre de faces et de sommets
    nbFaces   = size(faces,1);
    nbSommets = numel(facesParSommet);

    % Paramètres de segmentation
    anglePlane   = params.anglePlane;     % seuil d'angle (degrés)
    minPlaneSize = params.minPlaneSize;   % taille minimale du segment (nb faces)

    % Vecteur d’étiquettes : 0 = non segmenté / non plan
    idPlane = zeros(nbFaces,1);

    % Compteur de segments (1,2,3,...)
    seg = 0;

    % ============================================================
    % Parcours de tous les sommets : chaque sommet peut servir de seed
    % ============================================================
    for v0 = 1:nbSommets

        % 1) Le sommet seed doit être considéré "plan"
        if ~vertex_is_plane(v0, facesParSommet, matrice_angles, anglePlane)
            continue;
        end

        % 2) Faces incidentes au sommet v0
        F0 = facesParSommet{v0};

        % Si pas de faces autour, ou si toutes déjà segmentées -> on ignore
        if isempty(F0) || all(idPlane(F0) ~= 0)
            continue;
        end

        % On démarre un nouveau segment
        seg = seg + 1;

        % ========================================================
        % Region growing au niveau des sommets
        % stackV : pile de sommets à explorer (DFS)
        % visitedV : sommets déjà testés (éviter boucles)
        % ========================================================
        stackV = v0;
        visitedV = false(nbSommets,1);
        visitedV(v0) = true;

        % Masque des faces appartenant au segment courant
        facesSegMask = false(nbFaces,1);

        % ========================================================
        % Boucle de propagation : tant qu'il reste des sommets à explorer
        % ========================================================
        while ~isempty(stackV)

            % Pop (DFS) : on prend le dernier sommet ajouté
            v = stackV(end);
            stackV(end) = [];

            % 3) Ajouter toutes les faces incidentes à ce sommet dans le segment
            for f = facesParSommet{v}
                if idPlane(f) == 0
                    idPlane(f) = seg;
                    facesSegMask(f) = true;
                end
            end

            % 4) Récupérer toutes les faces du segment actuel
            facesSeg = find(facesSegMask);

            % 5) Récupérer tous les sommets de ces faces (voisinage "global" du segment)
            vertsSeg = unique(faces(facesSeg,:));

            % 6) Tester ces sommets : si sommet "plan", on l'ajoute à la pile
            for vv = vertsSeg'
                if ~visitedV(vv)

                    % On ne propage que vers les sommets eux aussi "plans"
                    if vertex_is_plane(vv, facesParSommet, matrice_angles, anglePlane)
                        stackV(end+1) = vv;
                    end

                    % Marquer comme visité (testé)
                    visitedV(vv) = true;
                end
            end
        end

        % ========================================================
        % Post-traitement : suppression des petits segments (bruit)
        % ========================================================
        if sum(idPlane == seg) < minPlaneSize
            idPlane(idPlane == seg) = 0;  % on annule l'étiquette
            seg = seg - 1;                % on "retire" ce segment du compteur
        end
    end
end
