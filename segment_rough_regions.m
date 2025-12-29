function idRough = segment_rough_regions(voisins, idPlane, idSmooth, matrice_angles, params)
% SEGMENT_ROUGH_REGIONS
% Segmentation des régions "rough" (rugueuses / bruitées) dans un maillage.

% Entrées :
%   voisins        : cell(nbFaces,1), voisins{f} = indices des faces voisines de f
%   idPlane        : (nbFaces x 1) étiquettes des régions planes (0 si non plane)
%   idSmooth       : (nbFaces x 1) étiquettes des régions smooth (0 si non smooth)
%   matrice_angles : (nbFaces x nbFaces) angles (en degrés) entre normales
%                    pour les faces voisines (NaN sinon)
%   params.angleRough   : seuil d'angle (degrés) pour connecter deux faces "rough"
%   params.minRoughSize : taille minimale (nb faces) pour garder un segment rough
%
% Sortie :
%   idRough : (nbFaces x 1)
%            idRough(f) = identifiant du segment rough auquel appartient la face f
%            0 si la face n'appartient à aucun segment rough

    % Nombre total de faces
    nbFaces = numel(voisins);

    % Paramètres
    angleRough   = params.angleRough;
    minRoughSize = params.minRoughSize;

    % Masques : faces déjà classées comme planes ou smooth
    isPlane  = (idPlane  ~= 0);
    isSmooth = (idSmooth ~= 0);

    % Faces restantes : ni plane ni smooth => candidates rough
    reste = (~isPlane) & (~isSmooth);

    % Initialisation des labels rough
    idRough = zeros(nbFaces,1);

    % Compteur de segments rough
    seg = 0;

    % ============================================================
    % Parcours de toutes les faces : chaque face restante peut être une seed
    % ============================================================
    for seed = 1:nbFaces

        % Si la face n'est pas dans "reste" ou déjà étiquetée rough, on saute
        if ~reste(seed) || idRough(seed) ~= 0
            continue;
        end

        % Nouveau segment rough
        seg = seg + 1;

        % Pile (DFS) de faces à explorer
        stackF = seed;

        % On affecte la seed au segment
        idRough(seed) = seg;

        % ========================================================
        % Region growing au niveau des faces (DFS)
        % ========================================================
        while ~isempty(stackF)

            % Pop : dernière face ajoutée
            f = stackF(end);
            stackF(end) = [];

            % Parcours des voisins topologiques de la face f
            for nb = voisins{f}

                % Conditions de base :
                % - le voisin doit être une face "reste"
                % - et non déjà affecté à un segment rough
                if ~reste(nb) || idRough(nb) ~= 0
                    continue;
                end

                % Angle entre normales (si NaN -> pas de valeur exploitable)
                ang = matrice_angles(f, nb);
                if isnan(ang)
                    continue;
                end

                % Critère flexible :
                % si l'angle est inférieur au seuil rough, on connecte
                if ang <= angleRough
                    idRough(nb) = seg;
                    stackF(end+1) = nb; %#ok<AGROW> (croissance dynamique OK ici)
                end
            end
        end

        % ========================================================
        % Post-traitement : suppression des segments rough trop petits
        % (souvent du bruit ou des artefacts isolés)
        % ========================================================
        if minRoughSize > 0 && sum(idRough == seg) < minRoughSize
            idRough(idRough == seg) = 0;  % annuler le segment
            seg = seg - 1;                % garder des labels compacts
        end
    end
end
