function idSmooth = segment_smooth_regions(sommets, faces, normales, facesParSommet, idPlane, params)
% SEGMENT_SMOOTH_REGIONS
% Segmentation des régions "smooth" (surfaces lisses / courbure régulière)
% sur un maillage triangulaire, en présence de bruit.

% Entrées :
%   sommets        : (nbSommets x 3) coordonnées 3D
%   faces          : (nbFaces x 3) indices sommets des triangles
%   normales       : (nbFaces x 3) normales unitaires des faces
%   facesParSommet : cell(nbSommets,1), faces incidentes à chaque sommet
%   idPlane        : (nbFaces x 1) étiquettes des zones planes (0 sinon)
%   params.alpha, params.beta       : poids pour rho12
%   params.curvThresh               : seuil de "boundary edge"
%   params.minSmoothSize            : taille minimale d’un segment smooth (nb faces)
%
% Sortie :
%   idSmooth : (nbFaces x 1) labels des segments smooth (0 si non smooth)

    % -----------------------------
    % Paramètres
    % -----------------------------
    alpha = params.alpha;
    beta  = params.beta;
    curvThresh    = params.curvThresh;
    minSmoothSize = params.minSmoothSize;

    nbFaces   = size(faces,1);
    nbSommets = size(sommets,1);

    % Masque des faces déjà planes
    isPlaneFace = (idPlane ~= 0);

    %% ============================================================
    % (A) Construction des arêtes + association arête -> faces
    %     + indexation edgesOfFace (3 arêtes par face)
    % ============================================================

    % Liste de toutes les arêtes (3 par face)
    allE = [faces(:,[1 2]);   % arête (1-2)
            faces(:,[2 3]);   % arête (2-3)
            faces(:,[3 1])];  % arête (3-1)

    % Tri des sommets dans chaque arête pour une représentation unique
    allE = sort(allE,2);

    % Pour chaque arête listée, on retient l’ID de la face correspondante
    faceId = [(1:nbFaces)'; (1:nbFaces)'; (1:nbFaces)'];

    % Tri lexicographique des arêtes -> arêtes identiques deviennent consécutives
    [Es, ord] = sortrows(allE);
    faceSorted = faceId(ord);

    % edgeVerts : liste unique des arêtes du maillage (paires [a b])
    % edgeFaces : faces adjacentes à chaque arête :
    %            [f1 f2] si arête interne, [f 0] si arête frontière (mesh boundary)
    edgeVerts = [];
    edgeFaces = [];

    k = 1;
    while k <= size(Es,1)
        kk = k;
        while kk <= size(Es,1) && all(Es(kk,:) == Es(k,:))
            kk = kk + 1;
        end

        cnt = kk - k;                    % nb d'occurrences de cette arête
        edgeVerts(end+1,:) = Es(k,:);    %#ok<AGROW>

        if cnt == 2
            % arête partagée par 2 faces -> arête interne
            edgeFaces(end+1,:) = [faceSorted(k) faceSorted(k+1)]; %#ok<AGROW>
        else
            % cnt = 1 typiquement -> arête frontière
            edgeFaces(end+1,:) = [faceSorted(k) 0]; %#ok<AGROW>
        end

        k = kk;
    end

    % edgesOfFace(f,:) = indices des 3 arêtes associées à la face f
    [~,~,edgeIdxAll] = unique(allE,'rows');
    edgesOfFace = reshape(edgeIdxAll, [nbFaces 3]);

    % edgesParSommet{v} = liste des arêtes incidentes au sommet v
    edgesParSommet = cell(nbSommets,1);
    for e = 1:size(edgeVerts,1)
        a = edgeVerts(e,1);
        b = edgeVerts(e,2);
        edgesParSommet{a}(end+1) = e;
        edgesParSommet{b}(end+1) = e;
    end
    for v = 1:nbSommets
        edgesParSommet{v} = unique(edgesParSommet{v});
    end

    %% ============================================================
    % (B) Calcul d'une "courbure robuste" sur les arêtes : rho12
    %     rho_e = theta/L  +  delta_rho via rhoFace
    % ============================================================

    nbEdges = size(edgeVerts,1);

    % Courbure locale par arête : rho_e = theta / L
    rho_e = zeros(nbEdges,1);
    for e = 1:nbEdges
        f1 = edgeFaces(e,1);
        f2 = edgeFaces(e,2);

        % Arête frontière : pas de 2ème face -> courbure infinie (séparation)
        if f2 == 0
            rho_e(e) = inf;
            continue;
        end

        % Angle entre normales (sécurisé numériquement)
        dp = max(min(dot(normales(f1,:), normales(f2,:)), 1), -1);
        theta = acos(dp); % radians

        % Longueur de l’arête
        L = norm(sommets(edgeVerts(e,1),:) - sommets(edgeVerts(e,2),:));

        % Courbure "theta/L" (éviter division par 0)
        rho_e(e) = theta / max(L, eps);
    end

    % Courbure moyenne par face : moyenne des courbures de ses 3 arêtes
    rhoFace = mean(rho_e(edgesOfFace), 2, "omitnan");
    rhoFace(isnan(rhoFace)) = 0;

    % Variation de courbure entre 2 faces adjacentes à l’arête : delta_rho
    delta_rho = zeros(nbEdges,1);
    for e = 1:nbEdges
        f1 = edgeFaces(e,1);
        f2 = edgeFaces(e,2);

        if f2 == 0
            delta_rho(e) = inf;
        else
            delta_rho(e) = abs(rhoFace(f1) - rhoFace(f2));
        end
    end

    % Courbure robuste combinée (pondérée)
    rho12 = alpha*rho_e + beta*delta_rho;

    % Une arête est "boundary" si la courbure dépasse le seuil
    isBoundaryEdge = (rho12 > curvThresh);

    %% ============================================================
    % (C) Region growing "smooth"
    % - seed = face non plane, non labellisée, et sans boundary edges
    % - propagation via les sommets "smooth"
    % ============================================================

    idSmooth = zeros(nbFaces,1);
    seg = 0;

    for seed = 1:nbFaces

        % Exclure faces planes ou déjà segmentées
        if isPlaneFace(seed) || idSmooth(seed) ~= 0
            continue;
        end

        % La seed doit être "smooth" : aucune de ses arêtes n'est boundary
        if any(isBoundaryEdge(edgesOfFace(seed,:)))
            continue;
        end

        % Nouveau segment smooth
        seg = seg + 1;
        idSmooth(seed) = seg;

        % Pile de sommets à explorer (propagation par sommets)
        stackV = faces(seed,:);
        visitedV = false(nbSommets,1);

        while ~isempty(stackV)

            % Pop sommet
            v = stackV(end);
            stackV(end) = [];

            % Si déjà visité, on saute
            if visitedV(v)
                continue;
            end
            visitedV(v) = true;

            % Critère sommet "smooth" :
            % aucune arête incidente à ce sommet n'est boundary
            if any(isBoundaryEdge(edgesParSommet{v}))
                continue;
            end

            % Ajouter toutes les faces incidentes à ce sommet
            for f = facesParSommet{v}

                % Exclure faces planes et faces déjà labellisées smooth
                if isPlaneFace(f) || idSmooth(f) ~= 0
                    continue;
                end

                % Affecter la face au segment smooth
                idSmooth(f) = seg;

                % Ajouter les 3 sommets de cette face à explorer
                stackV(end+1:end+3) = faces(f,:); 
            end
        end

        % Supprimer les segments smooth trop petits (souvent bruit/artefacts)
        if sum(idSmooth == seg) < minSmoothSize
            idSmooth(idSmooth == seg) = 0;
            seg = seg - 1;
        end
    end
end
