function voisins = compute_neighbors_edges(faces, nbFaces)
% COMPUTE_NEIGHBORS_EDGES
% Cette fonction détermine les faces voisines d’un maillage triangulaire
% en se basant sur le partage d’une arête commune.
%
% Deux faces sont considérées voisines si elles partagent exactement
% deux sommets (c’est-à-dire une arête).
%
% Entrées :
%   faces    : matrice (nbFaces x 3) contenant les indices des sommets
%              de chaque face triangulaire
%   nbFaces  : nombre total de faces
%
% Sortie :
%   voisins  : tableau de cellules (nbFaces x 1)
%              voisins{i} contient les indices des faces voisines de la face i

    % Initialisation du tableau des voisins
    voisins = cell(nbFaces,1);

    
    % 1. Construction de la liste de toutes les arêtes du maillage
    %    Chaque face triangulaire possède 3 arêtes
   

    % Arêtes sous forme de paires de sommets
    E = [faces(:,[1 2]);    % arête 1–2
         faces(:,[2 3]);    % arête 2–3
         faces(:,[3 1])];   % arête 3–1

    % Tri des sommets dans chaque arête
    % (permet d’avoir une représentation unique de chaque arête)
    E = sort(E,2);

    % 2. Association de chaque arête à sa face d’origine

    % Chaque face apparaît 3 fois (une fois par arête)
    faceId = [(1:nbFaces)'; (1:nbFaces)'; (1:nbFaces)'];

    
    % 3. Tri des arêtes pour regrouper les arêtes identiques
    

    % Tri lexicographique des arêtes
    [Es, ord] = sortrows(E);

    % Réorganisation des indices de faces selon ce tri
    faceSorted = faceId(ord);

    
    % 4. Détection des arêtes partagées par deux faces
   

    k = 1;
    while k <= size(Es,1)

        % Recherche du bloc d’arêtes identiques
        kk = k;
        while kk <= size(Es,1) && all(Es(kk,:) == Es(k,:))
            kk = kk + 1;
        end

        % Nombre de faces partageant cette arête
        cnt = kk - k;

        % Si exactement deux faces partagent l’arête → voisines
        if cnt == 2
            f1 = faceSorted(k);
            f2 = faceSorted(k+1);

            % Ajout réciproque des faces voisines
            voisins{f1}(end+1) = f2;
            voisins{f2}(end+1) = f1;
        end

        % Passage au groupe suivant d’arêtes
        k = kk;
    end

   
    % 5. Suppression des doublons dans les listes de voisins
    

    for f = 1:nbFaces
        voisins{f} = unique(voisins{f});
    end
end
