function facesParSommet = build_facesParSommet(faces, nbSommets)
% BUILD_FACESPARSOMMET
% Cette fonction construit, pour chaque sommet du maillage,
% la liste des faces (triangles) auxquelles il appartient.


    % Nombre total de faces (triangles)
    nbFaces = size(faces,1);

    % Initialisation d'un tableau de cellules :
    % une cellule par sommet
    facesParSommet = cell(nbSommets,1);

    % Parcours de toutes les faces du maillage
    for f = 1:nbFaces

        verts = faces(f,:);

        facesParSommet{verts(1)}(end+1) = f;

        facesParSommet{verts(2)}(end+1) = f;

        facesParSommet{verts(3)}(end+1) = f;
    end
end
