function normales = compute_normals(sommets, faces)
% COMPUTE_NORMALS
% Cette fonction calcule le vecteur normal unitaire de chaque face
% triangulaire d’un maillage STL.
%
% Entrées :
%   sommets : matrice (nbSommets x 3) contenant les coordonnées (x,y,z)
%             des sommets du maillage
%   faces   : matrice (nbFaces x 3) contenant les indices des sommets
%             de chaque face triangulaire
%
% Sortie :
%   normales : matrice (nbFaces x 3)
%              normales(i,:) = vecteur normal unitaire de la face i

    
    nbFaces = size(faces,1);

    normales = zeros(nbFaces, 3);

    for i = 1:nbFaces   
        v1 = sommets(faces(i,2),:) - sommets(faces(i,1),:);
        v2 = sommets(faces(i,3),:) - sommets(faces(i,1),:);
        n = cross(v1, v2);

        normales(i,:) = n / norm(n);
    end
end
