function matrice_angles = compute_angle_matrix(normales, voisins)

% COMPUTE_ANGLE_MATRIX
% Cette fonction calcule la matrice des angles entre les normales
% des faces voisines d’un maillage triangulaire.

% Nombre total de faces du maillage
nbFaces = size(normales,1);

% Initialisation de la matrice des angles avec NaN
matrice_angles = nan(nbFaces);

for i = 1:nbFaces
    for j = voisins{i}
        dp = dot(normales(i,:), normales(j,:));   
        dp = max(min(dp,1),-1);
        ang = acosd(dp);
        matrice_angles(i,j) = ang;
        matrice_angles(j,i) = ang;
    end
end
end
