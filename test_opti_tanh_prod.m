clear all; close all; clc;

nb_iterations = 1e+5;

H = [1 0 1 0 0 0;
     0 1 0 1 1 0;
     0 0 1 1 0 1];

[all_icy, all_ivx] = find(H); % On prend tous les indices des 1 de la matrice H
nb_ones = length(all_icy);
nb_voisins = size(H, 1)-1;
C2V_or = zeros(3,6);

V2C_or = [-4.00437270926381 0 6.24947167113775 0 0 0;
        0 8.66001138831218 0 8.06226376998732 -11.0602405375143 0;
        0 0 6.24947167113775 8.06226376998732 0 4.98786593622922];

%% METHODE 1

C2V = C2V_or;
V2C = V2C_or;
for n = 1:nb_iterations
        %% Boucle remplissage C2V
        
        for i_c = 1:length(all_icy)
            % i_c c'est un index pour un tableau d'index #MindBlown #JeanReviensPas
            
            icy = all_icy(i_c); % Notre index (i_c) marche pour all_icy, on en déduit notre index (icy) qui marche pour H
            ivx = all_ivx(i_c);
            % Mnt, on s'interesse sur la position H(icy, ivx)
            
            % On trouve ou sont les liens relies a ce noeud
            connecteds2c = find(H(icy, :)); % on repere les liens mon noeud de parite --> noeuds de var

            % On ne retient pas les liens sur lequel on est
            connecteds2c(connecteds2c == ivx) = [];
            values_connecteds2c = V2C(icy, connecteds2c);
            
%             C2V(icy, ivx) = 2*atanh(prod(tanh(values_connecteds2c/2)));
            C2V(icy, ivx) = prod(tanh(values_connecteds2c/2), 2);
        end     
end




%% METHODE 2

C2V = C2V_or;
V2C = V2C_or;
for n = 1:nb_iterations
        %% Boucle remplissage C2V
        
        optimized_atanh_matrix = zeros(nb_ones, size(H, 2));
        for i_c = 1:length(all_icy)
            % i_c c'est un index pour un tableau d'index #MindBlown #JeanReviensPas
            
            icy = all_icy(i_c); % Notre index (i_c) marche pour all_icy, on en déduit notre index (icy) qui marche pour H
            ivx = all_ivx(i_c);
            % Mnt, on s'interesse sur la position H(icy, ivx)
            

            % On trouve ou sont les liens relies a ce noeud
            connecteds2c = find(H(icy, :)); % on repere les liens mon noeud de parite --> noeuds de var

            % On ne retient pas les liens sur lequel on est
            connecteds2c(connecteds2c == ivx) = [];
            values_connecteds2c = V2C(icy, connecteds2c);
            
            optimized_atanh_matrix(i_c, 1:length(values_connecteds2c)) = values_connecteds2c;% chaque ligne est un tableau a partir duquel on doit faire 2atanhprodtanh
            
%             C2V(icy, ivx) = 2*atanh(prod(tanh(values_connecteds2c/2)));
        end     
        merdier = prod(tanh(optimized_atanh_matrix), 2);
        
end
