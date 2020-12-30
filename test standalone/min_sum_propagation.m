function y = min_sum_propagation( Lch, H, nb_iterations)

% Taille de la matrice H
nb_noeudsvar = size(H, 2);
nb_noeudsparite = size(H, 1);
% Reshape des LLRs si jamais + de 1 paquet
Lch_matrix = reshape(Lch, nb_noeudsvar, length(Lch)/nb_noeudsvar)';
% mise en forme de y
y = zeros(nb_noeudsparite*length(Lch)/nb_noeudsvar,1);


% icy = l'index Y de H (= index sur les noeuds de C_heck (pour parite))
% ivx = l'index X (V_ariable) de H
[all_icy, all_ivx] = find(H); % On prend tous les indices des 1 de la matrice H

for i_motdecode = 1:size(Lch_matrix, 1)

    mdc = Lch_matrix(i_motdecode, :); % mot de code

    LLR_repetes = repmat(mdc, nb_noeudsparite, 1); % on repete les LLRs pour pouvoir les mettre comme H (voir ligne V2C)

    %% It�ration n�0 avec tout C2V a 0 donc ca laisse passer que les LLRs

    C2V = sparse(size(H, 1), size(H, 2));
%     V2C = H.*LLR_repetes; % Observation comme des parites --> dans les noeuds de var
    V2C = sparse(size(H, 1), size(H, 2));
    % [y, x] = find(H);

    for n = 1:nb_iterations

%         for icy = cy % Ca marche pas du tout cette merde, ca bug

        %% Boucle remplissage V2C
        for i_v = 1:length(all_ivx)
            icy = all_icy(i_v); % i_v c'est un index pour un tableau d'index #MindBlown #JeanReviensPas
            ivx = all_ivx(i_v);
            % On s'interesse a la position H(icy, ivx) et on veut calculer V2C(icy, ivx) dans ce tour de boucle

            % Il y a un lien entre 2 noeuds ici
            connecteds2v = find(H(:, ivx)); % on repere les liens mon noeud de var --> noeuds de parites

            % On ne retient pas les liens sur lequel on est
            connecteds2v(connecteds2v == icy) = [];

            % LLR du bit du canal qui est arrive par ce noeud de var
            somme = mdc(ivx);
            for i_sum =1:length(connecteds2v)
                somme = somme + C2V(connecteds2v(i_sum), ivx);
            end
            V2C(icy, ivx) = somme;
        end

        %% Boucle remplissage C2V
        for i_c = 1:length(all_icy) % i_c c'est un index pour un tableau d'index #MindBlown #JeanReviensPas
            
            % On s'interesse a la position H(icy, ivx) et on veut calculer V2C(icy, ivx) dans ce tour de boucle
            icy = all_icy(i_c);
            ivx = all_ivx(i_c);
            
            % On trouve ou sont les liens relies a ce noeud de parite (check)
            connecteds2c = find(H(icy, :)); % on repere les liens noeuds de var --> mon noeud de parite

            % On ne retient pas le lien sur lequel on est
            connecteds2c(connecteds2c == ivx) = [];
            values_connecteds2c = V2C(icy, connecteds2c);

%             C2V(icy, ivx) = 2*atanh(prod(tanh(values_connecteds2c/2)));
            C2V(icy, ivx) = prod(sign(values_connecteds2c)) * min(abs(values_connecteds2c));
        end
              
        
        if(nb_iterations == 'ADAPTATIF')
                % On voit si C2V a assez converge
                final_LLRs_chap = C2VMatrix_to_final_LLRs(C2V, mdc, nb_noeudsparite)
                x_chap = double(final_LLRs_chap < 0)
                syndrome = x_chap.'*H;
                if(syndrome == 0)
                    n
                    break
                end 
        end
        
    end
    
    msg = C2VMatrix_to_final_LLRs(C2V, mdc, nb_noeudsparite)

    y((1:nb_noeudsparite)+nb_noeudsparite*(i_motdecode-1)) = msg(:);

end
