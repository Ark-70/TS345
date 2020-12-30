function [y, nb_iterations_done] = decoder_V2( Lch, H, MIN_SUM_ON)

% Taille de la matrice H
nb_noeudsvar = size(H, 2);
nb_noeudsparite = size(H, 1);

% icy = l'index Y de H (= index sur les noeuds de C_heck (pour parite))
% ivx = l'index X (V_ariable) de H
[all_icy, all_ivx] = find(H); % On prend tous les indices des 1 de la matrice H


mdc = Lch'; % mot de code

LLR_repetes = repmat(mdc, nb_noeudsparite, 1); % on repete les LLRs pour pouvoir les mettre comme H (voir ligne V2C)

%% Initialisation tout a 0

C2V = sparse(size(H, 1), size(H, 2));
V2C = sparse(size(H, 1), size(H, 2));

i_ite = 1;

final_LLRs_chap = zeros(1, nb_noeudsparite);
last_final_LLRs_chap = final_LLRs_chap;

while(1)
    
    % for icy = cy % Ca marche pas du tout ce truc, ca bug

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


        if(MIN_SUM_ON)
            C2V(icy, ivx) = prod(sign(values_connecteds2c)) * min(abs(values_connecteds2c));
        else % SUM_PRODUCT
            C2V(icy, ivx) = 2*atanh(prod(tanh(values_connecteds2c/2)));    
        end

    end


    %% STOP adaptatif pour le nombre d'itérations
    final_LLRs_chap = C2VMatrix_to_final_LLRs(C2V, mdc, nb_noeudsparite);

    % Si on a trouve le bon mdc, meme pas la peine d'essayer de converger plus
    x_chap = double(final_LLRs_chap < 0);
    syndrome = x_chap.'*H;
    if(sum(syndrome) == 0)
        nb_iterations_done = i_ite;
        'break !!';
        break
    end 

    if(sum(final_LLRs_chap ~= final_LLRs_chap)==0) % Si on a converge et qu'on n'avance plus
        nb_iterations_done = i_ite;
        'break !';
        break
    end
    %%
    
    i_ite = i_ite + 1;
end

y = C2VMatrix_to_final_LLRs(C2V, mdc, nb_noeudsparite);

end
