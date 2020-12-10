function y = decoder( Lch, H, nb_iterations)

% Taille de la matrice H
nb_noeudsvar = size(H, 2);
nb_noeudsparite = size(H, 1);
% Reshape des LLRs si jamais + de 1 paquet
Lch_matrix = reshape(Lch, nb_noeudsvar, length(Lch)/nb_noeudsvar)';
% mise en forme de y
y = zeros(nb_noeudsparite*length(Lch)/nb_noeudsvar,1);


% icy = l'index y de H (= index sur les noeuds de C ("check" pour parite))
% ivx = l'index x (variable) de H
[all_icy, all_ivx] = find(H); % On prend tous les indices des 1 de la matrice H


for i_motdecode = 1:size(Lch_matrix, 1)
    
    mdc = Lch_matrix(i_motdecode, :); % mot de code

    LLR_repetes = repmat(mdc, nb_noeudsparite, 1); % on repete les LLRs pour pouvoir les mettre comme H (voir ligne V2C)

    %% Itération n°0 avec tout C2V a 0 donc ca laisse passer que les LLRs
    
    C2V = zeros(size(H));
    V2C = H.*LLR_repetes; % Observation comme des parites --> dans les noeuds de var
    % [y, x] = find(H);    
    
    for n = 1:nb_iterations
        %% Boucle remplissage C2V
        n;
        full(C2V);
        full(V2C);
        nb_iterations;
%         for icy = cy % Ca marche pas du tout cette merde, monsieur
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
            
            C2V(icy, ivx) = 2*atanh(prod(tanh(values_connecteds2c/2)));
        end
        %% Boucle remplissage V2C
        for i_v = 1:length(all_ivx)
            icy = all_icy(i_v); % i_c c'est un index pour un tableau d'index #MindBlown #JeanReviensPas
            ivx = all_ivx(i_v);

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

    end
%     LLRs = sum(V2C, 1); % Je somme sur les colonnes
    LLRs = sum(C2V, 1); % Je somme sur les colonnes
    msg = LLRs + mdc;
    msg = msg(end-nb_noeudsparite+1:end)'; % car code systï¿½matique

    y((1:nb_noeudsparite)+nb_noeudsparite*(i_motdecode-1)) = msg(:);

end
