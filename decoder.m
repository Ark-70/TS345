function y = decoder( Lch, H )




%% premiere passe avec tout C2V a 0 donc ca laisse passer que les LLRs
nb_noeudsvar = size(H, 2);
nb_noeudsparite = size(H, 1);
Lch_matrix = reshape(Lch, nb_noeudsvar, length(Lch)/nb_noeudsvar)';

%% Je veux mettre tous les noeuds de parites (carres) dans les noeuds de var (ronds)
% sachant qu'on comprend tous les n_parites = 0, + sauf ceux qu'on rajoute
% qui sont les observations du canal

% Dans les 6 de Lc_matrix
% je dois refourguer Lc(1) dans V2C(1, 1) V2C(2, 1) V2C(3, 1) etc.,
% exceptes ceux qui sont a 0 dans H

% On a des msg de 3, codes avec redondance en mot de code de 6

% premier mot de code

% On veut dupliquer le premier bit du mot de code dans tous les noeuds de
% parite correspondant

% Apres on fait des +(XOR) pour avoir la valeur des noeuds de variables


% voir slide 43
y = zeros(nb_noeudsparite*length(Lch)/nb_noeudsvar,1);
for i_motdecode = 1:size(Lch_matrix, 1)
    
    mdc = Lch_matrix(i_motdecode, :);

    LLR_repetes = repmat(mdc, 3, 1);

    %% Initialisation V2C/C2V
    
    C2V = H.*LLR_repetes; % Observation comme des parites dans les noeuds de var
    V2C = C2V; % Noeuds de var dans les parites
    % [y, x] = find(H);
    nb_iterations = 4; % pour chaque mot de code
    
    for n = 1:nb_iterations
        %% Boucle remplissage C2V
        for icy = 1:nb_noeudsparite
            for ivx = 1:nb_noeudsvar
                if(H(icy, ivx) == 1)
                    % Il y a un lien entre 2 noeuds ici

                    % On trouve ou sont les liens relies a ce noeud
                    connecteds2c = find(H(icy, :)); % on repere les liens mon noeud de parite --> noeuds de var

                    % On ne retient pas les liens sur lequel on est
                    connecteds2c(connecteds2c == ivx) = [];

                    produit_de_tanh = 1;

                    for connected2c = connecteds2c
                        produit_de_tanh = produit_de_tanh * tanh(V2C(icy, connected2c)/2);
                    end
                    
                    C2V(icy, ivx) = 2*atanh(produit_de_tanh);
                end
            end
        end
        %% Boucle remplissage V2C
        for icy = 1:nb_noeudsparite
            for ivx = 1:nb_noeudsvar
                if(H(icy, ivx) == 1)
                    % Il y a un lien entre 2 noeuds ici
                    connecteds2v = find(H(:, ivx)); % on repere les liens mon noeud de var --> noeuds de parites

                    % On ne retient pas les liens sur lequel on est
                    connecteds2v(connecteds2v == icy) = [];

                    % LLR du bit du canal qui est arrive par ce noeud de var
                    somme = mdc(ivx);
                    for connected2v = connecteds2v
                        somme = somme + C2V(connected2v, ivx);
                    end
                end
            end
        end
    end
    LLRs = sum(V2C, 1); % Je somme sur les colonnes
    msg = LLRs + mdc;
    msg = msg(end-nb_noeudsparite+1:end)'; % car code systématique
    
    y((1:nb_noeudsparite)+nb_noeudsparite*(i_motdecode-1)) = msg(:);
    
end


