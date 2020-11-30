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
[all_icy, all_ivx] = find(H); % ne changent jamais
nb_iterations = 5; % pour chaque mot de code

for i_motdecode = 1:size(Lch_matrix, 1)

    mdc = Lch_matrix(i_motdecode, :);

    LLR_repetes = repmat(mdc, 3, 1);

    %% Initialisation V2C/C2V

    C2V = zeros(size(H)); % Observation comme des parites dans les noeuds de var
    V2C = H.*LLR_repetes; % Noeuds de var dans les parites
    % [y, x] = find(H);


    for n = 1:nb_iterations
        %% Boucle remplissage C2V
        n;
        full(C2V);
        full(V2C);
        nb_iterations;
%         for icy = cy % Ca marche pas du tout cette merde, monsieur
        for i_c = 1:length(all_icy)
            icy = all_icy(i_c); % i_c c'est un index pour un tableau d'index #MindBlown #JeanReviensPas
            ivx = all_ivx(i_c);
            % Il y a un lien entre 2 noeuds ici

            % On trouve ou sont les liens relies a ce noeud
            connecteds2c = find(H(icy, :)); % on repere les liens mon noeud de parite --> noeuds de var

            % On ne retient pas les liens sur lequel on est
            connecteds2c(connecteds2c == ivx) = [];

            C2V(icy, ivx) = 2*atanh(prod(tanh(connecteds2c/2)));
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
            for connected2v = connecteds2v
                somme = somme + C2V(connected2v, ivx);
            end
            V2C(icy, ivx) = somme;
        end

    end
    LLRs = sum(V2C, 1); % Je somme sur les colonnes
    msg = LLRs + mdc;
    msg = msg(end-nb_noeudsparite+1:end)'; % car code syst�matique

    y((1:nb_noeudsparite)+nb_noeudsparite*(i_motdecode-1)) = msg(:);

end
