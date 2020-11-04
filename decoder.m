function decoder( Lc, H )


%% Initialisation

matrix_C2V = sparse(size(H)); % matrice de 0
matrix_V2C = matrix_C2V;

[y, x] = find(H);

%% première passe avec tout C2V à 0 donc ca laisse passer que les LLRs
nb_noeudsparite = size(H, 2);
Lc_matrix = reshape(Lc, nb_noeudsparite, length(Lc)/nb_noeudsparite)';

%% Je veux mettre tous les noeuds de parités (carrés) dans les noeuds de var (ronds)
% sachant qu'on comprend tous les n_parités = 0, + sauf ceux qu'on rajoute
% qui sont les observations du canal

% Dans les 6 de Lc_matrix
% je dois refourguer Lc(1) dans V2C(1, 1) V2C(2, 1) V2C(3, 1) etc.,
% exceptés ceux qui sont à 0 dans H

% On a des msg de 3, codés avec redondance en mot de code de 6

% premier mot de code
premier_mdc = Lc_matrix(1, :);

% On veut dupliquer le premier bit du mot de code dans tous les noeuds de
% parité correspondant

% Après on fait des +(XOR) pour avoir la valeur des noeuds de variables
% find(H:, premier_mdc(1));

%% Init = V2C=LLR||0
% Itération = 1

LLR_repetes = repmat(premier_mdc, 3, 1);
matrix_V2C = H.*LLR_repetes;

matrix_C2V = 

for i = 1:
2*atanh*(tanh(LC1/2)*tanh(LC2/2)*tanh(LC3/2))
