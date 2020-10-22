[  ] = propagate_decod( H, Lc )


%%

[H] = alist2sparse('alist/DEBUG_6_3.alist');
[h, g] = ldpc_h2g(H) % on a 2 H, un H systématique et un H de base designé
% Les deux réponses sont les systématiques
% On utilise la G systématique et la H de base

matrix_C2V = sparse(size(H));
matrix_V2C = matrix_C2V;

%% Initialisation

matrix_C2V

end  % function
