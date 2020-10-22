function decoder( Lc, H )


%% Initialisation

matrix_C2V = sparse(size(H)); % matrice de 0
matrix_V2C = matrix_C2V;

[y, x] = find(H);

%% première passe avec tout C2V à 0 donc ca laisse passer que les LLRs

for i=1:length(x)
%     matrix_V2C(y, x) = 


end  % function
