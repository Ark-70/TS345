function [final_LLRs] = C2VMatrix_to_final_LLRs(C2V, mdc, nb_noeudsparite)
% Quand on a notre matrice C2V qui a converg�, on peut en d�duire nos LLRs
% finaux, avant d'en d�duire les bits par seuil de d�cision

    LLRs = sum(C2V, 1); % Je somme sur les colonnes
    final_LLRs = LLRs + mdc;
    final_LLRs = final_LLRs(end-nb_noeudsparite+1:end)'; % car code systematique


end

