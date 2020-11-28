clear
clc

%% Parametres
% -------------------------------------------------------------------------
addpath('src')
simulation_name = 'non_codee';

R = 1; % Rendement de la communication

n_trame = 110; % Nombre de paquets par trame
k = 3;% Nombre de bits par paquet
K = n_trame*k; % Nombre de bits de message par trame
N = K/R; % Nombre de bits cod�s par trame (cod�e)

M = 2; % Modulation BPSK <=> 2 symboles
phi0 = 0; % Offset de phase our la BPSK

EbN0dB_min  = 0; % Minimum de EbN0
EbN0dB_max  = 5; % Maximum de EbN0
EbN0dB_step = 1;% Pas de EbN0

nbr_erreur  = 100;  % Nombre d'erreurs � observer avant de calculer un BER
nbr_bit_max = 100e6;% Nombre de bits max � simuler
ber_min     = 1e-6; % BER min

EbN0dB = EbN0dB_min:EbN0dB_step:EbN0dB_max;     % Points de EbN0 en dB � simuler
EbN0   = 10.^(EbN0dB/10);% Points de EbN0 � simuler
EsN0   = R*log2(M)*EbN0; % Points de EsN0
EsN0dB = 10*log10(EsN0); % Points de EsN0 en dB � simuler

% -------------------------------------------------------------------------
%% Initialisation des vecteurs de r�sultats
ber = zeros(1,length(EbN0dB));
Pe = 0.5*erfc(sqrt(EbN0));
err_paquet = zeros(1,length(EbN0dB));


%% Pr�paration de l'affichage
figure(1)
h_ber = semilogy(EbN0dB,ber,'XDataSource','EbN0dB', 'YDataSource','ber');
hold all
ylim([1e-6 1])
grid on
xlabel('$\frac{E_b}{N_0}$ en dB','Interpreter', 'latex', 'FontSize',14)
ylabel('TEB','Interpreter', 'latex', 'FontSize',14)

%% Pr�paration de l'affichage en console
msg_format = '|   %7.2f  |   %9d   |  %9d | %2.2e |  %8.2f kO/s |   %8.2f kO/s |   %8.2f s |\n';

fprintf(      '|------------|---------------|------------|----------|----------------|-----------------|--------------|\n')
msg_header =  '|  Eb/N0 dB  |    Bit nbr    |  Bit err   |   TEB    |    Debit Tx    |     Debit Rx    | Tps restant  |\n';
fprintf(msg_header);
fprintf(      '|------------|---------------|------------|----------|----------------|-----------------|--------------|\n')


[H] = alist2sparse('alist/DEBUG_6_3.alist');
 %[h,g] = ldpc_h2g(H);
 g = [1 1 1 1 0 0; 0 1 0 0 1 0; 1 0 1 0 0 1];


%% Simulation
for i_snr = 1:length(EbN0dB)

    reverseStr = ''; % Pour affichage en console
    sigma2 = 1/(2*EsN0(i_snr));
    
    err_stat    = [0 0 0]; % vecteur r�sultat de stat_erreur
    nb_err_paquet = 0;
    
    n_frame = 0;
    T_rx = 0;
    T_tx = 0;
    general_tic = tic;
    
    
    
    while (err_stat(2) < nbr_erreur && err_stat(3) < nbr_bit_max)
        n_frame = n_frame + 1;
        
        %% Emetteur
        tx_tic = tic;                 % Mesure du d�bit d'encodage
        b      = randi([0,1],n_trame,k);    % G�n�ration du message al�atoire
        %b = zeros(K,1);
        
        %ON RAJOUTE UN ENCODEUR ICI
         
         c = mod(b * g,2);
         % x_code = mod(reshape(x_cod.', 660,1),2);
         x=1-2*c;
         
         T_tx   = T_tx+toc(tx_tic);    % Mesure du d�bit d'encodage
        
        %% Canal
        
        n     = sqrt(sigma2) * randn(size(x));
        y     = x + n; % Ajout d'un bruit gaussien
        
        %% Recepteur
        rx_tic = tic;                  % Mesure du d�bit de d�codage
        Lc      = 2*y/sigma2;   % D�modulation (retourne des LLRs)
        
        %% DECODEUR ICI AVEC Définition de v_proc et c_proc
% v_proc (resp. c_proc) prend en argument une matrice de contrainte 
%vers variable (resp. variable vers contrainte), et renvoie une matrice 
% de variable vers contrainte (resp. contrainte vers variable). On fait
% ainsi évoluer le décodeur BP par itérations successives de v_proc et
% c_proc.

        decision_finale = zeros(n_trame, k);
        for n_mess=1:n_trame
            c2v = zeros(size(H));
            v2c = zeros(size(H));

            
            for z=1:1%Nombre d'it�rations pour algo propagation de croyance
                for i=1:6 %ON CALCULE V2C 
                    connected2vi = find(H(:,i));
                    
                    for j=1:length(connected2vi)% calcul v_i -> c_j
                        cj = connected2vi(j);
                        connected = connected2vi(connected2vi ~= cj); 
                        v2c(cj,i) = Lc(n_mess,i);
                        for c = connected % loop c -> v_i, c != c_j 
                            v2c(cj,i) = v2c(cj,i) + c2v(c,i);
                            if isnan(v2c(cj,i))
                                v2c(cj,i)=0;
                            end
                        end
                    end
                end

                for i=1:3 %ON CALCULE C2V
                    %On cherche les noeuds de variables reliés aux noeuds de
                    %contraintes i
                    connected2ci = find(H(i,:));
                    
                    for j=1:length(connected2ci)% calcul c_i -> v_j
                        vj = connected2ci(j);
                        connected = connected2ci(connected2ci ~= vj);
                        produit = 1;
                        for c = connected
                            produit = produit * tanh(v2c(i,c)/2);
                        end
                        c2v(i,vj) = 2*atanh(produit);
                    end
                end

            end %A ce stade on a obtenu c2v et v2c définitifs.
            %Décodage de 1 message (composé de 2 messages de 3 bits)

            %Pour décoder le message, on somme de c2v + canal pour
            %et on regarde le signe
            decision = Lc(n_mess,:) + sum(c2v);

            decision_finale(n_mess,:) = decision(:,4:6)<0;
            %Il faut g�rer l'erreur paquet a ce niveau
            if decision_finale(n_mess,:)~=b(n_mess,:)
                nb_err_paquet = nb_err_paquet+1;
            end
        end
        
        %% Décision + erreur 
        
        %rec_b = double(Lc(1:K) < 0); % D�cision
        rec_b = decision_finale;
        %rec_b = reshape(rec_b, 330,1);
  
        T_rx    = T_rx + toc(rx_tic);  % Mesure du d�bit de d�codage
        
        err_stat(2) = err_stat(2) + sum(b(:) ~= rec_b(:));
        err_stat(3) = err_stat(3) + K;
        err_stat(1) = err_stat(2)/err_stat(3);
        
        %% Affichage du r�sultat
        if mod(n_frame,100) == 1
             msg = sprintf(msg_format,...
                EbN0dB(i_snr),         ... % EbN0 en dB
                err_stat(3),           ... % Nombre de bits envoy�s
                err_stat(2),           ... % Nombre d'erreurs observ�es
                err_stat(1),           ... % BER
                err_stat(3)/8/T_tx/1e3,... % D�bit d'encodage
                err_stat(3)/8/T_rx/1e3,... % D�bit de d�codage
                toc(general_tic)*(nbr_erreur - min(err_stat(2),nbr_erreur))/(min(err_stat(2),nbr_erreur))); % Temps restant
            fprintf(reverseStr);
            msg_sz =  fprintf(msg);
            reverseStr = repmat(sprintf('\b'), 1, msg_sz);
        end
        
    end
    
%     msg = sprintf(msg_format,...
%         EbN0dB(i_snr),         ... % EbN0 en dB
%         err_stat(3),           ... % Nombre de bits envoy�s
%         err_stat(2),           ... % Nombre d'erreurs observ�es
%         err_stat(1),           ... % BER
%         err_stat(3)/8/T_tx/1e3,... % D�bit d'encodage
%         err_stat(3)/8/T_rx/1e3,... % D�bit de d�codage
%         0); % Temps restant
%     fprintf(reverseStr);
%     msg_sz =  fprintf(msg);
%     reverseStr = repmat(sprintf('\b'), 1, msg_sz);
   
   ber(i_snr) = err_stat(1);
   err_paquet(i_snr) = nb_err_paquet/(err_stat(3)/K);
   
   refreshdata(h_ber);
   drawnow limitrate
    
    if err_stat(1) < ber_min
        break
    end
    
end
%fprintf('|------------|---------------|------------|----------|----------------|-----------------|--------------|\n')

%%
figure(1)
semilogy(EbN0dB,ber);
hold all
semilogy(EbN0dB,err_paquet);
xlim([0 10])
ylim([1e-6 1])
grid on
xlabel('$\frac{E_b}{N_0}$ en dB','Interpreter', 'latex', 'FontSize',14)
ylabel('TEB','Interpreter', 'latex', 'FontSize',14)

save(simulation_name,'EbN0dB','ber')




