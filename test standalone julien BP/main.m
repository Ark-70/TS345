clear all
clc

% FAIRE make DANS LA CONSOLE POUR CREER LA FONCTION LDPC_H2G

%% Constantes USER

ENCODE_DECODE_ON = 1;
RANDOM_ON = 1;
NOISE_ON = 1;
SAVE_IN_FILE_ON = 1;
CHOIX_DU_CODE = 1; % entre 1 et 3
nb_iterations = 1;
USE_JULIEN_BP = 1;

%% Constantes  GRAPH

EbN0dB_min  = 0; % Minimum de EbN0
EbN0dB_max  = 10; % Maximum de EbN0
EbN0dB_step = 1;% Pas de EbN0

nbr_err_min  = 50; % Nombre d'erreurs a observer avant de calculer un BER
nbr_err_min_mauvais_snr = 300; % Nombre d'erreurs a observer sur un mauvais snr (dï¿½fini en dessous)
seuil_nbr_bits_requis = 1e5; % Si on ne dï¿½passe pas ce nombre de bits, on doit au moins avoir 200 erreurs
nbr_bit_max = 1e7; % Nombre de bits max a simuler
ber_min     = 1e-6; % BER min

%% Parametres calcules ou figes pour le TP
% -------------------------------------------------------------------------
addpath('src');
codes_path = ["DEBUG_6_3", "CCSDS_64_128", "MACKAY_504_1008"];
save_name = strcat(codes_path(CHOIX_DU_CODE), "_nb_ite_", int2str(nb_iterations));
construct_full_path = strcat('alist/', codes_path(CHOIX_DU_CODE), '.alist')
H = alist2sparse(construct_full_path); % on lit H dans un json-like version matlab
[h, g] = ldpc_h2g(H); % donne h et g systematiques (necessite le compilateur C de matlab add-on MinGW)
[height, width] = size(h);

% load('code2.mat')

R = 1-(rank(full(h))/width); % rendement de la communication

% Pour ce TP, on prend 1_msg = 1_paquet
pqt_par_trame = 1; % Nombre de paquets par trame
bit_par_msg = height;
bit_par_pqt = bit_par_msg;% Nombre de bits par paquet

K = pqt_par_trame*bit_par_pqt; % Nombre de bits de message par trame
N = K/R; % Nombre de bits codes par trame (codï¿½e)

M = 2; % Modulation BPSK <=> 2 symboles
phi0 = 0; % Offset de phase our la BPSK

EbN0dB = EbN0dB_min:EbN0dB_step:EbN0dB_max;     % Points de EbN0 en dB ï¿½ simuler
EbN0   = 10.^(EbN0dB/10);% Points de EbN0 ï¿½ simuler
EsN0   = R*log2(M)*EbN0; % Points de EsN0
EsN0dB = 10*log10(EsN0); % Points de EsN0 en dB ï¿½ simuler

% -------------------------------------------------------------------------

%% Construction du modulateur
mod_psk = comm.PSKModulator(...
    'ModulationOrder', M, ... % BPSK
    'PhaseOffset'    , phi0, ...
    'SymbolMapping'  , 'Gray',...
    'BitInput'       , true);

%% Construction du demodulateur
demod_psk = comm.PSKDemodulator(...
    'ModulationOrder', M, ...
    'PhaseOffset'    , phi0, ...
    'SymbolMapping'  , 'Gray',...
    'BitOutput'      , true,...
    'DecisionMethod' , 'Log-likelihood ratio');

%% Construction du canal AWGN
awgn_channel = comm.AWGNChannel(...
    'NoiseMethod', 'Signal to noise ratio (Es/No)',...
    'EsNo',EsN0dB(1),...
    'SignalPower',1);

%% Construction de l'objet evaluant le TEB
stat_erreur = comm.ErrorRate(); % Calcul du nombre d'erreur et du BER

%% Initialisation des vecteurs de resultats
ber = zeros(1,length(EbN0dB));
per = zeros(1,length(EbN0dB));
err_paquet = zeros(1,length(EbN0dB));
Pe = qfunc(sqrt(2*EbN0));

%% Preparation de l'affichage
figure(1)
h_ber = semilogy(EbN0dB,ber,'XDataSource','EbN0dB', 'YDataSource','ber');
hold all
ylim([1e-6 1])
grid on
xlabel('$\frac{E_b}{N_0}$ en dB','Interpreter', 'latex', 'FontSize',14)
ylabel('TEB','Interpreter', 'latex', 'FontSize',14)

%% Prï¿½paration de l'affichage en console
msg_format = '|   %7.2f  |   %9d   |  %9d | %2.2e |  %8.2f kO/s |   %8.2f kO/s |   %8.2f s |\n';

fprintf(      '|------------|---------------|------------|----------|----------------|-----------------|--------------|\n')
msg_header =  '|  Eb/N0 dB  |    Bit nbr    |  Bit err   |   TEB    |    Debit Tx    |     Debit Rx    | Tps restant  |\n';
fprintf(msg_header);
fprintf(      '|------------|---------------|------------|----------|----------------|-----------------|--------------|\n')


%% Simulation
for i_snr = 1:length(EbN0dB)
    reverseStr = ''; % Pour affichage en console
    awgn_channel.EsNo = EsN0dB(i_snr);% Mise a jour du EbN0 pour le canal

    stat_erreur.reset; % reset du compteur d'erreur
    err_stat    = [0 0 0]; % vecteur rï¿½sultat de stat_erreur
    nb_err_paquets = 0;
    taux_err_paquets = 0;

    demod_psk.Variance = awgn_channel.Variance;

    n_frame = 0;
    T_rx = 0;
    T_tx = 0;
    general_tic = tic;


%     nb_err_bits_verite_terrain = 0

    while ((err_stat(3)<seuil_nbr_bits_requis || err_stat(2)<nbr_err_min) && err_stat(2) < nbr_err_min_mauvais_snr && err_stat(3) < nbr_bit_max)
        n_frame = n_frame + 1;

        %% Emetteur
        tx_tic = tic;                 % Mesure du dï¿½bit d'encodage

        if(RANDOM_ON)
            b    = randi([0,1],K,1);    % Gï¿½nï¿½ration du message alï¿½atoire
        else
            b    = ones(K,1);    % Gï¿½nï¿½ration du message alï¿½atoire
        end

        if(ENCODE_DECODE_ON)
            code = encoder(g, b); % encodage LDPC
        else 
            code = b;
        end
            
        x      = step(mod_psk,  code); % Modulation BPSK


        T_tx   = T_tx+toc(tx_tic);    % Mesure du dï¿½bit d'encodage

        %% Canal

        if(NOISE_ON)
            y     = step(awgn_channel,x); % Ajout d'un bruit gaussien
        else
            y = x;
        end

        %% Recepteur
        rx_tic = tic;                  % Mesure du dï¿½bit de dï¿½codage
        Lch      = step(demod_psk,y);   % Dï¿½modulation (retourne des LLRs)

        if(ENCODE_DECODE_ON)
            if(USE_JULIEN_BP)
                y = BeliefProp_JULIEN(h, Lch, nb_iterations);
            else
                y = decoder(Lch, h, nb_iterations); % CEST PAS DU JULIEN    
            end
        else
            y = Lch;
        end

        if(USE_JULIEN_BP)
            rec_b = y';
        else
            rec_b = double(y(1:K) < 0); % Dï¿½cision CEST PAS DU JULIEN    
        end


        T_rx    = T_rx + toc(rx_tic);  % Mesure du dï¿½bit de dï¿½codage

        err_stat   = step(stat_erreur, b, rec_b); % Comptage des erreurs binaires et des taux

        nb_err_paquets = nb_err_paquets + any(rec_b~=b); % Comptage des erreurs paquets//motdecode
        taux_err_paquets = nb_err_paquets /n_frame;

        %% Affichage du rï¿½sultat
        if mod(n_frame,100) == 1
            msg = sprintf(msg_format,...
                EbN0dB(i_snr),         ... % EbN0 en dB
                err_stat(3),           ... % Nombre de bits envoyï¿½s
                err_stat(2),           ... % Nombre d'erreurs observï¿½es
                err_stat(1),           ... % BER
                err_stat(3)/8/T_tx/1e3,... % Dï¿½bit d'encodage
                err_stat(3)/8/T_rx/1e3,... % Dï¿½bit de dï¿½codage
                toc(general_tic)*(nbr_err_min - min(err_stat(2),nbr_err_min))/(min(err_stat(2),nbr_err_min))); % Temps restant
            fprintf(reverseStr);
            msg_sz =  fprintf(msg);
            reverseStr = repmat(sprintf('\b'), 1, msg_sz);
        end

    end

    msg = sprintf(msg_format,...
        EbN0dB(i_snr),         ... % EbN0 en dB
        err_stat(3),           ... % Nombre de bits envoyï¿½s
        err_stat(2),           ... % Nombre d'erreurs observï¿½es
        err_stat(1),           ... % BER
        err_stat(3)/8/T_tx/1e3,... % Dï¿½bit d'encodage
        err_stat(3)/8/T_rx/1e3,... % Dï¿½bit de dï¿½codage
        0); % Temps restant
    fprintf(reverseStr);
    msg_sz =  fprintf(msg);
    reverseStr = repmat(sprintf('\b'), 1, msg_sz);

    ber(i_snr) = err_stat(1);
    per(i_snr) = nb_err_paquets/(err_stat(3)/height);
    refreshdata(h_ber);
    drawnow limitrate

    if err_stat(1) < ber_min
        break
    end

end
fprintf('|------------|---------------|------------|----------|----------------|-----------------|--------------|\n')

%%

if(SAVE_IN_FILE_ON && ENCODE_DECODE_ON && RANDOM_ON && NOISE_ON)
    path = get_free_path(save_name);
    save(path,'EbN0dB','ber', 'per');
    disp(strcat('saved on : ',path));
else
    warning("data is not saved because one of the User Constants (ON/OFF variables) was not appropriate. Save workspace now with save('datafilename.mat')"); 
    warning("Les données n'ont pas été sauvegardées parce qu'une des constantes utilisateurs (variables ON/OFF) semble non appropriée. Sauvegardez maintenant toutes vos variables avec save('datafilename.mat')"); 
end

figure(1)
semilogy(EbN0dB,ber);
hold all;
semilogy(EbN0dB,per);
xlim([0 10])
ylim([1e-6 1])
grid on
xlabel('$\frac{E_b}{N_0}$ en dB','Interpreter', 'latex', 'FontSize',14)
ylabel('TEB','Interpreter', 'latex', 'FontSize',14)


