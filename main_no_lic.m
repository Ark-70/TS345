clear
clc

%% Parametres
% -------------------------------------------------------------------------
addpath('src')
simulation_name = 'non_codee';

R = 1; % Rendement de la communication

pqt_par_trame = 1; % Nombre de paquets par trame
bit_par_pqt   = 330;% Nombre de bits par paquet
K = pqt_par_trame*bit_par_pqt; % Nombre de bits de message par trame
N = K/R; % Nombre de bits codes par trame (codee)

M = 2; % Modulation BPSK <=> 2 symboles
phi0 = 0; % Offset de phase our la BPSK

EbN0dB_min  = 10; % Minimum de EbN0
EbN0dB_max  = 10; % Maximum de EbN0
EbN0dB_step = 1;% Pas de EbN0

nbr_erreur  = 100;  % Nombre d'erreurs e observer avant de calculer un BER
nbr_bit_max = 100e6;% Nombre de bits max e simuler
ber_min     = 1e-6; % BER min

EbN0dB = EbN0dB_min:EbN0dB_step:EbN0dB_max;     % Points de EbN0 en dB e simuler
EbN0   = 10.^(EbN0dB/10);% Points de EbN0 a simuler
EsN0   = R*log2(M)*EbN0; % Points de EsN0
EsN0dB = 10*log10(EsN0); % Points de EsN0 en dB a simuler

% -------------------------------------------------------------------------
%% Initialisation des vecteurs de resultats
ber = zeros(1,length(EbN0dB));
Pe = 0.5*erfc(sqrt(EbN0));

%% Preparation de l'affichage
figure(1)
h_ber = semilogy(EbN0dB,ber,'XDataSource','EbN0dB', 'YDataSource','ber');
hold all
ylim([1e-6 1])
grid on
xlabel('$\frac{E_b}{N_0}$ en dB','Interpreter', 'latex', 'FontSize',14)
ylabel('TEB','Interpreter', 'latex', 'FontSize',14)

%% Preparation de l'affichage en console
msg_format = '|   %7.2f  |   %9d   |  %9d | %2.2e |  %8.2f kO/s |   %8.2f kO/s |   %8.2f s |\n';

fprintf(      '|------------|---------------|------------|----------|----------------|-----------------|--------------|\n')
msg_header =  '|  Eb/N0 dB  |    Bit nbr    |  Bit err   |   TEB    |    Debit Tx    |     Debit Rx    | Tps restant  |\n';
fprintf(msg_header);
fprintf(      '|------------|---------------|------------|----------|----------------|-----------------|--------------|\n')

[H] = alist2sparse('alist/DEBUG_6_3.alist');
g = ldpc_h2g(H); % on a 2 H, un h systématique et un H de base bien design é

%% Simulation
for i_snr = 1:length(EbN0dB)
    reverseStr = ''; % Pour affichage en console
    sigma2 = 1/(2*EsN0(i_snr));

    err_stat    = [0 0 0]; % vecteur resultat de stat_erreur

    n_frame = 0;
    T_rx = 0;
    T_tx = 0;
    general_tic = tic;
    while (err_stat(2) < nbr_erreur && err_stat(3) < nbr_bit_max)
        n_frame = n_frame + 1;

        %% Mise en commun Emetteur/Recepteur
        % On utilise la G systématique pour générer et la H de base pour décoder


        %% Emetteur
        tx_tic = tic;                 % Mesure du debit d'encodage
        b      = randi([0,1],K,1);    % Generation du message aleatoire
%         b = zeros(K,1);

        code = encoder(g, b); % encodage LDPC

        x      = 1 - 2*code; % Modulation BPSK
        T_tx   = T_tx+toc(tx_tic);    % Mesure du debit d'encodage

        %% Canal
        n     = sqrt(sigma2) * randn(size(x));
        
        y     = x + n; % Ajout d'un bruit gaussien

        %% Recepteur
        rx_tic = tic;                  % Mesure du debit de decodage
        Lch      = 2*y/sigma2;         % Demodulation (retourne des LLRs)

        y = decoder(Lch, H);

%         rec_b = double(y) < 0);   % Decision
        rec_b = double(y(1:K) < 0);   % Decision
        T_rx    = T_rx + toc(rx_tic);  % Mesure du debit de decodage

        err_stat(2) = err_stat(2) + sum(b(:) ~= rec_b(:));
        err_stat(3) = err_stat(3) + K;
        err_stat(1) = err_stat(2)/err_stat(3);

        %% Affichage du resultat
        if mod(n_frame,1) == 0
            msg = sprintf(msg_format,...
                EbN0dB(i_snr),         ... % EbN0 en dB
                err_stat(3),           ... % Nombre de bits envoyes
                err_stat(2),           ... % Nombre d'erreurs observees
                err_stat(1),           ... % BER
                err_stat(3)/8/T_tx/1e3,... % Debit d'encodage
                err_stat(3)/8/T_rx/1e3,... % Debit de decodage
                toc(general_tic)*(nbr_erreur - min(err_stat(2),nbr_erreur))/(min(err_stat(2),nbr_erreur))); % Temps restant
            fprintf(reverseStr);
            msg_sz =  fprintf(msg);
            reverseStr = repmat(sprintf('\b'), 1, msg_sz);
        end

    end

    msg = sprintf(msg_format,...
        EbN0dB(i_snr),         ... % EbN0 en dB
        err_stat(3),           ... % Nombre de bits envoyes
        err_stat(2),           ... % Nombre d'erreurs observees
        err_stat(1),           ... % BER
        err_stat(3)/8/T_tx/1e3,... % Debit d'encodage
        err_stat(3)/8/T_rx/1e3,... % Debit de decodage
        0); % Temps restant
    fprintf(reverseStr);
    msg_sz =  fprintf(msg);
    reverseStr = repmat(sprintf('\b'), 1, msg_sz);

    ber(i_snr) = err_stat(1);
    refreshdata(h_ber);
    drawnow limitrate

    if err_stat(1) < ber_min
        break
    end

end
fprintf('|------------|---------------|------------|----------|----------------|-----------------|--------------|\n')

%%
figure(1)
semilogy(EbN0dB,ber);
hold all
xlim([0 10])
ylim([1e-6 1])
grid on
xlabel('$\frac{E_b}{N_0}$ en dB','Interpreter', 'latex', 'FontSize',14)
ylabel('TEB','Interpreter', 'latex', 'FontSize',14)

save(simulation_name,'EbN0dB','ber')
