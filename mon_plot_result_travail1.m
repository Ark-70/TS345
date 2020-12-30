%% clean
% close all
% clc
% clear all

%% Files

addpath('rapport');
figure(1);
code_name = "DEBUG_6_3";
bridge = "_nb_ite_";
tab_ite = 1:5;
bornes_SNR = [0, 10]

%% Automated preparation of graph loop

nb_courbes = length(tab_ite);

% These are the default colors in matlab
C = [0, 0.4470, 0.7410;
0.8500, 0.3250, 0.0980;
0.4940, 0.1840, 0.5560;
0.4660, 0.6740, 0.1880;
0.3010, 0.7450, 0.9330;
0.6350, 0.0780, 0.1840]

tab = strings(1, nb_courbes)
for i=1:nb_courbes
    tab(i) = strcat(code_name, bridge, num2str(tab_ite(i)))
end
tab

%% Loop plot graph
for i=1:length(tab)
    
    load(tab(i))
    mabite = semilogy(EbN0dB,ber, 'color', C(i,:),'LineWidth', 1);
    
    if(i==1)
        hold on; 
    end
    semilogy(EbN0dB,per, '--', 'color', C(i,:),'LineWidth', 1);
    
end

%% Setup graph

xlim([bornes_SNR(1) bornes_SNR(2)])
ylim([1e-6 2])
grid on
xlabel('$\frac{E_b}{N_0}$ en dB','Interpreter', 'latex', 'FontSize',14)
ylabel('TEB','Interpreter', 'latex', 'FontSize',14)
    
legs = {'TEB 1 it�ration', 'TEP 1 it�ration', 'TEB 2 it�rations', 'TEP 2 it�rations', 'TEB 3 it�rations', 'TEP 3 it�rations', 'TEB 4 it�rations', 'TEP 4 it�rations', 'TEB 5 it�rations', 'TEP 5 it�rations'} 
legend(legs, 'Interpreter', 'latex', 'FontSize',8);
%legend({'TEB (non cod\''e)'}, 'Interpreter', 'latex', 'FontSize',14);