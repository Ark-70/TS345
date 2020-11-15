%% clean
close all
clc
clear all
%%

addpath('rapport');
figure(1);
load("1ite.mat");
semilogy(EbN0dB,ber);
hold all
semilogy(EbN0dB,per,'--');
xlim([-2 10])
ylim([1e-6 2])
grid on
xlabel('$\frac{E_b}{N_0}$ en dB','Interpreter', 'latex', 'FontSize',14)
ylabel('TEB','Interpreter', 'latex', 'FontSize',14)


tab = ["2ite.mat", "3ite.mat", "4ite.mat", "5ite.mat"];
for val = tab
    load(val)
    semilogy(EbN0dB,ber);
    semilogy(EbN0dB,per, '--');
end
    
legs = {'TEB 1 it�ration', 'TEP 1 it�ration', 'TEB 2 it�rations', 'TEP 2 it�rations', 'TEB 3 it�rations', 'TEP 3 it�rations', 'TEB 4 it�rations', 'TEP 4 it�rations', 'TEB 5 it�rations', 'TEP 5 it�rations'} 
legend(legs, 'Interpreter', 'latex', 'FontSize',8);
%legend({'TEB (non cod\''e)'}, 'Interpreter', 'latex', 'FontSize',14);