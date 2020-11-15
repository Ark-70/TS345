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
    
legs = {'TEB 1 itération', 'TEP 1 itération', 'TEB 2 itérations', 'TEP 2 itérations', 'TEB 3 itérations', 'TEP 3 itérations', 'TEB 4 itérations', 'TEP 4 itérations', 'TEB 5 itérations', 'TEP 5 itérations'} 
legend(legs, 'Interpreter', 'latex', 'FontSize',8);
%legend({'TEB (non cod\''e)'}, 'Interpreter', 'latex', 'FontSize',14);