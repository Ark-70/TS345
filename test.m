clc, clear all, close all

EbN0dB = -2:1:3;     % Points de EbN0 en dB � simuler
truc = [1 10 100 1000 10000 100000]
% figure(1);
semilogy(EbN0dB, truc);

hold all;

semilogy(EbN0dB, truc*2);
