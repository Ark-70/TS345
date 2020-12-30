clear 
close all
clc


[H] = alist2sparse('alist/DEBUG_6_3.alist');
[h, g] = ldpc_h2g(H)