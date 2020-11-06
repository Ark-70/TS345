function msg_vect = encoder(g, b)
%
%     test = 1:42
%     test2 = reshape(test, 6, 7)'
%     pour etre franc je suis pas sur de mon encodeur
%     reshape(test2', size(test2, 1)*size(test2, 2), 1)  
% 

    p_len = size(g, 1); % on veut des paquets de longueur 3 = n de g
    b = [b ; zeros(mod(length(b), p_len))]; % on s'assure d'avoir un vecteur de b divisible par la longueur de paquets
    paquets = reshape(b, p_len, []);
%     nb_repet_g = ceil(length(b)/size(g, 2));
%     g_fat = repmat(g, 1, nb_repet_g);
    msg = mod(paquets'*g, 2);
    msg_vect = reshape(msg', [], 1); % c'est dégueu les non-' mais jpp

end
