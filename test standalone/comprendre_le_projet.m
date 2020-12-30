%% Setting up

CHOIX_DU_CODE = 3; % entre 1 et 3

addpath('src');
codes_path = ["DEBUG_6_3", "CCSDS_64_128", "MACKAY_504_1008"];

construct_full_path = strcat('alist/', codes_path(CHOIX_DU_CODE), '.alist')
H = alist2sparse(construct_full_path); % on lit H dans un json-like version matlab
[h, g] = ldpc_h2g(H); % donne h et g systematiques (necessite le compilateur C de matlab add-on MinGW)


%% Si je comprends bien et selon le code qui marche (etonnemment) de ce projet
    figure();

    % Ca c'est la matrice de parite, construite avec precaution et tout
    % bien designe par les chercheurs (on devine une procedure iterative pour la construire)
    subplot(3, 2, 1);
    spy(H); title('H'); 
    
    % Sauf qu'elle correspond pas (ou pas toujours) a une matrice génératrice systematique
    % Du coup on va devoir en faire un peu de la bouillie avec des permutations de 
    % colonnes pour qu'elle corresponde a une G systematique
    subplot(3, 2, 2);
    spy(h); title('h'); 
    
    % Par ailleurs on a notre G systematique
    truc = subplot(3, 2, [3, 4]);
    spy(g); title('g'); 
    
    % On voit bien que g*H ça marche pas bien, ca fait pas 0
    subplot(3, 2, 5);
    spy( mod(g*H', 2) ); title("g*H'"); 
    
    % En revanche g*h marche bien, ca fait bien 0 
    subplot(3, 2, 6);
    spy( mod(g*h', 2) ); title("g*h'"); 
    
    % la H donnée du code (6, 3) correspond deja a une matrice g systematique
    % donc ca ne changera rien pour ce code
    
    % BONUS : En cours nous on voyait comment passer d'une G systematique a
    % une H correspondante qui faisait apparaitre l'identite, mais ca changeait 
    %les dimensions de la matrice. Or, ici, dans les fonctions sombres, il doit 
    % y avoir une contrainte pour garder les memes dimensions peut-etre, 
    % ce qui fait surement qu'on voit pas H comme une partie Identite dedans