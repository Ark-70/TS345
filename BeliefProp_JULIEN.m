function b = BeliefProp_JULIEN(H, Lc)
%Belief propagation (LLR en entrée,bits en sortie)

iterations = 4;

[m,n] = size(H);


% Lc de taille 1xn;


[y,x] = find(H); %indices des valeurs non nulles
pos = find(H); %donne la liste à parcourir
%full(H(7)); donne la valeur à la position de la liste
%find(H(1,:)) %donne les indices de H a sommer

v2c = sparse(m,n);
c2v = v2c;

%LLRtemp = ones(length(x));

%% INITIALISATION


% remplir v2c avec LLR initiaux
for i = 1:length(x)
    v2c(pos(i)) = Lc(x(i));
end

v2cInit = v2c;

% itérer sur c2v une fois


for i = 1:length(x)
    indices = logical(H(y(i),:));
    indices(x(i)) = false;
    %indices = find(H(y(i),:)); %[1 5] indices en x
    %indices(indices == x(i)) = [];
    %valeurs = full(v2c(y(i),indices));
    valeurs = v2c(y(i),indices);
    %LLRtemp = 2*atanh(prod(tanh(valeurs/2)));
    LLRtemp = prod(tanh(valeurs/2));
    c2v(pos(i)) = LLRtemp;

end
    c2v = 2*atanh(c2v);


%% BOUCLE PRINCIPALE

for cpt=1:iterations
    
    %% FILL V2C

    for i = 1:length(x)
        indices = logical(H(:,x(i)));
        indices(y(i)) = false;
        %indices = find(H(:,x(i))); %[1 3] indices en y
        %indices(indices == y(i)) = [];
        %valeurs = full(c2v(indices,x(i)));
        valeurs = c2v(indices,x(i));
        LLRtemp = sum(valeurs);
        v2c(pos(i)) = LLRtemp;

    end
    v2c = v2c + v2cInit;
    


    %% FILL C2V


    for j = 1:length(x)
        indices = logical(H(y(j),:));
        indices(x(j)) = false;
        %indices = find(H(y(j),:)); %[1 5] indices en x
        %indices(indices == x(j)) = [];
        %valeurs = full(v2c(y(j),indices));
        valeurs = v2c(y(j),indices);
        %LLRtemp = 2*atanh(prod(tanh(valeurs/2)));
        LLRtemp = prod(tanh(valeurs/2));
        c2v(pos(j)) = LLRtemp;
    end
    c2v = 2*atanh(c2v);


end

%% DECISION

b_temp = full(sum(c2v)) + Lc.';
b = double(b_temp(:,n-m+1:n)<0);





end






