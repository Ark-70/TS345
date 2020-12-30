function newpath = get_free_path(save_name)


    if(exist( strcat(save_name, '.mat'), 'file') ~= 2) % Si c'est dispo sans les numéros i en plus
        newpath = save_name;
        
    else
        
        i = 2;
        newpath = strcat(save_name, '__', int2str(i), '.mat');
        
        while(exist(newpath, 'file') == 2) % Si c'est toujours pas dispo
            i = i + 1;
            newpath = strcat(save_name,'__',int2str(i), '.mat');
        end
    end
    
    newpath; % return
    
end

