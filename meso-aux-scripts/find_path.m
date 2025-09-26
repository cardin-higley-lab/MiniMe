function found_path = find_path(subfolder)
%find_path find path out of known places we store data
%   Searches known data directories (stored in possible_input_paths.txt) 
%   for the one that kinds the subfolder, which will be a unique key, and 
%   return the corresponding directory

% read possible inputs
    possible_paths = readlines("possible_input_paths.txt");
    found_path = "";
    for p = 1:length(possible_paths)
        % test inputs with subfolder
        if isfolder(fullfile(possible_paths(p),subfolder))
            % return the corresponding correct input
            found_path = possible_paths(p);
            break
        end
    end
end