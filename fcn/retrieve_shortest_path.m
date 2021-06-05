
function [path] = retrieve_shortest_path(s,d,SPL,B)

path_length = SPL(s,d);
path = nan(path_length+1,1);
path(1) = s;
for ind = 2:length(path);
    s = B(s,d);
    path(ind) = s;    
end




