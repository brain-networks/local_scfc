function nav = navigate(C, MS, max_hops)

% Function to perform navigation routing in a connecitivty matrix C, guided by
% the metric space MS.
% Inputs:
% C - square connectivity matrix
% MS - square symmetric matrix of distances between nodes in the metric
% space
% max_hops (optinal) - limits the maximum number of hops of navigation
% paths
% Output:
% nav - structure containing:
% .inputs - User specific inputs
% .paths - cell matrix of navigation paths
% .num_hops - number of hops of navigation paths
% .pl_MS - navigation path lengths in metric space
% .failed_paths - binary matrix of failed paths (failure = 1)
%
% Caio Seguin, University of Melbourne, 2017

    if nargin == 2
        max_hops = length(C);
    end

    components = get_components(C);

    if (length(unique(components)) > 1)
        fprintf('Warning: disconnected connectivity matrix\n');
    end

    n = size(C, 1);
    
    nav.inputs.C = C;
    nav.inputs.MS = MS;
    nav.inputs.max_hops = max_hops;
    
    nav.paths = cell(n);
    nav.num_hops = zeros(n);
    nav.pl_MS = zeros(n);
    nav.failed_paths = zeros(n); nav.failed_paths(1:n+1:end) = 1;
    
    for i = 1:n
        for j = 1:n
            if (i ~= j)
                
                curr_node = i;
                last_node = curr_node;
                target = j;
                num_hops = 0;
                pl_MS = 0;
                nav.paths{i,j} = curr_node;
                
                while (curr_node ~= target)
                    
                    neighbors = find(C(curr_node,:) ~= 0);
                    [~, min_index] = min(MS(target, neighbors));
                    next_node = neighbors(min_index);
                    
                    if isempty(next_node) || next_node == last_node || num_hops > max_hops

                        num_hops = 0;
                        pl_MS = 0;
                        nav.failed_paths(i,j) = 1;
                        
                        break
                    
                    end
                    
                    nav.paths{i,j} = [nav.paths{i,j} next_node];
                    num_hops = num_hops + 1;
                    pl_MS = MS(curr_node, next_node) + pl_MS;
                    
                    last_node = curr_node;
                    curr_node = next_node;
                
                end

                nav.num_hops(i,j) = num_hops;
                nav.pl_MS(i,j) = pl_MS;
                
            end
        end
    end

end