        
function [f] = nodal_areas(topology, coordinates)
    num_elements = size(topology, 1);
    num_nodes = size(coordinates, 1); 
    f = zeros(num_nodes, 1); 
    element_areas = zeros(num_elements, 1);
    
    % Calcolo dell'area per ogni elemento una volta sola
    for j = 1:num_elements
        nodes = topology(j, :);
        x = coordinates(nodes, 1);
        y = coordinates(nodes, 2);
        element_areas(j) = 0.5 * abs(det([1, x(1), y(1); 1, x(2), y(2); 1, x(3), y(3)]));
    end
    
    % Assegnazione dell'area a ciascun nodo
    for j = 1:num_elements
        nodes = topology(j, :);
        for k = 1:3
            f(nodes(k)) = f(nodes(k)) + element_areas(j) / 3;
        end
    end
end