function [H, M] = compute_local_matrices(topology, coordinates)
    num_nodes = size(coordinates, 1);
    num_elements = size(topology, 1);
    H = zeros(num_nodes); % Global stiffness matrix
    M = zeros(num_nodes); % Global mass matrix
    
    for e = 1:num_elements
        element_nodes = topology(e, :);
     % Node coordinates for the current element
        xi = coordinates(element_nodes(1), 1);
        yi = coordinates(element_nodes(1), 2);
        xj = coordinates(element_nodes(2), 1);
        yj = coordinates(element_nodes(2), 2);
        xm = coordinates(element_nodes(3), 1);
        ym = coordinates(element_nodes(3), 2);
       % Calculate the area of the triangle (Delta)
        delta = 0.5 * det([1, xi, yi; 1, xj, yj; 1, xm, ym]);
       % Calculate coefficients for the stiffness matrix
        ai = xj*ym - xm*yj;
        bi = yj - ym;
        ci = xm - xj;
        aj = xm*yi - xi*ym;
        bj = ym - yi;
        cj = xi - xm;
        am = xi*yj - xj*yi;
        bm = yi - yj;
        cm = xj - xi;
        % Local stiffness matrix Hloc
        Hloc = (1/(4 *delta)) * [bi*bi+ci*ci, bi*bj+ci*cj, bi*bm+ci*cm; 
                                 bj*bi+cj*ci, bj*bj+cj*cj, bj*bm+cj*cm; 
                                 bm*bi+cm*ci, bm*bj+cm*cj, bm*bm+cm*cm];
        
        % Local mass matrix Mloc
        Mloc = (delta / 12) * [2, 1, 1; 1, 2, 1;  1, 1, 2];
        
        % Assemble local matrices into global matrices
        for i = 1:3
            for j = 1:3
                H(element_nodes(i), element_nodes(j)) = H(element_nodes(i), element_nodes(j)) + Hloc(i, j);
                M(element_nodes(i), element_nodes(j)) = M(element_nodes(i), element_nodes(j)) + Mloc(i, j);
            end
        end
    end
end
        
