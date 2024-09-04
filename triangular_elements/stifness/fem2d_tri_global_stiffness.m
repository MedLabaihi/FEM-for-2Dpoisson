function K = fem2d_tri_global_stiffness(coords, ien)
% Generate the global stiffness matrix for $-\nabla^2 u$ in 2D case
% using linear basis functions and triangular elements.
% [IN]  coords : n * 2 matrix, n is the number of grid points,
%                each row contains the x and y coordinates of a grid point.
% [IN]  ien    : m * 3 matrix, m is the number of triangular elements,
%                each row contains the IDs of the three vertex points of the element,
%                listed in counter-clockwise order.
% [OUT] K      : n * n Global stiffness matrix.
    
    n = size(coords, 1);
    m = size(ien, 1);
    
    K = sparse(n, n);
    
    for i_elem = 1 : m
        elem_vertex_ids = ien(i_elem, :);
        vertex_coords = coords(elem_vertex_ids, :)';
        
        % Compute the unit stiffness matrix for the current element
        k = fem2d_tri_unit_stiffness(vertex_coords);
        
        % Assemble the local stiffness matrix into the global matrix
        K(elem_vertex_ids, elem_vertex_ids) = K(elem_vertex_ids, elem_vertex_ids) + k;
    end
end
