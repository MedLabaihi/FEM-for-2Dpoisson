function b = fem2d_tri_global_load(coords, ien)
% Assemble the global right-hand side vector for the equation -?u = f
% using linear basis functions and triangular elements.
%
% INPUT:
%   coords : n x 2 matrix, where n is the number of grid points.
%            Each row contains the x and y coordinates of a grid point.
%   ien    : m x 3 matrix, where m is the number of triangular elements.
%            Each row contains the vertex point ids of the element, listed in counter-clockwise order.
%
% OUTPUT:
%   b      : n x 1 vector representing the right-hand side load vector.

    n = size(coords, 1); % Number of grid points
    m = size(ien, 1);    % Number of elements
    
    b = zeros(n, 1);     % Initialize the RHS vector
    
    % Loop over each triangular element
    for i_elem = 1 : m
        elem_vertex_ids = ien(i_elem, :);          % Get vertex ids for the current element
        vertex_coords = coords(elem_vertex_ids, :)'; % Get the coordinates of the vertices
        
        % Compute the local load vector for the current element
        ub = fem2d_tri_rhs_load(vertex_coords);
        
        % Assemble the local load vector into the global vector
        b(elem_vertex_ids) = b(elem_vertex_ids) + ub;
    end
end
