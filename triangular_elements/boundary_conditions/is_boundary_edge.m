function is_boundary = is_boundary_edge(coords)
% Check if the edge defined by the coordinates is on the boundary of the domain.
    x_coords = coords(1, :);
    y_coords = coords(2, :);
    is_boundary = any(x_coords == 0 | x_coords == 1 | y_coords == 0 | y_coords == 1);
end