function [x, y] = local_to_global_coords(xi, eta, geo_coord)
% LOCAL_TO_GLOBAL_COORDS transforms local coordinates (xi, eta) to global coordinates (x, y)
% in a quadrilateral element based on the element's geometric coordinates.
%
% [IN]  xi, eta       : Local coordinates in the range [-1, 1]
% [IN]  geo_coord     : 2 * 4 matrix of geometric coordinates of the element's vertices,
%                        first row contains x coordinates, second row contains y coordinates,
%                        vertices should be listed in counter-clockwise order
% [OUT] x, y          : Global coordinates corresponding to the local coordinates (xi, eta)

    % Compute the global x coordinate from the local (xi, eta) using bilinear shape functions
    x = geo_coord(1, 1) * 0.25 * (1.0 - xi) * (1.0 - eta);
    x = x + geo_coord(1, 2) * 0.25 * (1.0 + xi) * (1.0 - eta);
    x = x + geo_coord(1, 3) * 0.25 * (1.0 + xi) * (1.0 + eta);
    x = x + geo_coord(1, 4) * 0.25 * (1.0 - xi) * (1.0 + eta);
    
    % Compute the global y coordinate from the local (xi, eta) using bilinear shape functions
    y = geo_coord(2, 1) * 0.25 * (1.0 - xi) * (1.0 - eta);
    y = y + geo_coord(2, 2) * 0.25 * (1.0 + xi) * (1.0 - eta);
    y = y + geo_coord(2, 3) * 0.25 * (1.0 + xi) * (1.0 + eta);
    y = y + geo_coord(2, 4) * 0.25 * (1.0 - xi) * (1.0 + eta);
end
