function [x, y] = fem2d_tri_local_to_global(xi, eta, geo_coord)
% Transform local coordinates (xi, eta) to global coordinates (x, y) in a triangular element.
%
% The transformation is based on the shape functions of the triangular element:
% x(xi, eta) = x1 * (1 - xi - eta) + x2 * xi + x3 * eta
% y(xi, eta) = y1 * (1 - xi - eta) + y2 * xi + y3 * eta
%
% INPUTS:
%   xi         - Local coordinate xi
%   eta        - Local coordinate eta
%   geo_coord  - 2x3 matrix of the geometric coordinates of the triangle vertices
%                The first row contains x-coordinates, and the second row contains y-coordinates.
%                Vertices should be listed in counter-clockwise order.
%
% OUTPUTS:
%   x          - Global x-coordinate at (xi, eta)
%   y          - Global y-coordinate at (xi, eta)

    % Extract geometric coordinates
    x1 = geo_coord(1, 1); x2 = geo_coord(1, 2); x3 = geo_coord(1, 3);
    y1 = geo_coord(2, 1); y2 = geo_coord(2, 2); y3 = geo_coord(2, 3);

    % Compute global coordinates from local coordinates
    x = x1 * (1 - xi - eta) + x2 * xi + x3 * eta;
    y = y1 * (1 - xi - eta) + y2 * xi + y3 * eta;
end
