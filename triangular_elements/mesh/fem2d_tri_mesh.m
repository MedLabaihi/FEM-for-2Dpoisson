function [coords, ien, bgp] = fem2d_tri_mesh(x, y)
% Generate a triangular mesh for the 2D finite element method on the unit square [0, 1] x [0, 1].
%
% INPUTS:
%   x  - Vector of x-coordinates of the grid points.
%   y  - Vector of y-coordinates of the grid points.
%
% OUTPUTS:
%   coords - N x 2 matrix of grid point coordinates (x, y), where N is the number of grid points.
%   ien    - M x 3 matrix of element connectivity, where M is the number of triangular elements.
%            Each row contains 3 vertex indices of the element in counter-clockwise order.
%   bgp    - Vector of boundary grid point indices.
%
% NOTE:
%   The grid is assumed to be regular, with the nodes listed in counter-clockwise order.

    nx = length(x); % Number of grid points in x direction
    ny = length(y); % Number of grid points in y direction
    
    N = nx * ny; % Total number of grid points
    M = 2 * (nx - 1) * (ny - 1); % Total number of triangular elements
    
    coords = zeros(N, 2); % Initialize grid point coordinates matrix
    ien    = zeros(M, 3); % Initialize element connectivity matrix
    bgp    = zeros(2 * (nx + ny - 2), 1); % Initialize boundary grid points vector
    
    ibgp   = 0; % Boundary grid point index
    
    for iy = 1:ny
        for ix = 1:nx
            % Calculate the index of the current grid point
            icoord = (iy - 1) * nx + ix;
            coords(icoord, :) = [x(ix), y(iy)];
            
            % Process the elements in the current grid
            if (iy < ny) && (ix < nx)
                % Element index
                ielem  = (iy - 1) * (nx - 1) + ix;
                ielem1 = 2 * ielem - 1;
                ielem2 = 2 * ielem;
                
                % Indices of the vertices of the current element
                icoord1 = icoord;
                icoord2 = icoord + 1;
                icoord3 = icoord + 1 + nx;
                icoord4 = icoord + nx;
                
                % Define the connectivity of the two triangular elements
                ien(ielem1, :) = [icoord1, icoord2, icoord3];
                ien(ielem2, :) = [icoord1, icoord3, icoord4];
            end
            
            % Identify boundary grid points
            if (iy == 1) || (iy == ny) || (ix == 1) || (ix == nx)
                ibgp = ibgp + 1;
                bgp(ibgp) = icoord;
            end
        end
    end
    
    % Trim unused space in boundary grid points vector
    bgp = bgp(1:ibgp);
end
