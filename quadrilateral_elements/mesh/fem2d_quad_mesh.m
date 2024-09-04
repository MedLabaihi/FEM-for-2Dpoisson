function [coords, ien, bgp] = fem2d_quad_mesh(x, y)
% Generate a rectangular mesh for 2D finite element analysis.
% [IN]  x, y   : Vectors defining grid points in x and y directions
% [OUT] coords : N x 2 matrix, where N is the number of grid points.
%                Each row contains the x and y coordinates of a grid point.
% [OUT] ien    : M x 4 matrix, where M is the number of quadrilateral elements.
%                Each row contains 4 vertex point IDs of the element, in 
%                counterclockwise order.
% [OUT] bgp    : Vector containing the IDs of the grid points on the boundary.
	
	nx = length(x);
	ny = length(y);
	
	N = nx * ny;               % Total number of grid points
	M = (nx - 1) * (ny - 1);   % Total number of quadrilateral elements
	
	coords = zeros(N, 2);      % Initialize coordinates matrix
	ien    = zeros(M, 4);      % Initialize element connectivity matrix
	bgp    = zeros(2 * (nx + ny - 2), 1); % Initialize boundary grid points vector
	
	ibgp   = 0;                % Index for boundary grid points
	for iy = 1 : ny
		for ix = 1 : nx
			% Compute the index for the current grid point
			icoord = (iy - 1) * nx + ix;
			coords(icoord, :) = [x(ix), y(iy)];
			
			% Handle the element that is to the top-right of the current grid point
			if (iy < ny) && (ix < nx)
				% Compute the index for the current element
				ielem = (iy - 1) * (nx - 1) + ix;
				
				% Compute the IDs of the 4 vertices in the element, counterclockwise order
				icoord1 = icoord;
				icoord2 = icoord + 1;
				icoord3 = icoord + 1 + nx;
				icoord4 = icoord + nx;
				
				% Assign the vertex IDs to the element connectivity matrix
				ien(ielem, :) = [icoord1, icoord2, icoord3, icoord4];
			end
			
			% Check if the current grid point is on the boundary
			if (iy == 1) || (iy == ny) || (ix == 1) || (ix == nx)
				ibgp = ibgp + 1;
				bgp(ibgp) = icoord;
			end
		end
	end
end
