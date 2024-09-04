addpath('../boundary_conditions')
function bb = fem2d_tri_integral_g(x, i, j)
% Perform the line integral \int_{\partial \Omega} g * \phi_j \, ds
% along the boundary edge from node i to node j of a triangular element.
%
% INPUTS:
%   x   - 2x3 matrix of geometric coordinates of the triangle vertices.
%         First row contains x-coordinates, second row contains y-coordinates.
%         Vertices should be listed in counter-clockwise order.
%   i, j - Indices of the nodes on the boundary edge.
%
% OUTPUT:
%   bb  - 3x1 matrix containing the result of the line integral.

    % Gaussian quadrature points and weights
    qx = [-0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 0.8611363115940526];
    qw = [0.6521451548625461, 0.6521451548625461, 0.3478548451374538, 0.3478548451374538];
    n_quadrature = length(qx);
    
    % Compute boundary edge length
    dx = x(1, j) - x(1, i);
    dy = x(2, j) - x(2, i);
    l = sqrt(dx^2 + dy^2);
    
    % Pre-compute semi-lengths
    semi_bma = l / 2;
    semi_bpa = l / 2;
    
    % Initialize integral results
    bi = 0;
    bj = 0;
    
    % Perform Gaussian quadrature integration
    for iq = 1:n_quadrature
        % Convert quadrature point from [-1, 1] to [0, l]
        t = semi_bma * qx(iq) + semi_bpa;
        t_over_l = t / l;
        
        % Compute g(x, y) at the quadrature point
        g_x = x(1, i) + t_over_l * dx;
        g_y = x(2, i) + t_over_l * dy;
        gxy = poisson2d_robin_bc_g(g_x, g_y);
        
        % Accumulate the integral results
        bi = bi + ((1 - t_over_l) * gxy) * qw(iq);
        bj = bj + (t_over_l * gxy) * qw(iq);
    end
    
    % Assign results to the output vector
    bb = zeros(3, 1);
    bb(i) = semi_bma * bi;
    bb(j) = semi_bpa * bj;
end
