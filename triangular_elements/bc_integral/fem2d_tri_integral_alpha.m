addpath('../boundary_conditions')
function k2 = fem2d_tri_integral_alpha(x, i, j)
% Perform the line integral $\int_{\partial \Omega} \alpha * \phi_i * \phi_j d s$ 
% on the boundary edge x(i)-->x(j) of a given triangular element.
%
% INPUTS:
%   x    - 2x3 matrix of the geometric coordinates of the triangle's vertices.
%          The first row contains x-coordinates, and the second row contains y-coordinates.
%          Vertices should be listed in counter-clockwise order.
%   i, j - Indices of the nodes on the boundary edge.
%
% OUTPUT:
%   k2   - 3x3 matrix containing the result of the line integral, with entries:
%          k2(i, i), k2(i, j), k2(j, i), k2(j, j).

    % Gauss quadrature points and weights for 1D integration
    qx = [-0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 0.8611363115940526];
    qw = [ 0.6521451548625461, 0.6521451548625461,  0.3478548451374538, 0.3478548451374538];
    n_quadrature = numel(qx);

    % Calculate the length of the boundary edge
    dx = x(1, j) - x(1, i);
    dy = x(2, j) - x(2, i);
    l = sqrt(dx^2 + dy^2);
    semi_bma = l / 2;
    semi_bpa = l / 2;
    
    % Initialize result matrix
    k2 = zeros(3, 3);
    fii = 0; fij = 0; fjj = 0;

    % Perform 1D integration using Gauss quadrature
    for iq = 1 : n_quadrature
        % Map quadrature points from [-1, 1] to [0, l]
        t = semi_bma * qx(iq) + semi_bpa;
        t_over_l = t / l;
        
        % Compute the coordinate on the line segment
        alpha_x = x(1, i) + t_over_l * dx;
        alpha_y = x(2, i) + t_over_l * dy;
        alpha = poisson2d_robin_bc_alpha(alpha_x, alpha_y);
        
        % Update integral values for basis functions
        fii = fii + ((1 - t_over_l)^2 * alpha) * qw(iq);
        fij = fij + ((1 - t_over_l) * t_over_l * alpha) * qw(iq);
        fjj = fjj + (t_over_l^2 * alpha) * qw(iq);
    end
    
    % Populate the result matrix
    k2(i, i) = semi_bma * fii;
    k2(i, j) = semi_bma * fij;
    k2(j, i) = semi_bma * fij;
    k2(j, j) = semi_bma * fjj;
end
