addpath('../boundary_conditions')
function bb = fem2d_quad_bc_integral(x, i, j)
% FEM2D_QUAD_BC_INTEGRAL performs the line integral 
% \int_{\partial \Omega} g * \phi_j d s on the boundary edge x(i)-->x(j)
% of a given quadrilateral element.
%
% [IN]  x    : 2 * 4 matrix, the geometric coordinates of the element's nodes,
%              first row is x coordinates, second row is y coordinates,
%              points should be in counter clockwise order
% [IN]  i, j : The index of the nodes on the boundary edge
% [OUT] bb   : 4 * 1 matrix, contains the result of the line integral

    % Error handling: ensure i and j are valid indices
    if i < 1 || i > 4 || j < 1 || j > 4 || i == j
        error('Indices i and j must be within the range 1 to 4 and must not be equal.');
    end
    
    % Gauss-Legendre quadrature points and weights (4-point quadrature)
    qx = [-0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 0.8611363115940526];
    qw = [ 0.6521451548625461, 0.6521451548625461,  0.3478548451374538, 0.3478548451374538];
    n_quadrature = length(qx);
    
    % Compute the length of the boundary edge
    dx = x(1, j) - x(1, i);
    dy = x(2, j) - x(2, i);
    l  = sqrt(dx * dx + dy * dy);
    semi_bma = l / 2;
    semi_bpa = l / 2;
    
    % Initialize the result vector
    bb = zeros(4, 1);
    b_i = 0; b_j = 0;
    
    % Perform the quadrature integration
    for iq = 1:n_quadrature
        % Convert from [-1, 1] to [0, l] to evaluate linear basis functions
        t = semi_bma * qx(iq) + semi_bpa;
        t_over_l = t / l;
        
        % Evaluate the function g(x, y) at the quadrature point on the boundary edge
        g_x = x(1, i) + t_over_l * dx;
        g_y = x(2, i) + t_over_l * dy;
        gxy = poisson2d_robin_bc_g(g_x, g_y);
        
        % Accumulate the contributions to the integral
        b_i = b_i + ((1 - t_over_l) * gxy) * qw(iq);
        b_j = b_j + (t_over_l       * gxy) * qw(iq);
    end
    
    % Store the results in the output vector
    bb(i) = semi_bma * b_i;
    bb(j) = semi_bpa * b_j;
end
