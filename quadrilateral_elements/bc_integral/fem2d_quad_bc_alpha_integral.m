addpath('../boundary_conditions')
function k2 = fem2d_quad_bc_alpha_integral(x, i, j)
% FEM2D_QUAD_BC_ALPHA_INTEGRAL performs the line integral 
% \int_{\partial \Omega} \alpha * \phi_i * \phi_j d s on the boundary edge x(i)-->x(j)
% of a given quadrilateral element.
%
% [IN]  x    : 2 * 4 matrix, the geometric coordinates of the element's nodes,
%              first row is x coordinates, second row is y coordinates,
%              points should be in counter clockwise order
% [IN]  i, j : The index of the nodes on the boundary edge
% [OUT] k2   : 4 * 4 matrix, contains the result of the line integral

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
    
    % Initialize the result matrix
    k2 = zeros(4, 4);
    integral_fii = 0; integral_fij = 0; integral_fjj = 0;
    
    % Perform the quadrature integration
    for iq = 1:n_quadrature
        % Convert from [-1, 1] to [0, l] to evaluate linear basis functions
        t = semi_bma * qx(iq) + semi_bpa;
        t_over_l = t / l;
        
        % Evaluate alpha at the quadrature point on the boundary edge
        alpha_x = x(1, i) + t_over_l * dx;
        alpha_y = x(2, i) + t_over_l * dy;
        alpha = poisson2d_robin_bc_alpha(alpha_x, alpha_y);
        
        % Accumulate the contributions to the integral
        integral_fii = integral_fii + ((1 - t_over_l) * (1 - t_over_l) * alpha) * qw(iq);
        integral_fij = integral_fij + ((1 - t_over_l) * t_over_l * alpha) * qw(iq);
        integral_fjj = integral_fjj + (t_over_l * t_over_l * alpha) * qw(iq);
    end
    
    % Store the results in the output matrix
    k2(i, i) = semi_bma * integral_fii;
    k2(i, j) = semi_bma * integral_fij;
    k2(j, i) = semi_bma * integral_fij;
    k2(j, j) = semi_bma * integral_fjj;
end
