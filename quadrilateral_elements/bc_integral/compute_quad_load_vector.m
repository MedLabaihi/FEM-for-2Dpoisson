addpath('../jacobian')
addpath('../boundary_conditions')
addpath('../bc_integral')
function ub = compute_quad_load_vector(x)
% COMPUTE_QUADRILATERAL_LOAD_VECTOR computes the integral of the load function 
% f(x, y) over a quadrilateral element and incorporates boundary contributions.
%
% [IN]  x  : 2 * 4 matrix of the geometric coordinates of the element's vertices,
%            1st row contains x coordinates, 2nd row contains y coordinates,
%            vertices should be listed in counter-clockwise order.
% [OUT] ub : 4 * 1 vector containing the integral result

    % Gauss-Legendre quadrature weights and points
    w = [0.6521451548625461, 0.6521451548625461, 0.3478548451374538, 0.3478548451374538];
    q = [-0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 0.8611363115940526];
    
    % Initialize the load vector
    ub = zeros(4, 1);
    
    % Integrate the load function f over the element using Gauss quadrature
    for ix = 1:4
        for iy = 1:4
            [f_x, f_y] = local_to_global_coords(q(ix), q(iy), x); % Updated function name
            f_xy = poisson2d_rhs_f(f_x, f_y); % Evaluate f(x, y)
            
            dtm = fem2d_quad_jacobian_det(q(ix), q(iy), x); % Compute determinant of the Jacobian
            
            % Evaluate shape functions
            N = [0.25 * (1.0 - q(ix)) * (1.0 - q(iy)), ...
                 0.25 * (1.0 + q(ix)) * (1.0 - q(iy)), ...
                 0.25 * (1.0 + q(ix)) * (1.0 + q(iy)), ...
                 0.25 * (1.0 - q(ix)) * (1.0 + q(iy))];
            
            % Accumulate the integral results
            for j = 1:4
                ub(j) = ub(j) + dtm * f_xy * N(j) * w(ix) * w(iy);
            end
        end
    end
    
    % Check boundary conditions and add boundary contributions
    edges = [1 2; 2 3; 3 4; 4 1]; % List of edges
    
    for e = 1:4
        i = edges(e, 1);
        j = edges(e, 2);
        
        if is_boundary_condition_met(x, i, j) % Updated function name
            ub = ub + fem2d_quad_bc_integral(x, i, j); % Updated function name
        end
    end
end
