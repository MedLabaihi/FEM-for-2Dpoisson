addpath('../jacobian')
addpath('../boundary_conditions')
function ub = fem2d_tri_rhs_load(x)
% Compute the vector of integrals for the linear basis functions on a 
% triangular element with coordinates given in the matrix x.
% This function uses Gauss quadrature to integrate the term:
% \int_{\Omega^e} f * \phi_i d\Omega, and adds contributions from
% boundary integrals if the triangle is on the boundary.
%
% INPUT:
%   x  - 2x3 matrix of the geometric coordinates of the element's vertices.
%        The first row contains x-coordinates, the second row contains y-coordinates.
%        Vertices should be listed in counter-clockwise order.
%
% OUTPUT:
%   ub - 3x1 vector, integral results for each basis function.

    % Gauss quadrature points and weights for 2D integration
    qx = [1.0/3.0, 0.6, 0.2, 0.2];
    qy = [1.0/3.0, 0.2, 0.2, 0.6];
    w  = [-27.0/96.0, 25.0/96.0, 25.0/96.0, 25.0/96.0];
    n_quadrature = numel(w);
    
    % Initialize the result vector
    ub = zeros(3, 1);
    
    % Compute the integral over the element's area
    for iq = 1 : n_quadrature
        [f_x, f_y] = fem2d_tri_local_to_global(qx(iq), qy(iq), x);
        f_xy = poisson2d_rhs_f(f_x, f_y);
        
        dtm = fem2d_tri_jacobian_det(qx(iq), qy(iq), x);
        
        % Evaluate the basis functions at quadrature point
        N = [1.0 - qx(iq) - qy(iq), qx(iq), qy(iq)];
        
        % Update the integral results
        ub = ub + dtm * f_xy * N' * w(iq);
    end
    
    % Add boundary contributions if applicable
    if is_boundary_edge(x(:, [1, 2]))
        ub = ub + fem2d_tri_integral_g(x, 1, 2);
    end
    
    if is_boundary_edge(x(:, [2, 3]))
        ub = ub + fem2d_tri_integral_g(x, 2, 3);
    end
    
    if is_boundary_edge(x(:, [3, 1]))
        ub = ub + fem2d_tri_integral_g(x, 3, 1);
    end
end