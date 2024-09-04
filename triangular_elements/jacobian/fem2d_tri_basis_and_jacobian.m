function [sh, dtm] = fem2d_tri_basis_and_jacobian(qx, qy, x)
% Compute the basis function information and the Jacobian determinant for a quadrature point
% in a triangular element.
%
% INPUTS:
%   qx, qy - Coordinates of the quadrature point in the local (xi, eta) space.
%   x      - 2x3 matrix of the geometric coordinates of the triangle vertices.
%            The first row contains x-coordinates, and the second row contains y-coordinates.
%            Vertices should be listed in counter-clockwise order.
%
% OUTPUTS:
%   sh     - 3x3 matrix containing the shape function information at the quadrature point:
%            First row: dN/dx, Second row: dN/dy, Third row: N (basis functions values).
%   dtm    - Determinant of the Jacobian matrix.

    % Basis function derivatives with respect to xi and eta
    d_N1_d_xi  = -1.0; d_N1_d_eta = -1.0;
    d_N2_d_xi  =  1.0; d_N2_d_eta =  0.0;
    d_N3_d_xi  =  0.0; d_N3_d_eta =  1.0;
    
    % Compute partial derivatives of x and y with respect to xi and eta
    d_x_d_xi  = d_N1_d_xi * x(1, 1) + d_N2_d_xi * x(1, 2) + d_N3_d_xi * x(1, 3);
    d_x_d_eta = d_N1_d_eta * x(1, 1) + d_N2_d_eta * x(1, 2) + d_N3_d_eta * x(1, 3);
    d_y_d_xi  = d_N1_d_xi * x(2, 1) + d_N2_d_xi * x(2, 2) + d_N3_d_xi * x(2, 3);
    d_y_d_eta = d_N1_d_eta * x(2, 1) + d_N2_d_eta * x(2, 2) + d_N3_d_eta * x(2, 3);
    
    % Compute the determinant of the Jacobian matrix
    dtm = d_x_d_xi * d_y_d_eta - d_y_d_xi * d_x_d_eta;
    
    % Evaluate the partial derivatives of basis functions with respect to x and y
    sh(1, 1) = (d_N1_d_xi *  d_y_d_eta + d_N1_d_eta * -d_y_d_xi) / dtm;
    sh(2, 1) = (d_N1_d_xi * -d_x_d_eta + d_N1_d_eta *  d_x_d_xi) / dtm;
    sh(1, 2) = (d_N2_d_xi *  d_y_d_eta + d_N2_d_eta * -d_y_d_xi) / dtm;
    sh(2, 2) = (d_N2_d_xi * -d_x_d_eta + d_N2_d_eta *  d_x_d_xi) / dtm; 
    sh(1, 3) = (d_N3_d_xi *  d_y_d_eta + d_N3_d_eta * -d_y_d_xi) / dtm;
    sh(2, 3) = (d_N3_d_xi * -d_x_d_eta + d_N3_d_eta *  d_x_d_xi) / dtm;
    
    % Evaluate the basis functions at the quadrature point
    sh(3, 1) = 1.0 - qx - qy;
    sh(3, 2) = qx;
    sh(3, 3) = qy;
end
