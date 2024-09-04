function dtm = fem2d_quad_jacobian_det(qx, qy, x)
% FEM2D_QUAD_JACOBIAN_DET computes the determinant of the Jacobian
% at a given quadrature point in a quadrilateral element.
%
% [IN]  qx, qy   : The given quadrature point (to evaluate the function value)
% [IN]  x        : 2 * 4 matrix, the geometric coordinates of the element,
%                  first row is x coordinates, second row is y coordinates,
%                  points should be in counter-clockwise order
% [OUT] dtm      : The determinant of the Jacobian
	
	% Check the input sizes
	if numel(qx) ~= 1 || numel(qy) ~= 1
		error('qx and qy should be scalars.');
	end
	if size(x, 1) ~= 2 || size(x, 2) ~= 4
		error('x should be a 2x4 matrix.');
	end
	
	% Evaluate the partial derivatives of basis functions w.r.t. xi and eta at quadrature point
	d_N1_d_xi  = -0.25 * (1.0 - qy);
	d_N1_d_eta = -0.25 * (1.0 - qx);
	d_N2_d_xi  =  0.25 * (1.0 - qy);
	d_N2_d_eta = -0.25 * (1.0 + qx);
	d_N3_d_xi  =  0.25 * (1.0 + qy);
	d_N3_d_eta =  0.25 * (1.0 + qx);
	d_N4_d_xi  = -0.25 * (1.0 + qy);
	d_N4_d_eta =  0.25 * (1.0 - qx);
	
	% Compute the partial derivatives between x, y and xi, eta
	d_x_d_xi  = d_N1_d_xi  * x(1, 1) + d_N2_d_xi  * x(1, 2) + d_N3_d_xi  * x(1, 3) + d_N4_d_xi  * x(1, 4);
	d_x_d_eta = d_N1_d_eta * x(1, 1) + d_N2_d_eta * x(1, 2) + d_N3_d_eta * x(1, 3) + d_N4_d_eta * x(1, 4);
	d_y_d_xi  = d_N1_d_xi  * x(2, 1) + d_N2_d_xi  * x(2, 2) + d_N3_d_xi  * x(2, 3) + d_N4_d_xi  * x(2, 4);
	d_y_d_eta = d_N1_d_eta * x(2, 1) + d_N2_d_eta * x(2, 2) + d_N3_d_eta * x(2, 3) + d_N4_d_eta * x(2, 4);
	
	% Compute the determinant of Jacobian
	dtm = d_x_d_xi * d_y_d_eta - d_y_d_xi * d_x_d_eta;
end
