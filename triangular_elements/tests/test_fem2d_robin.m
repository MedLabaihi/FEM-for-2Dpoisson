addpath('../jacobian')
addpath('../boundary_conditions')
addpath('../bc_integral')
addpath('../stiffness')
addpath('../mesh')
function u = test_fem2d_robin(h)
% Test 2D FEM solver with Robin boundary conditions, h is the grid point interval

    if (h > 0)
        x = 0 : h : 1; x(end) = 1;
        y = x;
    else
        x = [0 : 0.05 : 0.09, 0.1 : 0.1 : 0.8, 0.9 : 0.05 : 1.0];
        y = [0 : 0.05 : 0.19, 0.2 : 0.1 : 0.7, 0.8 : 0.05 : 1.0];
    end
    
    nx = size(x, 2);
    ny = size(y, 2);
    
    % Generate the square mesh
    [coords, ien, bgp] = fem2d_tri_mesh(x, y);
    
    % Generate the global stiffness matrix
    K = fem2d_tri_global_stiffness(coords, ien);
    
    % Generate the right hand side
    b = fem2d_tri_rhs_load(coords, ien);
    
    % Apply the Robin boundary condition
    n_bgp = max(size(bgp));
    for i = 1 : n_bgp
        bgp_id = bgp(i);
        % Calculate alpha and g at the boundary node
        bv_x = coords(bgp_id, 1);
        bv_y = coords(bgp_id, 2);
        alpha = poisson2d_robin_bc_alpha(bv_x, bv_y);
        g = poisson2d_robin_bc_g(bv_x, bv_y);
        
        % Modify global stiffness matrix and right-hand side
        % Assuming boundary node's contribution is proportional to alpha
        K(bgp_id, :) = 0;
        K(bgp_id, bgp_id) = alpha;
        b(bgp_id) = g;
    end
    
    % Solve the linear system
    u = K \ b;
    
    % Plot the result
    xm = reshape(coords(:, 1), nx, ny);
    ym = reshape(coords(:, 2), nx, ny);
    zm = reshape(u, nx, ny);
    
    subplot(1, 2, 1);
    contourf(xm, ym, zm, 'ShowText','on');
    title('Contour Plot');
    
    subplot(1, 2, 2);
    trisurf(ien, xm, ym, zm);
    title('Surface Plot');
end
