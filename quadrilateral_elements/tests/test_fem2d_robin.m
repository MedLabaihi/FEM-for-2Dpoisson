addpath('../jacobian')
addpath('../boundary_conditions')
addpath('../bc_integral')
addpath('../stiffness')
addpath('../mesh')
function u = test_fem2d_robin(h)
% Test 2D FEM solver with Robin boundary conditions
% h is the grid point interval

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
    [coords, ien, bgp] = fem2d_quad_mesh(x, y);
    
    % Generate the global stiffness matrix
    K = fem2d_quad_global_stiff(coords, ien);
    
    % Generate the right-hand side
    b = fem2d_quad_global_rhs_load(coords, ien);
    
    % Apply the Robin boundary condition
    % Assume the Robin boundary condition uses alpha(x, y) and g(x, y)
    n_bgp = max(size(bgp));
    for i = 1 : n_bgp
        bgp_id = bgp(i);
        bv_x = coords(bgp_id, 1);
        bv_y = coords(bgp_id, 2);
        
        % Modify K and b for Robin BC
        alpha_val = poisson2d_robin_bc_alpha(bv_x, bv_y);
        g_val = poisson2d_robin_bc_g(bv_x, bv_y);
        
        if alpha_val > 0
            K(bgp_id, :) = 0;
            K(bgp_id, bgp_id) = alpha_val;
            b(bgp_id) = g_val;
        end
    end
    
    % Solve the linear system
    u = K \ b;
    
    % Plot the result
    xm = reshape(coords(:, 1), nx, ny);
    ym = reshape(coords(:, 2), nx, ny);
    zm = reshape(u, nx, ny);
    
    subplot(1, 2, 1);
    contourf(xm, ym, zm, 'ShowText', 'on');
    title('Robin BC - Contour Plot');
    
    subplot(1, 2, 2);
    surf(xm, ym, zm);
    title('Robin BC - Surface Plot');
end
