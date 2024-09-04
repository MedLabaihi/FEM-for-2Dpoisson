function flag = is_boundary_condition_met(x, i, j)
    % Helper function to check if the edge (i, j) meets the boundary condition
    flag = (x(1, i) == 0 && x(1, j) == 0) || ... % x == 0
           (x(1, i) == 1 && x(1, j) == 1) || ... % x == 1
           (x(2, i) == 0 && x(2, j) == 0) || ... % y == 0
           (x(2, i) == 1 && x(2, j) == 1);    % y == 1
end