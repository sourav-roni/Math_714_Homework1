% Name: Sourav Pal
% Date: 09/05/2020
% This function reflects major steps as discussed in class and present in
% the book. Uses formulae from the reference book
function x_iter = VCycle(step, all_As, all_restrictions, b, x0 , pre_steps, post_steps, n_direct_sol)
	% Coefficient matrix for current run
    A = cell2mat(all_As(step));
	
	% If the problem is small enough, solve it directly
    % This is the base case for the recursive function
	n = size(b, 1);
	if (n <= n_direct_sol)
        % Using backslash 
		x_iter = A \ b;
		return;
	end

	% Pre recursion iteration
	x_iter = gauss_sidel(A, b, pre_steps, x0);
	
	R = cell2mat(all_restrictions(step));
	
	% Compute residual
	r   = b - A * x_iter;
    
    % Coarsen
	r_H = R * r;
	
	% Initilize, here size is changing at each recursion, so 
    % can't do a simple initalization, need to pass it
	x0  = zeros(size(R, 1), 1);
    
    % Recursion to next level
	e_H = VCycle(step + 1, all_As, all_restrictions, r_H, x0, pre_steps, post_steps, n_direct_sol);
	
	% Error correction, 4 because of 2D, definition from book
	x_iter = x_iter + 4 * transpose(R) * e_H;
	
	% Post recursion iteration
    % After this step, the value of x_iter is returned to the caller of the
    % recursive function
	x_iter = gauss_sidel(A, b, post_steps, x_iter);
end