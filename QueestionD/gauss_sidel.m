% Name: Sourav Pal
% Date: 09/05/2020
%Implementation of Gauss Sidel Iterative method
function x = gauss_sidel(A, b, max_iter, x0)
    % Note A is a matrix and b and x are vectors here
	n = size(A, 1);
    % Initialization
    x = x0;
    % Used for convergence of the iterative method
    iter_tol = 1e-10;

	for i=1:max_iter
		x_prev = x;	
		for j = 1 : n
			x(j) = (b(j) - A(j, :) * x + A(j,j) * x(j))/A(j,j);
        end		
		err = max(abs(x_prev - x));
        if err < iter_tol
            break
        end
	end
end