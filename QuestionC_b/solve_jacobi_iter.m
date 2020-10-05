% Name: Sourav Pal
% Date: 09/05/2020
function [max_error, temp_diff, iters] = solve_jacobi_iter(nx, ny)

% Define number of grid points in the two dimensions and their appropriate
% spacings
hx = 1/(nx + 1);
hy = 1/(ny + 1);

% Generate grid points in X and Y dimensions
X = 0:hx:1;
Y = 0:hy:1;
xs = X;
ys = Y;

% Meshgrid
[X, Y] = meshgrid(Y,X);

% Allocate the u matrix to all zeroes
u = zeros(nx+2, ny+2);
unew = zeros(nx+2, ny+2);
actual = zeros(nx+2, ny+2);

% Dirichlet B.C.
for j=1:ny+2
    u(1,j) = cos(2*pi*hy*(j-1));
    unew(1,j) = cos(2*pi*hy*(j-1));
end

% Dirichlet B.C.
u(nx+2,:) = 0;
unew(nx+2,:) = 0;

% Number of iterations for Jacobi
maxiter = 1000000;

factor = (hx^2 * hy^2) / (2 * (hx^2 + hy^2)) ;

% Tolerance of Convergence for Jacobi
tol = 1.0000e-10;

for iter=1:maxiter
    for i=2:(nx+1)
        for j=1:(ny+2)
            if j==1
                unew(i,j) = factor * ( u(i+1,j)/hx^2 + u(i-1,j)/hx^2 + 2*u(i,j+1)/hy^2 );
            elseif j==ny+2
                unew(i,j) = factor * ( u(i+1,j)/hx^2 + u(i-1,j)/hx^2 + 2*u(i,j-1)/hy^2);
            else
                unew(i,j) = factor * ( ((u(i+1,j) + u(i-1,j)) / hx^2 ) + ((u(i,j+1) + u(i,j-1)) / hy^2) );
            end
        end
    end
    temp_diff_mat = unew - u;
    temp_diff_vec = reshape(temp_diff_mat, 1, []);
    temp_diff = norm(temp_diff_vec);
    if temp_diff < tol
        iters = iter;
        break
    end
    u = unew;
end

exact = @(x,y) (-(cosh(2*pi)./sinh(2*pi)).*sinh(2*pi.*x) + cosh(2*pi.*x)).*cos(2*pi.*y);
for i=1:nx+2
    for j=1:ny+2
        actual(i,j) = exact(xs(i), ys(j));
    end
end

max_error = max(max(abs(actual - u)));
