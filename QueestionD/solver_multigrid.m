% Name: Sourav Pal
% Date: 09/05/2020
% Implementation of Multigrid solver for A * x = b 
% Uses formulae from the reference book

nx = 9;
ny = 12;

% Form the matrix representing the minus 2 D Laplacian
[A,b] = form_matrix(nx, ny);

% Saving it for VCycle compute as updating it while creating restriction
% operators
origA = A;

% Size of the matrix
size_A = nx*(ny+2);

% inital guess of all zeros
x  = zeros(size_A, 1);

% Tolerance value
tol = 1e-10;

% Solve directly for small sizes
n_direct_sol = 49; % For smaller than this restrictions can't be formed, will results in index exceeds

% Pre and post recursion iterative steps
pre_steps = 3; post_steps = 3;	

% Rough estimate of the largest grid points to start with
n = floor(sqrt(size_A));

% Construct restriction operators below
% Store them because they will be required hierarchically during recursive
% calls in the VCycle
all_As = {};
all_restrictions = {};

% Lowest coarse level is refers to the most fine grid
coarse_level = 1;

% First entry is the original matrix
all_As(coarse_level) = {A};

% We make use of the following
% In the two-dimensional case, the stencil for the FW averaging is given by
% 1/16[ 1 2 1 , 2 4 2 , 1 2 1 ] from the reference book chapter

% While construction, we move from finer to coarser grid and stop the
% construction at a point where we reach the capacity to solve directly
% which is referred by the variable n_direct_sol

while(1)
    % One level coarser, number of points on one-axis and the entire
    next_n = floor(n/2);
    next_total_N = next_n^2;
    
    % Start with all zeros
    R = zeros(next_total_N, size_A);
    
    % grid point counter
    k = 0;
    
    for i = 2 : 2 : n
        for j = 2 : 2 : n
            k = k + 1;
            next_k_grid = j + n*(i-1);

            R(k, next_k_grid - n - 1) = 1/16;   R(k, next_k_grid - n   )  = 1/8;  R(k, next_k_grid - n + 1) = 1/16;

            R(k, next_k_grid - 1)     = 1/8;    R(k, next_k_grid    )     = 1/4;  R(k, next_k_grid + 1)     = 1/8;

            R(k, next_k_grid + n - 1) = 1/16;   R(k, next_k_grid + n   )  = 1/8;  R(k, next_k_grid + n + 1) = 1/16;
        end
    end

    all_restrictions(coarse_level) = {R};

    % Get next coefficient matrix, 4 because of 2D, definition from book
    A = R * A * 4 * transpose(R);

    %Next level is more coarser as we are reducing the number of grid
    %points
    coarse_level = coarse_level + 1;
    
    all_As(coarse_level) = {A};
    
    n = next_n;
    size_A = next_total_N;
    if size_A <= n_direct_sol
        break
    end
end

% Now run the VCycle and keep a VCycle counter
vc_count  = 0;
max_vcycles = 1000;
for i=1:max_vcycles
    vc_count = vc_count + 1;
    x  = VCycle(1, all_As, all_restrictions, b, x, pre_steps, post_steps, n_direct_sol);
    
    % Compute final residue and find error in terms of 2-norm
    % Note here we use the original A matrix, we saved it earlier!!
    residue  = b - origA * x;
    residue_norm = norm(residue, 2);
    
    %normalize by norm of b
    residue_norm = residue_norm/norm(b);
    
    if residue_norm < tol
        break
    end
end

fprintf('Total VCycle count required= %d\n', vc_count);
fprintf('2-norm error= %e\n', residue_norm);