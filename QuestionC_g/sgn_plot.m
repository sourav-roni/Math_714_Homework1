% Signum function plots
% In order to compute error in this case, we will use a fine mesh 
% and compute the error for coarser grids and then plot the error 
% in a log-log plot compared with the delta difference in the grid points

% We beging with setting the number of points in the most coarse grid below
nx = 3;
ny = 3;

% We will compute the iterative solution for 10 other grids with increasing fineness
% We also take care such that the coarser grids are indeed subset of the
% most fine grid

% Variable to account for most fine grid
finest_grid = 5;

% Function to compute the number of grid points in-between the edges 
% k accounts for the coarseness of the grid, the higher the value, the
% finer the grid
% n is the the number of points int the most coarse grid
num_pts = @(n,k) (2^k*(n+1)-1);

% Compute iterative solution for most fine grid
max_nx = num_pts(nx, finest_grid);
max_ny = num_pts(ny, finest_grid);
[fin_u, fin_diff, fin_iter] = sgn_jacobi_iter(max_nx, max_ny);

% Storage variables
all_errors = zeros(1,finest_grid+1);
all_diffs = zeros(1, finest_grid+1);
all_iters = zeros(1,finest_grid+1);

for k=0:finest_grid
    [curr_u, curr_diff, curr_iter] = sgn_jacobi_iter(num_pts(nx, k), num_pts(ny, k));    % Get values for coarse grid
    to_skip_pow = finest_grid - k;
    fin_u_downsample = fin_u(1:2^to_skip_pow:max_nx+2, 1:2^to_skip_pow:max_ny+2);        % Downsample from fine grid to compare with coarser grid
    size(fin_u_downsample)
    size(curr_u)
    all_errors(1, k+1) = max(max(abs(fin_u_downsample - curr_u)));
    all_diffs(1, k+1) = curr_diff;
    all_iters(1, k+1) = curr_iter;
end

delh = zeros(1,finest_grid+1);
for i=1:finest_grid+1
    delh(1,i) = 1/(num_pts(nx, i-1)+1);
end

l_delh = log(delh);
l_err = log(all_errors);
plot(l_delh, l_err)
xlabel('log h')
ylabel('log max-norm error')

