% This code is inspired from convergence benchmark code provided in class

function [Delta, f] = form_matrix(nx, ny)
dx = 1/(nx+1);
dy = 1/(ny+1);

x = 0:dx:1;
y = 0:dy:1;

% Building the 1D Dirichelet minus Laplacian matrix
e = ones(nx,1);
Asp = spdiags([-e 2*e -e], -1:1, nx, nx);
  
% Building the 1D Neumann minus Laplacian matrix
e = ones(ny+2,1);
Bsp = spdiags([-e 2*e -e], -1:1, ny+2, ny+2);
Bsp(1,1:3) = [1.5 -2 0.5 ]*dy;
Bsp(end,end-2:end) = [ 0.5 -2 1.5 ]*dy;
    
% creating the identities (here carefull with the boundaries)
I_A = speye(ny+2,ny+2);   I_A(1,1) = 0;   I_A(end,end) = 0;
I_B = speye(nx,nx);
    
% assembling the 2D minus Laplacian
Delta = kron(Asp/dx^2,I_A) + kron(I_B,Bsp/dy^2);
    
% writing the source term
f = cos(2*pi*y);    f(1) = 0;   f(end) = 0;
e_1 = zeros(nx,1);   e_1(1) = 1;
F = kron(e_1, f).';
f = F(:)/dx^2;

end