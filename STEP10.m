clc;
close all;
clear;

nx = 51; 
ny = 26; 
nit_a = 100; 
nit_b = 1000; 
nit_c = 10000; 
xmax = 2; 
ymax = 1;
dx = xmax/(nx-1); 
dy = ymax/(ny-1);
x = 0:dx:xmax;
y = 0:dy:ymax;

%init pressure field P
P = zeros(nx, ny);

%init source term b
b = zeros(nx, ny);
b(floor(nx/4), floor(ny/4)) = 100;
b(floor(nx*3/4), floor(ny*3/4)) = -100;

%solvr each test case
P_a = solve_poisson(nx, ny, nit_a, dx, dy, P, b);
P_b = solve_poisson(nx, ny, nit_b, dx, dy, P, b);
P_c = solve_poisson(nx, ny, nit_c, dx, dy, P, b);

[X, Y] = meshgrid(x, y);

figure;
subplot(1, 3, 1);
surf(X, Y, P_a');
title('Pressure field after 100 iterations');
xlabel('x');
ylabel('y');
zlabel('Pressure');

subplot(1, 3, 2);
surf(X, Y, P_b');
title('Pressure field after 1000 iterations');
xlabel('x');
ylabel('y');
zlabel('Pressure');

subplot(1, 3, 3);
surf(X, Y, P_c');
title('Pressure field after 10000 iterations');
xlabel('x');
ylabel('y');
zlabel('Pressure');

%function to solve the Poisson equation using the provided pseudocode logic
%from thr class notes
function P = solve_poisson(nx, ny, nit, dx, dy, P, b)
    for nt = 1:nit
        Pd = P;
        for i = 2:nx-1
            for j = 2:ny-1
                P(i, j) = 0.25 * (Pd(i+1, j) + Pd(i-1, j) + Pd(i, j+1) + Pd(i, j-1)) - b(i, j) * dx^2;  % Δx = Δy
            end
        end
        P(2:nx-1, 1) = P(2:nx-1, 2);  % ∂P/∂y @ y = 0
        P(2:nx-1, ny) = P(2:nx-1, ny-1);  % ∂P/∂y @ y = ymax
    end
end

