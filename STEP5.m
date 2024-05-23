clc
close all;
clear;

test_cases = {'a', 'b', 'c'};
C_values = [1.1, 1.0, 0.5];

nx = 81;          %spatial grid points (nx)
ny = 81;          %spatial grid points (ny)
xmax = 2;         %max distance in x
ymax = 2;         %max distance in y
dx = xmax / (nx - 1); %gridin x
dy = ymax / (ny - 1); %grid in y

%for loop for looking through each test case
for k = 1:length(test_cases)
    nt = 100;
    c = 1;
    sigma = C_values(k); 
    dt = sigma * dx / (2*c);
    CFL = 2* c * dt / dx;
    
    u = ones(ny, nx); 
    v = ones(ny, nx); 
    
    %init conditions given for step 5
    u(0.5/dy : 1/dy+1, 0.5/dx : 1/dx+1) = 2;
    v(0.5/dy : 1/dy+1, 0.5/dx : 1/dx+1) = 2;
    
    %ENFORCE BOUNDARY CONDITIONS
    u(:, [1, nx]) = 1;
    u([1, ny], :) = 1;
    v(:, [1, nx]) = 1;
    v([1, ny], :) = 1;
    
    un = u;
    vn = v;
    
    for n = 1:nt
        un = u;
        vn = v;
        for i = 2:(ny-1)
            for j = 2:(nx-1)
                u(i, j) = un(i, j) - c * dt / dx * (un(i, j) - un(i, j-1)) - c * dt / dy * (vn(i, j) - vn(i-1, j));
                v(i, j) = vn(i, j) - c * dt / dx * (un(i, j) - un(i, j-1)) - c * dt / dy * (vn(i, j) - vn(i-1, j));
            end
        end
        
        %ENFORCE BOUNDRY CONDITION
        u(:, [1, nx]) = 1;
        u([1, ny], :) = 1;
        v(:, [1, nx]) = 1;
        v([1, ny], :) = 1;
        
        figure(k);
        surf(linspace(0, xmax, nx), linspace(0, ymax, ny), u);
        shading interp;
        title(['2D Linear Convection - Test Case ', test_cases{k}, ' at Time Step ', num2str(n)]);
        xlabel('X axis');
        ylabel('Y axis');
        zlabel('Wave amplitude');
        drawnow;
        pause(1.50);
    end
    
    fprintf('Test Case (%s): Δx = %.4f, Δt = %.4f, CFL = %.4f\n', test_cases{k}, dx, dt, CFL);
end