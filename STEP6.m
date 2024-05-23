clc;
close all;
clear;

% Given parameters for all test cases
nx = 101; % Number of grid points in x
ny = 101; % Number of grid points in y
nt = 161; % Number of time steps
xmax = 2; % Max distance in x
ymax = 2; % Max distance in y
dx = xmax / (nx - 1); % Grid spacing in x
dy = ymax / (ny - 1); % Grid spacing in y

% Define mesh grid for entire domain
[X, Y] = meshgrid(linspace(0, xmax, nx), linspace(0, ymax, ny));

% Test cases
test_cases_CFL = [1.1, 1, 0.5]; % C' values for each test case

% Time steps at which to plot
plot_steps = [1, 81, 161];

for test_case = 1:length(test_cases_CFL)
    CFL = test_cases_CFL(test_case);
    
    % Calculate time step based on CFL condition
    dt = CFL * dx / 2; % Corrected time step calculation
    
    % Initialize u and v fields with initial conditions
    u = ones(ny, nx);
    v = ones(ny, nx);
    
    % Apply initial condition for u and v
    for i = 1:nx
        for j = 1:ny
            if X(j, i) >= 0.5 && X(j, i) <= 1 && Y(j, i) >= 0.5 && Y(j, i) <= 1
                u(j, i) = 2;
                v(j, i) = 2;
            else
                u(j, i) = 1;
                v(j, i) = 1;
            end
        end
    end
    
    % Enforce constant boundary conditions
    u(:, 1) = 1; u(:, nx) = 1; u(1, :) = 1; u(ny, :) = 1;
    v(:, 1) = 1; v(:, nx) = 1; v(1, :) = 1; v(ny, :) = 1;

    % Time stepping
    for n = 1:nt
        un = u;
        vn = v;
        for i = 2:ny-1
            for j = 2:nx-1
                u(i,j) = un(i,j) - un(i,j) * dt/dx * (un(i,j) - un(i,j-1)) - vn(i,j) * dt/dy * (un(i,j) - un(i-1,j));
                v(i,j) = vn(i,j) - un(i,j) * dt/dx * (vn(i,j) - vn(i,j-1)) - vn(i,j) * dt/dy * (vn(i,j) - vn(i-1,j));
            end
        end
        
        % Enforce the boundary conditions again
        u(:, 1) = 1; u(:, nx) = 1; u(1, :) = 1; u(ny, :) = 1;
        v(:, 1) = 1; v(:, nx) = 1; v(1, :) = 1; v(ny, :) = 1;
        
        % Visualization at specific time steps
        if ismember(n, plot_steps)
            figure;
            surf(X, Y, u, 'EdgeColor', 'none');
            shading interp;
            title(['2D Non-Linear Convection at time step ', num2str(n-1), ' with CFL = ', num2str(CFL)]);
            xlabel('X');
            ylabel('Y');
            zlabel('Velocity field u');
            drawnow;
        end
    end
    
    % Print the Δx, Δt, and CFL for each test case
    fprintf('Test case %d: Δx = %.3f, Δt = %.5f, CFL = %.1f\n', test_case, dx, dt, CFL);
end

