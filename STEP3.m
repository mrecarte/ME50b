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
dy = ymax / (ny - 1); % Grid spacing in y, should be the same as dx

% Test cases
test_cases_CFL = [1.1, 1, 0.5]; % C' values for each test case

for test_case = 1:length(test_cases_CFL)
    CFL = test_cases_CFL(test_case);
    
    % Calculate time step based on CFL condition
    dt = CFL * dx / 2; % Corrected time step calculation
    
    % Initialize u and v fields with initial conditions
    u = ones(ny, nx); % Initialize u to 1 everywhere
    v = ones(ny, nx); % Initialize v to 1 everywhere
    
    % Apply initial condition for u and v
    [X, Y] = meshgrid(linspace(0, xmax, nx), linspace(0, ymax, ny));
    u(X >= 0.5 & X <= 1 & Y >= 0.5 & Y <= 1) = 2;
    v(X >= 0.5 & X <= 1 & Y >= 0.5 & Y <= 1) = 2;
    
    % Enforce constant boundary conditions
    u(:, [1, nx]) = 1;
    u([1, ny], :) = 1;
    v(:, [1, nx]) = 1;
    v([1, ny], :) = 1;

    % Time stepping
    for n = 1:nt
        un = u;
        vn = v;
        for i = 2:(ny-1)
            for j = 2:(nx-1)
                % Apply update formulas for u and v here
                % Note: These need to be updated to include your model's specific update logic.
            end
        end
        
        % Enforce the boundary conditions
        u(:, [1, nx]) = 1;
        u([1, ny], :) = 1;
        v(:, [1, nx]) = 1;
        v([1, ny], :) = 1;
        
        % Visualization
        if any(isnan(u), 'all') || any(isinf(u), 'all')
            warning('Solution array contains NaN or Inf, skipping plot for time step %d', n);
        else
            surf(X, Y, u, 'EdgeColor', 'none');
            shading interp;
            title(['2D Non-Linear Convection at time step ', num2str(n), ' with CFL = ', num2str(CFL)]);
            xlabel('X');
            ylabel('Y');
            zlabel('Velocity field u');
            xlim([0 xmax]);
            ylim([0 ymax]);
            zlim([1 2.5]); % Adjust the limits as necessary based on your data range
            drawnow;
        end
    end
    
    % Print the dx, dt, and CFL for each test case
    fprintf('Test case %d: Δx = %.3f, Δt = %.5f, CFL = %.1f\n', test_case, dx, dt, CFL);
end


