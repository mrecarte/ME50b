clc;
close all;
clear;

vis = 0.1;
xmax = 2;
ymax = 2;
nx = 51;
ny = 51;

% Test cases to loop through
test_cases = {'a', 'b', 'c'};
C_values = [0.6, 0.5, 0.25];
time_steps_to_plot = [1, 50, 100]; % Time steps at which to plot

for k = 1:length(test_cases)
    C = C_values(k);
    nt = 101;
    dx = xmax / (nx - 1);
    dy = ymax / (ny - 1);
    dt = C * (dx^2) / (2 * vis);

    x = linspace(0, xmax, nx);
    y = linspace(0, ymax, ny);
    u = ones(nx, ny);
    v = ones(nx, ny); % Assuming similar dynamics for v

    % Apply initial conditions
    for i = 1:nx
        for j = 1:ny
            if x(i) >= 0.5 && x(i) <= 1 && y(j) >= 0.5 && y(j) <= 1
                u(i, j) = 2; % u-field
                v(i, j) = 2; % v-field
            end
        end
    end

    [X, Y] = meshgrid(x, y);

    for n = 1:nt
        un = u;
        vn = v; % Assuming similar dynamics for v

        for i = 2:nx-1
            for j = 2:ny-1
                u(i, j) = un(i, j) + vis * dt / dx^2 * (un(i+1, j) - 2*un(i, j) + un(i-1, j)) ...
                          + vis * dt / dy^2 * (un(i, j+1) - 2*un(i, j) + un(i, j-1));
            end
        end

        % Enforce boundary conditions
        u(1,:) = 1; u(nx,:) = 1; u(:,1) = 1; u(:,ny) = 1;

        % Plot at specific time steps
        if ismember(n, time_steps_to_plot)
            figure; % Open a new figure window for each plot
            surf(X, Y, u', 'EdgeColor', 'none');
            title(sprintf('2D Diffusion - Test Case (%s) at Time Step %d', test_cases{k}, n));
            xlabel('X');
            ylabel('Y');
            zlim([1 2.5]);
            ylim('auto'); % This automatically scales the y-axis
            drawnow;
            pause(0.8); % This will keep the plot displayed for 0.8 seconds
        end
    end

    fprintf('Test Case (%s): Δx = %.4f, Δy = %.4f, Δt = %.4f', test_cases{k}, dx, dy, dt);
end

