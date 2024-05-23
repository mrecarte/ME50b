function [u, v, p, X, Y] = lid_driven_cavity(nu, rho, nx, ny, nt, dt, nit)
    dx = 1 / (nx - 1);
    dy = 1 / (ny - 1);
    x = linspace(0, 1, nx);
    y = linspace(0, 1, ny);
    u = zeros(nx, ny);
    v = zeros(nx, ny);
    p = zeros(nx, ny);
    b = zeros(nx, ny);

    % Boundary conditions
    u(:, end) = 1; % Lid velocity

    % Helper function for updating pressure field
    function b = build_up_b(rho, dt, u, v, dx, dy)
        [nx, ny] = size(u);
        b = zeros(nx, ny);
        for i = 2:nx-1
            for j = 2:ny-1
                b(i, j) = rho * (1 / dt * ((u(i+1, j) - u(i-1, j)) / (2*dx) + (v(i, j+1) - v(i, j-1)) / (2*dy)) - ...
                          ((u(i+1, j) - u(i-1, j)) / (2*dx))^2 - 2 * ((u(i, j+1) - u(i, j-1)) / (2*dy) * (v(i+1, j) - v(i-1, j)) / (2*dx)) - ...
                          ((v(i, j+1) - v(i, j-1)) / (2*dy))^2);
            end
        end
    end

    % Helper function for solving pressure Poisson equation
    function p = pressure_poisson(p, dx, dy, b, nit)
        [nx, ny] = size(p);
        pn = zeros(nx, ny);
        for q = 1:nit
            pn = p;
            for i = 2:nx-1
                for j = 2:ny-1
                    p(i, j) = (((pn(i+1, j) + pn(i-1, j)) * dy^2 + (pn(i, j+1) + pn(i, j-1)) * dx^2) / ...
                              (2 * (dx^2 + dy^2)) - dx^2 * dy^2 / (2 * (dx^2 + dy^2)) * b(i, j));
                end
            end
            % Boundary conditions for pressure
            p(:, end) = p(:, end-1); % dp/dy = 0 at y = 1
            p(:, 1) = p(:, 2); % dp/dy = 0 at y = 0
            p(1, :) = p(2, :); % dp/dx = 0 at x = 0
            p(end, :) = p(end-1, :); % dp/dx = 0 at x = 1
        end
    end

    % Time-stepping loop
    for n = 1:nt
        un = u;
        vn = v;
        b = build_up_b(rho, dt, un, vn, dx, dy);
        p = pressure_poisson(p, dx, dy, b, nit);
        
        % Update velocity fields using discretized Navier-Stokes equations
        for i = 2:nx-1
            for j = 2:ny-1
                u(i, j) = un(i, j) ...
                    - un(i, j) * dt / dx * (un(i, j) - un(i-1, j)) ...
                    - vn(i, j) * dt / dy * (un(i, j) - un(i, j-1)) ...
                    - dt / (2 * rho * dx) * (p(i+1, j) - p(i-1, j)) ...
                    + nu * (dt / dx^2 * (un(i+1, j) - 2 * un(i, j) + un(i-1, j)) ...
                    + dt / dy^2 * (un(i, j+1) - 2 * un(i, j) + un(i, j-1)));
                
                v(i, j) = vn(i, j) ...
                    - un(i, j) * dt / dx * (vn(i, j) - vn(i-1, j)) ...
                    - vn(i, j) * dt / dy * (vn(i, j) - vn(i, j-1)) ...
                    - dt / (2 * rho * dy) * (p(i, j+1) - p(i, j-1)) ...
                    + nu * (dt / dx^2 * (vn(i+1, j) - 2 * vn(i, j) + vn(i-1, j)) ...
                    + dt / dy^2 * (vn(i, j+1) - 2 * vn(i, j) + vn(i, j-1)));
            end
        end
        
        % Boundary conditions for velocity
        u(:, 1) = 0;
        u(:, end) = 1;
        u(1, :) = 0;
        u(end, :) = 0;
        
        v(:, 1) = 0;
        v(:, end) = 0;
        v(1, :) = 0;
        v(end, :) = 0;
    end
    
    % Generate meshgrid for plotting
    [X, Y] = meshgrid(0:dx:1, 0:dy:1);
end
