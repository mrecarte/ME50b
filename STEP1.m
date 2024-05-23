clc
close all;
clear;

%test cases a, b, c
test_cases = {'a', 'b', 'c'};  
nt_values = [201, 101, 100];   %# of time steps (nt)
nx = 201;                      %# of spatial grid points (nx)
xmax = 2.0;                    %max distance (Xmax)
c = 0.5;                       %wave speed (c)
x = linspace(0, xmax, nx);     %# of spatial grid points (nx)

for k = 1:length(nt_values)
    nt = nt_values(k);
    dt = 2.0 / (nt - 1);       %Δt
    dx = xmax / (nx - 1);      %Δx
    CFL = c * dt / dx;         %courant number (CFL)

    %initial conditions
    u = ones(1, nx);
    u(0.5 <= x & x < 1) = 2;   % Square wave initial condition
    u(1) = 1;                  % Boundary condition at x = 0
    u(end) = 1;                % Boundary condition at x = xmax

    figure;

    for n = 1:nt
        un = u;
        for i = 2:nx-1
            %finite difference scheme for convection
            u(i) = un(i) - c * dt / dx * (un(i) - un(i-1));
        end
        u(1) = 1;   %ENFORCE PERIODIC BOUNDARY CONDITION at x = 0
        u(end) = 1; %ENFORCE PERIODIC BOUNDARY CONDITION at x = xmax

        %plot at each time step
        plot(x, u, 'b-', 'LineWidth', 2);
        ylim([-3, 5]);
        title(sprintf('1D Linear Convection - Test Case (%s) at Time Step %d', test_cases{k}, n));
        xlabel('x');
        ylabel('u');
        pause(3.0);
    end

    fprintf('Test Case (%s): Δx = %.4f, Δt = %.4f, CFL = %.4f\n', test_cases{k}, dx, dt, CFL);
end
