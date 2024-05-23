clc;           
close all;       
clear;          

%test cases a, b, c
test_cases = {'a', 'b', 'c'};
nt_values = [301, 201, 101]; %# of time steps (nt)
nx = 101;                    %# of spatial grid points (nx)
xmax = 2.0;                  %max distance (Xmax)

for k = 1:length(test_cases)
    nt = nt_values(k); 
    dt = xmax / (nt - 1);    %Δt
    dx = xmax / (nx - 1);    %Δx
    x = linspace(0, xmax, nx);  %# of spatial grid points (nx)

    %initial condition
    u = 2 * sin(pi * x / xmax);

    umax = max(u);
    CFL = umax * dt / dx; 

    fprintf('Initial Test Case (%s): Δx = %.4f, Δt = %.4f, Umax = %.4f, CFL = %.4f\n', test_cases{k}, dx, dt, umax, CFL); 
    
    figure;
    for n = 1:nt
        un = u;
        for i = 2:nx-1
            %nonlinear convection term, using u as wave speed
            u(i) = un(i) - un(i) * dt / dx * (un(i) - un(i-1));
        end
        
        u(1) = 0;   % ENFORCE BOUNDARY CONDITION
        u(end) = 0; % ENFORCE BOUNDARY CONDITION

        %results at each time step
        plot(x, u, 'b-', 'LineWidth', 2);
        ylim([min(u) - 0.5, max(u) + 0.5]);
        title(sprintf('Inviscid Burgers Equation - Test Case (%s) at Time Step %d', test_cases{k}, n));
        xlabel('x');
        ylabel('u');
        pause(0.005);
        drawnow; 
    end
end

