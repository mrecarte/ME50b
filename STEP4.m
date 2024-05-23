clc;
close all;
clear;

%test cases and parameters
test_cases = {'a', 'b', 'c'};
viscosities = [0.1, 0.01, 0.07];
nt_values = [151, 151, 115];
nx_values = [151, 151, 101];
tmax = 0.5;
xmax = 2 * pi;

%loop over test cases
for k = 1:length(test_cases)

    nu = viscosities(k);
    nt = nt_values(k);
    nx = nx_values(k);
    dt = tmax / (nt - 1);
    dx = xmax / (nx - 1);
    x = linspace(0, xmax, nx);
    
    % Use symbolic math to define the analytical solution and its derivative
    syms x_sym nu_sym t_sym
    phi_sym = exp(-(x_sym-4*t_sym)^2/(4*nu_sym*(t_sym+1))) + exp(-(x_sym-4*t_sym-2*pi)^2/(4*nu_sym*(t_sym+1)));
    dphi_dx_sym = diff(phi_sym, x_sym);
    u_analytical_sym = -2*nu_sym*(dphi_dx_sym/phi_sym) + 4;
    u_analytical = matlabFunction(u_analytical_sym, 'Vars', [nu_sym, t_sym, x_sym]);
    
    %initial condiion
    u = u_analytical(nu, 0, x);
    
    %plot initial condition
    figure;
    plot(x, u, 'o');
    hold on;
    
    for n = 1:nt-1
        un = u;
        
        for i = 2:nx-1
            u(i) = un(i) - un(i) * dt / dx * (un(i) - un(i-1)) + nu * dt / dx^2 * (un(i+1) - 2*un(i) + un(i-1));
        end
        
        %enforce periodic boundary conditions
        u(1) = u(nx); %u(0) = u(2*pi)
        u(nx) = u(1); %u(2*pi) = u(0)
       
        plot(x, u, 'b-', 'LineWidth', 0.1);
        ylim([min(u)-1.4, max(u)+1.8]);
        title(sprintf('1D Viscous Burgers'' Equation - Test Case (%s) at Time Step %d', test_cases{k}, n));
        xlabel('x');
        ylabel('u');
        pause(0.005);
    end
    
    %analytical solution to plot
    u_analytical_final = u_analytical(nu, (nt-1)*dt, x);
    plot(x, u_analytical_final, 'r-', 'LineWidth', 1.5);
    legend('Numerical', 'Analytical');
    
    fprintf('Test Case (%s): Δx = %.4f, Δt = %.4f, ν = %.4f\n', test_cases{k}, dx, dt, nu);
end




