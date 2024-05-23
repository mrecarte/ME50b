clc;
close all;
clear;

nx = 51;
ny = 26;
xmax = 2;
ymax = 1;
dx = xmax / (nx - 1);
dy = ymax / (ny - 1);
x = 0:dx:xmax; %x grid
y = 0:dy:ymax; %y grid
[X, Y] = meshgrid(x, y);

test_cases = [100, 1000, 10000];
for t = 1:length(test_cases)
    nit = test_cases(t);
    
    P = zeros(nx, ny);
    P(end, :) = y; %Dirichlet BC: P = y at x = 2
    P(1, :) = 0; %Dirichlet BC: P = 0 at x = 0

    for itt = 1:nit
        pd = P;
        for i = 2:nx-1
            for j = 2:ny-1
                P(i, j) = 0.25 * (pd(i+1, j) + pd(i-1, j) + pd(i, j+1) + pd(i, j-1));
                %applies when Δx = Δy
            end
        end

        %Neumann BC
        P(2:nx-1, 1) = P(2:nx-1, 2); % ∂P/∂y @ y = 0
        P(2:nx-1, ny) = P(2:nx-1, ny-1); % ∂P/∂y @ y = ymax
    end

    figure
    surf(X, Y, P')
    title(['Numerical Solution (nit = ', num2str(nit), ')'])
    xlabel('x')
    ylabel('y')
    zlabel('P')
end

%analytical Solution
P_coeffs = zeros(nx, ny);
series = zeros(nx, ny);
P_analytical = zeros(nx, ny);

for m = 1:2:99 %m odd from 1-99
    for i = 1:nx
        for j = 1:ny
            P_coeffs(i, j) = (sinh(m * pi * x(i)) * cos(m * pi * y(j))) / ...
                             (((m * pi)^2) * sinh(2 * pi * m));
            series(i, j) = series(i, j) + P_coeffs(i, j);
        end
    end
end

for i = 1:nx
    for j = 1:ny
        P_analytical(i, j) = (x(i) / 4) - 4 * series(i, j);
    end
end

%Extra credit (1): Generate a surface plot of the analytical solution (use m= 100 to approximate)
figure
surf(X, Y, P_analytical')
title('Analytical Solution (m=100)')
xlabel('x')
ylabel('y')
zlabel('P')

%Extra credit (2): Generate curves of p vs. x at fixed values of y for test cases a and b

%(a): nit = 100
nit = 100;
P_nit100 = zeros(nx, ny); %reset P_nit100

%for nit = 100
P = zeros(nx, ny);
P(end, :) = y; %Dirichlet BC: P = y at x = 2
P(1, :) = 0; %Dirichlet BC: P = 0 at x = 0

for itt = 1:nit
    pd = P;
    for i = 2:nx-1
        for j = 2:ny-1
            P(i, j) = 0.25 * (pd(i+1, j) + pd(i-1, j) + pd(i, j+1) + pd(i, j-1));
        end
    end
    P(2:nx-1, 1) = P(2:nx-1, 2); % ∂P/∂y @ y = 0
    P(2:nx-1, ny) = P(2:nx-1, ny-1); % ∂P/∂y @ y = ymax
end
P_nit100 = P;

%test case (b): nit = 1000
nit = 1000;
P_nit1000 = zeros(nx, ny); %reset P_nit1000

%nit = 1000
P = zeros(nx, ny);
P(end, :) = y; %Dirichlet BC: P = y at x = 2
P(1, :) = 0; %Dirichlet BC: P = 0 at x = 0

for itt = 1:nit
    pd = P;
    for i = 2:nx-1
        for j = 2:ny-1
            P(i, j) = 0.25 * (pd(i+1, j) + pd(i-1, j) + pd(i, j+1) + pd(i, j-1));
        end
    end
    P(2:nx-1, 1) = P(2:nx-1, 2); % ∂P/∂y @ y = 0
    P(2:nx-1, ny) = P(2:nx-1, ny-1); % ∂P/∂y @ y = ymax
end
P_nit1000 = P;

y_indices = [1, 5, 10, 15, 20, 26]; %fixed values of j

%test case (a) nit = 100
for idx = y_indices
    figure
    plot(x, P_analytical(:, idx), '-')
    hold on
    plot(x, P_nit100(:, idx), 'o')
    title(['Comparison of Numerical and Analytical Solutions (nit = 100, j=', num2str(idx), ')'])
    xlabel('x')
    ylabel('P')
    legend(['Analytical, j=', num2str(idx)], ['Numerical, j=', num2str(idx)], 'Location', 'eastoutside')
    hold off
end

%test case (b) nit = 1000
for idx = y_indices
    figure
    plot(x, P_analytical(:, idx), '-')
    hold on
    plot(x, P_nit1000(:, idx), 'o')
    title(['Comparison of Numerical and Analytical Solutions (nit = 1000, j=', num2str(idx), ')'])
    xlabel('x')
    ylabel('P')
    legend(['Analytical, j=', num2str(idx)], ['Numerical, j=', num2str(idx)], 'Location', 'eastoutside')
    hold off
end
