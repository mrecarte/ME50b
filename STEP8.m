clc;
close all;
clear;

vis = 0.1;
xmax = 2;
ymax = 2;
nx = 21;
ny = 21;
nt = 51;
tmax = 0.5;

dx = xmax / (nx - 1);
dy = ymax / (ny - 1);
dt = tmax / (nt - 1);

x = linspace(0, xmax, nx);
y = linspace(0, ymax, ny);
u = ones(ny, nx);
v = ones(ny, nx);

%init conditions
for i = 1:nx
    for j = 1:ny
        if x(i) >= 0.5 && x(i) <= 1 && y(j) >= 0.5 && y(j) <= 1
            u(j, i) = 2;
            v(j, i) = 2;
        else
            u(j, i) = 1;
            v(j, i) = 1;
        end
    end
end

times_to_plot = [1, round(0.25 / dt) + 1, nt]; % t=0, t=0.25, t=0.5

for n = 1:nt
    un = u;
    vn = v;
    for i = 2:nx-1
        for j = 2:ny-1
            u(j, i) = un(j, i) - un(j, i) * (dt/dx) * (un(j, i) - un(j, i-1)) - un(j, i) * (dt/dy) * (un(j, i) - un(j-1, i)) ...
                      + vis * (dt/dx^2) * (un(j, i+1) - 2*un(j, i) + un(j, i-1)) + vis * (dt/dy^2) * (un(j+1, i) - 2*un(j, i) + un(j-1, i));
        end
    end

    %ENFORCEBOUNDARY CONDITIONS
    u(:, 1) = 1; u(:, end) = 1; u(1, :) = 1; u(end, :) = 1;

    if ismember(n, times_to_plot)
        figure;
        surf(x, y, u', 'EdgeColor', 'none');
        title(sprintf('2D Viscous Burgers - u-field at Time Step %d', n));
        xlabel('X');
        ylabel('Y');
        zlabel('U field');
        drawnow;
    end
end
