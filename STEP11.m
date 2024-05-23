clc;
close all;
clear;

%Ghia et al. data for Re=100
Ghia_y_100 = [0.0000, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813, 0.4531, 0.5000, 0.6172, 0.7344, 0.8516, 0.9531, 0.9609, 0.9688, 0.9766, 1.0000];
Ghia_u_100 = [0.0000, -0.03717, -0.04192, -0.04775, -0.06434, -0.10150, -0.15662, -0.21090, -0.20581, -0.13641, 0.00332, 0.23151, 0.68717, 0.73722, 0.78871, 0.84123, 1.0000];
Ghia_x_100 = [0.0000, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563, 0.2266, 0.2344, 0.5000, 0.8047, 0.8594, 0.9063, 0.9453, 0.9531, 0.9609, 0.9688, 1.0000];
Ghia_v_100 = [0.0000, 0.09233, 0.10091, 0.10890, 0.12317, 0.16077, 0.17507, 0.17527, 0.05454, -0.24533, -0.22445, -0.16914, -0.10313, -0.08864, -0.07391, -0.05906, 0.0000];

%Ghia et al. data for Re=400
Ghia_y_400 = [0.0000, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813, 0.4531, 0.5000, 0.6172, 0.7344, 0.8516, 0.9531, 0.9609, 0.9688, 0.9766, 1.0000];
Ghia_u_400 = [0.0000, -0.08186, -0.09266, -0.10338, -0.14612, -0.24299, -0.32726, -0.17119, -0.11477, 0.02135, 0.16256, 0.29093, 0.55892, 0.61756, 0.68439, 0.75837, 1.0000];
Ghia_x_400 = [0.0000, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563, 0.2266, 0.2344, 0.5000, 0.8047, 0.8594, 0.9063, 0.9453, 0.9531, 0.9609, 0.9688, 1.0000];
Ghia_v_400 = [0.0000, 0.18360, 0.19713, 0.20920, 0.22965, 0.28124, 0.30203, 0.30174, 0.05186, -0.38598, -0.44993, -0.23827, -0.22847, -0.19254, -0.15663, -0.12146, 0.0000];

%load results
results_files = {'Results_case_1.mat', 'Results_case_2.mat', 'Results_case_3.mat'};
results = cell(1, length(results_files));
for i = 1:length(results_files)
    data = load(results_files{i});
    if isfield(data, 'u') && isfield(data, 'v') && isfield(data, 'p') && isfield(data, 'X') && isfield(data, 'Y')
        results{i} = struct('u', data.u, 'v', data.v, 'p', data.p, 'X', data.X, 'Y', data.Y);
    else
        error('Missing fields in the loaded data for test case %d', i);
    end
end

for i = 1:length(results)
    disp(['Test Case ', num2str(i)]);
    disp('u:');
    disp(results{i}.u(1:5, 1:5));
    disp('v:');
    disp(results{i}.v(1:5, 1:5));
    disp('p:');
    disp(results{i}.p(1:5, 1:5));
end

for i = 1:length(results)
    res = results{i};
    u = res.u;
    v = res.v;
    p = res.p;
    X = res.X;
    Y = res.Y;
    
    %velocity magnitude plot
    VEC = sqrt(u.^2 + v.^2);
    figure;
    contourf(X, Y, VEC', 'LineColor', 'none');
    colorbar;
    hold on;
    skip = 5;
    quiver(Y(1:skip:end), X(1:skip:end), u(1:skip:end), v(1:skip:end), 4, 'w');
    xlim([0 1]);
    ylim([0 1]);
    xlabel('x');
    ylabel('y');
    title(sprintf('Fluid velocity field (magnitude) for test case (%d)', i));
    
    %pressure field plot
    figure;
    contourf(X, Y, p', 'LineColor', 'none');
    colorbar;
    hold on;
    quiver(Y(1:skip:end), X(1:skip:end), u(1:skip:end), v(1:skip:end), 4, 'w');
    xlim([0 1]);
    ylim([0 1]);
    xlabel('x');
    ylabel('y');
    title(sprintf('Pressure field with velocity vectors for test case (%d)', i));
    
    %u-velocity profile along x=0.5
    figure;
    if mod(size(u, 1), 2) == 0
        plot(Y(:, size(u, 1)/2), u(:, size(u, 1)/2), 'k');
    else
        plot(Y(:, ceil(size(u, 1)/2)), u(:, ceil(size(u, 1)/2)), 'k');
    end
    xlabel('y');
    ylabel('u(y) at x=0.5');
    grid on;
    title(sprintf('u(y) - Comparison to Ghia et al. for test case (%d)', i));
    hold on;
    if i == 2
        plot(Ghia_y_100, Ghia_u_100, 'r*');
        legend('u(y)', 'Ghia solution (Re=100)', 'Location', 'northwest');
    elseif i == 3
        plot(Ghia_y_400, Ghia_u_400, 'r*');
        legend('u(y)', 'Ghia solution (Re=400)', 'Location', 'northwest');
    else
        legend('u(y)', 'Location', 'northwest');
    end
    
    %v-velocity profile along y=0.5
    figure;
    if mod(size(v, 2), 2) == 0
        plot(X(size(v, 2)/2, :), v(size(v, 2)/2, :), 'k');
    else
        plot(X(ceil(size(v, 2)/2), :), v(ceil(size(v, 2)/2), :), 'k');
    end
    xlabel('x');
    ylabel('v(x) at y=0.5');
    grid on;
    title(sprintf('v(x) - Comparison to Ghia et al. for test case (%d)', i));
    hold on;
    if i == 2
        plot(Ghia_x_100, Ghia_v_100, 'r*');
        legend('v(x)', 'Ghia solution (Re=100)', 'Location', 'northwest');
    elseif i == 3
        plot(Ghia_x_400, Ghia_v_400, 'r*');
        legend('v(x)', 'Ghia solution (Re=400)', 'Location', 'northwest');
    else
        legend('v(x)', 'Location', 'northwest');
    end
end

