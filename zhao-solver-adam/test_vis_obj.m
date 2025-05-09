close all
clf;

sd = 40097;
method_type = 101;

%%
rng(sd);
num_line_total = 10; % number of lines
num_pt = 100; % number of points for each line
t_max = 0.5;
[events_total, orientations, v_gt, w_gt] = generate_data(num_line_total, num_pt, t_max);

x = linspace(-0.2, 0.2, 1001);
y = linspace(-0.2, 0.2, 1001);
z = [-0.2, -0.1, 0, 0.1, 0.2];
nx = numel(x);
ny = numel(y);
nz = numel(z);
[X,Y] = meshgrid(x, y);

%%
for num_line = 1:5
    num_line
    events = events_total(1:num_line);
    tic
    slice_set = evaluate_obj_batch(events, x, y, z, method_type);
    toc
    
    % visualization
    figure(1); clf
    for k = 1:nz
        z_k = z(k);
        data = log10(slice_set(:, :, k));
        Z = z_k * ones(size(X));
        % Display image with colormap
        hold on;
        h = surf(X, Y, Z, 'CData', data, 'FaceColor', 'texturemap', 'EdgeColor', 'none');
        colormap(parula);
       % view(3);
    end
    axis on
    colorbar;
    view(3);
    xlabel('$\omega_{x}$ (m/s)', 'Interpreter', 'latex');
    ylabel('$\omega_{y}$ (m/s)', 'Interpreter', 'latex');
    zlabel('$\omega_{z}$ (m/s)', 'Interpreter', 'latex');
    grid on;
    %axis equal;
    % title('Landscape of the Objective Function (log10 scale)');
    hold off;
    set(gca, 'XTick', z);
    set(gca, 'YTick', z);
    set(gca, 'ZTick', z);
    view(-25.6851, 15.586);
    
    filename = sprintf('landscape_%d_%d_%d.png', sd, method_type, num_line);
    exportgraphics(gcf, filename);
end