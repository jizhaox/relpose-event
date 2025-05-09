function visualize_events(events)

figure(1); hold on
n_line = numel(events);
colors = lines(n_line);
for i = 1:n_line
    event_group = events{i};
    n_pt = numel(event_group);
    x = zeros(n_pt, 1);
    y = zeros(n_pt, 1);
    t = zeros(n_pt, 1);
    for j = 1:n_pt
        evt = event_group{j};
        img = evt.img;
        x(j) = img(1);
        y(j) = img(2);
        t(j) = evt.t;
    end
    plot3(x, y, t, 'o', 'Color', colors(i, :)); hold on
end
xlabel('x');
ylabel('y');
zlabel('time')