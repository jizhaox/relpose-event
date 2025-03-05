function [data_events, data_orientations, v_gt, w_gt] = generate_data(n_line, n_pt, t_max)
% author: Ji Zhao
% email: zhaoji84@gmail.com
% 2024/05/30

if nargin < 1
    n_line = 2;
end
if nargin < 2
    n_pt = 5;
end
if nargin < 3
    t_max = 0.5;
end
if n_line < 2
    warning('At leaset 2 lines should be synthesized!')
end
if n_pt < 5
    warning('At leaset 5 points should be synthesized!')
end

%  linear and angular velocities
v_gt = (rand(3, 1) - 0.5)*10;
w_gt = (rand(3, 1) - 0.5)*0.25;

data_orientations = cell(n_line, 1);
data_events = cell(n_line, 1);
for i = 1:n_line
    [orientations, events] = generate_line_event(v_gt, w_gt, t_max, n_pt);
    data_orientations{i} = orientations;
    data_events{i} = events;
end

% normalize linear velocity since its scale is unobservable
tmp = norm(v_gt);
if tmp > 1e-10
    v_gt = v_gt / tmp;
end
