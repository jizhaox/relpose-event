clear;

%% parameters for scene generation
num_line = 5; % number of lines
num_pt = 100; % number of points for each line
t_max = 0.5;

%% help function
disp('run relpose_event() for help.');

%% method type
% The recommended method types are 103 and 303.
% 101: incidence + exact
% 102: incidence + approximation
% 103: incidence + cascade
% 301: coplanarity + exact
% 302: coplanarity + approximation
% 303: coplanarity + cascade
method_type = 103;

%% scene generation & egomotion estimation
w0 = []; % intialization
[events, orientations, v_gt, w_gt] = generate_data(num_line, num_pt, t_max);
[w_est, v_est] = relpose_event(events, w0, method_type);

%% performance evaluation
err_w = norm(w_est - w_gt) / (norm(w_est) + norm(w_gt));
err_v = rad2deg(acos(abs(v_est'*v_gt)/norm(v_est)/norm(v_gt)));
% output
format long
disp('angular velocity: ground truth and estimation');
disp([w_gt, w_est])
disp('linear velocity: ground truth and estimation (with scale ambiguity)');
disp([v_gt, v_est])
disp('error for angular velocity')
disp(err_w)
disp('error for linear velocity (unit: degree)')
disp(err_v)
