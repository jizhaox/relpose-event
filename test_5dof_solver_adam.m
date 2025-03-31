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
[w_est, v_est, line_struct_all, obj, tm] = relpose_event(events, w0, method_type);

%% performance evaluation
err_w = evaluate_ang_error(w_est, w_gt);
err_v = evaluate_lin_error(v_est, v_gt);
% angular re-projection error
err_reproj = evaluate_angular_reproj_error(events, w_est, v_est, line_struct_all);

%% output
format long
disp('angular velocity: ground truth and estimation');
disp([w_gt, w_est])
disp('linear velocity: ground truth and estimation (with scale ambiguity)');
disp([v_gt, v_est])
disp('error of angular velocity using ground truth')
disp(err_w)
disp('error of linear velocity using ground truth (unit: degree)')
disp(err_v)
disp('error of linear velocity using angular re-projection (unit: degree)')
disp(max(err_reproj))
disp('runtime of angular and linear velocities estimation (unit: microsecond):');
disp(tm)
