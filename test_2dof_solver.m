% N-point linear solver for event cameras [1]
% author: Ji Zhao
% email: zhaoji84@gmail.com
% 2024/05/30
% Reference:
% [1] Ji Zhao, Banglei Guan, Zibin Liu, and Laurent Kneip.
%     Full-DoF Egomotion Estimation for Event Cameras Using Geometric Solvers.
%     IEEE/CVF Conference on Computer Vision and Pattern Recognition (CVPR), 2025.

%% parameters for scene generation
num_line = 5; % number of lines
num_pt = 20; % number of points for each line
t_max = 0.5;
[events, orientations, v_gt, w_gt]  = generate_data(num_line, num_pt, t_max);

%% help function
disp('run npt_event_solver_cop() for help.');

%% linear velocity estimation
[v_sol, line_struct_all, tm] = npt_event_solver_cop(orientations, events);

%% performance evaluation
err_v = evaluate_lin_error(v_sol, v_gt);
% angular re-projection error
err_reproj = evaluate_angular_reproj_error(events, w_gt, v_sol, line_struct_all);

%% output
format long
disp('N-point linear solver for event cameras:')
disp('linear velocity: ground truth and estimation (with scale ambiguity)');
disp([v_gt, v_sol])
disp('error of linear velocity using ground truth (unit: degree)')
disp(err_v)
disp('error of linear velocity using angular re-projection (unit: degree)')
disp(max(err_reproj))
disp('runtime (unit: microsecond):');
disp(tm)

