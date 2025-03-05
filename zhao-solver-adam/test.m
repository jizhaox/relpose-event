clear;
num_line = 5; % number of lines
num_pt = 100; % number of points for each line
t_max = 0.5;
test_num = 100;
times = zeros(test_num, 1);
errors = zeros(test_num, 1);
w0 = [0;0;0]; % initialization
num_success = 0;

% hyperparameters for optimiziers
% Do not change them unless you know what you are doing
options = struct('Display', 'final','MaxFunEvals', 1e6,'MaxIter',1e6, 'TolFun', 1e-6, 'TolX', 1e-6, ...
    'stepSize',0.001,'beta1',0.9,'beta2',0.999,"epsilon",1e-9,"nEpochSize",3);
parames = struct('scale_lambda', 1e6, 'diff_var', 1e-6);

% help function
disp('run relpose_event() for help.');
%relpose_event();

format long
for i = 1:1
    [events, orientations, v_gt, w_gt] = generate_data(num_line, num_pt, t_max);
    w_gt, v_gt

    tic
    %% input parameters
    %[w, v] = relpose_event(events, w0, 303, parames, options)
    %% default parameters
    [w, v, objValue, t] = relpose_event(events, [], 101)
    toc
    
    err_w = norm(w - w_gt) / (norm(w) + norm(w_gt))
    err_v = rad2deg(acos(abs(v'*v_gt)/norm(v)/norm(v_gt)))
end
