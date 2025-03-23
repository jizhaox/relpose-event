function [v_sol, R_l_sol, u_l_sol] = npt_event_solver_cop(data_orientations, data_events)
% N-point linear solver for event cameras developed in [1].
% It estimate linear velocity given camera orientations.
% author: Ji Zhao
% email: zhaoji84@gmail.com
% 2024/08/13
% Reference:
% [1] Ji Zhao, Banglei Guan, Zibin Liu, and Laurent Kneip.
%     Full-DoF Egomotion Estimation for Event Cameras Using Geometric Solvers.
%     IEEE/CVF Conference on Computer Vision and Pattern Recognition (CVPR), 2025.

SMALL_VALUE = 1e-4;

n = numel(data_orientations);
R_l_sol = cell(n, 1);
u_l_sol = cell(n, 1);
d_sol = cell(n, 1);
m_sol = cell(n, 1);
D = zeros(n, 3);
for i = 1:n
    orientations = data_orientations{i};
    events = data_events{i};

    d = line_dir(events, orientations);
    d_sol{i} = d;
    B = construct_coef_mat(events, orientations, d);
    [U,S,V] = svd(B);
    
    if(norm(V(4:6, end)) > SMALL_VALUE)
        x = V(:, end);
    else
        x = V(:, end) + V(:, end-1);
    end

    v = x(1:3);
    m = x(4:6);
    s = norm(m);
    v = v/s;
    m = m/s;
    m_sol{i} = m;

    e1_l = d;
    e2_l = -m;
    e3_l = cross(e1_l, e2_l);
    R_l = [e1_l, e2_l, e3_l];
    u_l = R_l'*v;
    R_l_sol{i} = R_l;
    u_l_sol{i} = u_l;
    D(i, :) = (u_l(2)*e3_l - u_l(3)*e2_l)';
end

[U,S,V] = svd(D);
v_sol = V(:, end);
v_sol = v_sol/norm(v_sol);

function d = line_dir(events, orientations)
n = numel(events);
B = zeros(n, 3);
for j = 1:n
    n_j = events{j}.plane_nml;
    R_j = orientations{j};
    n_j_p = R_j*n_j;
    B(j, :) = n_j_p';
end
[U,S,V] = svd(B);
d = V(:, end);
d = d/norm(d);

function B = construct_coef_mat(events, orientations, d)
n = numel(events);
B = zeros(n, 6);
for j = 1:n
    t_j = events{j}.t;
    tmp = [events{j}.img; 1];
    f_j = tmp/norm(tmp);
    R_j = orientations{j};
    f_j_p = R_j*f_j;
    B(j, :) = [t_j*cross(f_j_p, d); f_j_p]';
end
