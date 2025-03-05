function [orientations, events] = generate_line_event(v_gt, w_gt, t_max, n, is_display)
% author: Ji Zhao
% email: zhaoji84@gmail.com
% 2024/05/30

if nargin < 5
    is_display = false;
end

% choose n random time stamps from [-t_max/2, t_max/2]
t_all = t_max * sort(rand(n, 1)-0.5);

% construct 3D line
P = [0; 0; 1] + rand(3, 1)*5;
d = rand(3, 1) - 0.5;
d = d / norm(d);

% generate camera poses, sampled points, and events
translations = cell(n, 1);
orientations = cell(n, 1);
events = cell(n, 1);
for j = 1:n
    dt = t_all(j);
    R_j = expmap(w_gt * dt);
    t_j = v_gt * dt;
    orientations{j} = R_j;
    translations{j} = t_j;

    len = (rand(1) - 0.5)*5;
    pt = P + d*len;
    pt_new = R_j'*(pt - t_j);
    img = pt_new(1:2)/pt_new(3);

    len2 = len + 5;
    pt2 = P + d*len2;
    pt_new2 = R_j'*(pt2 - t_j);
    img2 = pt_new2(1:2)/pt_new2(3);
    
    line_dir = [img-img2; 0];
    nml_flow = [-line_dir(2); line_dir(1)];
    plane_nml = cross(pt_new, line_dir);
    plane_nml = plane_nml/norm(plane_nml);
%     plane_nml = cross(pt_new, pt_new2); % another method

    events{j} = struct('t', t_all(j), ...
        'img', img, 'nml_flow', nml_flow, ...
        'plane_nml', plane_nml);
end

if is_display
    % check the incidence relationship (Eq. (2))
    m = cross(P, d);
    [err_max] = check_incidence(d, m, events, orientations, translations);
    disp(['maximal residual is ' num2str(err_max)]);
    % check the incidence relationship (Eq. (4)(7))
    [R_l, u_l] = re_frame(d, P, [0;0;0], v_gt);
    [err_max2] = check_incidence2(R_l, u_l, events, orientations);
    disp(['maximal residual is ' num2str(err_max2)]);
    % check surface normals: method 1
    err = check_line_normal(events, orientations);
    disp(['minimal singular value of surface normal matrix is ' num2str(err)]);
    % check surface normals: method 2
    err2 = check_line_normal2(events, orientations, d);
    disp(['minimal singular value of surface normal matrix is ' num2str(err2)]);
end

end

function err = check_line_normal(events, orientations)
    n = numel(events);
    line_nml_mat = zeros(n, 3);
    for j = 1:n
        plane_nml = events{j}.plane_nml;
        R_j = orientations{j};
        line_nml_p = R_j * plane_nml;
        line_nml_mat(j, :) = line_nml_p';
    end
    S = svd(line_nml_mat);
    err = min(S(:));
end

function err = check_line_normal2(events, orientations, d)
    n = numel(events);
    line_nml_mat = zeros(n, 3);
    for j = 1:n
        tmp = [events{j}.img; 1];
        f_j = tmp/norm(tmp);
        R_j = orientations{j};
        line_nml_p = cross(R_j * f_j, d);
        line_nml_mat(j, :) = line_nml_p';
    end
    S = svd(line_nml_mat);
    err = min(S(:));
end

function [R_l, u_l] = re_frame(d, P, T0, v)
    % check Eq.(4)
    e1_l = d;
    P_p = point_line_perp_foot(T0, P, d);
    e3_l = T0 - P_p;
    s = 1/norm(e3_l);
    e3_l = s*e3_l;
    e2_l = -cross(e1_l, e3_l);
    R_l = [e1_l, e2_l, e3_l];
    % apply Eq.(3). 
    % The scale of camera velocity should be consistent with the line
    u_l = R_l'*(s*v);
end

function P_p = point_line_perp_foot(P0, P1, V1)
    P_p = P1 + dot(P0-P1,V1)*V1;
end

function [err_max, err] = check_incidence(d, m, events, orientations, translations)
    n = numel(events);
    err = zeros(n, 1);
    for j = 1:n
        tmp = [events{j}.img; 1];
        f_j = tmp/norm(tmp);
        R_j = orientations{j};
        t_j = translations{j};
        f_j_p = R_j * f_j;
        m_j_p = cross(t_j, f_j_p);
        err(j) = d'*m_j_p + m'*f_j_p;
    end
    err_max = max(abs(err(:)));
end

function [err_max, err, err2] = check_incidence2(R_l, u_l, events, orientations)
    e1_l = R_l(:, 1);
    e2_l = R_l(:, 2);
    e3_l = R_l(:, 3);
    uy_l = u_l(2);
    uz_l = u_l(3);
    n = numel(events);
    err = zeros(n, 1);
    for j = 1:n
        dt = events{j}.t;
        tmp = [events{j}.img; 1];
        f_j = tmp/norm(tmp);
        R_j = orientations{j};
        f_j_p = R_j * f_j;
        m_j_p = cross(dt*R_l*u_l, f_j_p);
        err(j) = e1_l'*m_j_p - f_j_p'*e2_l;
    end
    A = construct_coef_mat(events, orientations);
    x = [uz_l*e2_l - uy_l*e3_l; e2_l];
    err2 = A*x;

    err_max = max(abs([err; err2]));
end