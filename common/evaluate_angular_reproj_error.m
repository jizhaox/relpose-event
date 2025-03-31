function err = evaluate_angular_reproj_error(events_all, w, v, line_struct_all)
% angular reprojection error
% Reference:
% https://github.com/laurentkneip/opengv/blob/master/src/sac_problems/relative_pose/CentralRelativePoseSacProblem.cpp#L185

% determine the scale factor between line moment and translation
n_line = numel(events_all);
err_all = zeros(n_line, 1);
s_all = zeros(n_line, 1);
for i = 1:n_line
    events= events_all{i};
    n_pt = numel(events);
    % extract line structure
    d = line_struct_all{i}(:, 1);
    m = -line_struct_all{i}(:, 2);
    a_coef = zeros(n_pt, 1);
    b_coef = zeros(n_pt, 1);
    for j = 1:n_pt
        tmp = [events{j}.img; 1];
        f_j = tmp/norm(tmp);

        dt = events{j}.t;
        R_j = expmap(w * dt);
        t_j = v * dt;

        f_j_p = R_j * f_j;
        m_j_p = cross(t_j, f_j_p);

        a_coef(j) = m'*f_j_p;
        b_coef(j) = d'*m_j_p;
    end
    s = -(a_coef'*b_coef) / (a_coef'*a_coef);
    err_all(i) = max(abs(a_coef*s+b_coef));
    s_all(i) = s;
end

err = [];
k = 0;
for i = 1:n_line
    events= events_all{i};
    n_pt = numel(events);
    % extract line structure
    s = s_all(i);
    d = line_struct_all{i}(:, 1);
    m = -line_struct_all{i}(:, 2) * s;
    pt = cross(d, m) / norm(d);
    for j = 1:n_pt
        tmp = [events{j}.img; 1];
        f_j = tmp/norm(tmp);

        dt = events{j}.t;
        R_j = expmap(w * dt);
        t_j = v * dt;
        f_j_p = R_j * f_j;
        
        triang_pt = triangulate(pt, d, t_j, f_j_p);
        theta1 = evaluate_lin_error(triang_pt-pt, d);
        theta2 = evaluate_lin_error(triang_pt-t_j, f_j_p);
        k = k + 1;
        err(k) = theta1 + theta2;
    end
end
err = err(:);
