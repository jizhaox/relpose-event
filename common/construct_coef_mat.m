function A = construct_coef_mat(events, orientations)
% author: Ji Zhao
% email: zhaoji84@gmail.com
% 2024/05/30

n = numel(events);
A = zeros(n, 6);
for j = 1:n
    t_j = events{j}.t;
    tmp = [events{j}.img; 1];
    f_j = tmp/norm(tmp);
    R_j = orientations{j};
    f_j_p = R_j*f_j;
    A(j, :) = [t_j*f_j_p; f_j_p]';
end