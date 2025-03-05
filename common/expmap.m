function R = expmap(omega)
% exponential map

theta = norm(omega); % rotation angle
if theta < 1e-10
    R = eye(3);
    return;
end
% skew-symmetry matrix K
K = [0 -omega(3) omega(2); 
     omega(3) 0 -omega(1); 
     -omega(2) omega(1) 0];

% Rodrigues formula
R = eye(3) + (sin(theta) / theta) * K + ((1 - cos(theta)) / (theta^2)) * K^2;
