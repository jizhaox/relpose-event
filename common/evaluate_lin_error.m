function theta = evaluate_lin_error(a, b)

cos_theta = abs(dot(a, b)) / (norm(a) * norm(b));
cos_theta = max(-1, min(1, cos_theta));
theta_rad = acos(cos_theta);
theta = rad2deg(theta_rad);
