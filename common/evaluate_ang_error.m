function e = evaluate_ang_error(a, b)

a = a(:);
b = b(:);
e = norm(a-b)/(norm(a) + norm(b));
