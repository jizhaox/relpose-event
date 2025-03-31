function [triang_pt, pt1, pt2, lambda] = triangulate(p1, f1, p2, f2)
% triangulation for two rays (p1, f1) and (p2, f2)
% Reference:
% https://github.com/laurentkneip/opengv/blob/master/src/triangulation/methods.cpp#L66

t12 = p2 - p1;
A = [f1'*f1, -f1'*f2; f2'*f1, -f2'*f2];
b = [t12'*f1; t12'*f2];
lambda = A\b;
pt1 = p1 + lambda(1)*f1;
pt2 = p2 + lambda(2)*f2;
triang_pt = (pt1 + pt2) / 2;
