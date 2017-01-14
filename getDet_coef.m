function [fdet, r, s] = getDet_coef(v1, v2, v3)
A = [v1;v2];
AtA = A*A';
c = A*v3';
a1 = AtA(1,1); a2 = AtA(1,2); b1 = AtA(2,1); b2 = AtA(2,2);
c1 = c(1); c2 = c(2);

detAtA = det(AtA);
fdet=abs(detAtA);
r=(c1*b2-c2*b1)*detAtA;
s=(a1*c2-a2*c1)*detAtA;