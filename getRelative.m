function [R, t] = getRelative(Kl, Kr, c1, c2)
%compute the essential matrix E
p1 = inv(Kl)*c1;
p1 = p1./repmat(p1(3,:),3,1);
p2 = inv(Kr)*c2;
p2 = p2./repmat(p2(3,:),3,1);
E = eightpoint(p1, p2);
%extract R and t
[R, t] = decomposeE_new(E, p1(:,1), p2(:,1));
