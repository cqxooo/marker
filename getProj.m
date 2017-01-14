function [P1 P2] = getProj(left, right, E, K, P1)
p1 = inv(K)*left.point;
p1 = p1./repmat(p1(3,:),3,1);
p2 = inv(K)*right.point;
p2 = p2./repmat(p2(3,:),3,1);
[R, t] = decomposeE_new(E, p1(:,1), p2(:,1));

if nargin == 4
    P2 = [R t];
    P1 = [eye(3) zeros(3,1)];
else
    P2 = [R*P1(:,1:3) R*P1(:,end)+t];
end