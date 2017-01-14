function E = eightpointE(p1, p2)
x1 = p1';
x2 = p2';

A = [x2(:,1).*x1(:,1), x2(:,1).*x1(:,2), x2(:,1),...
     x2(:,2).*x1(:,1), x2(:,2).*x1(:,2), x2(:,2),...
     x1(:,1), x1(:,2), ones(length(x1),1)];
[u d v] = svd(A);
v=v(:,end);
v = v/norm(v);
E = reshape(v,3,3)';
[u d v] = svd(E);
dd = (d(1,1)+d(2,2))/2;
d = diag([dd dd 0]);
E = u*d*v';



 