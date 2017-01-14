function Depth = getDepth(P, X)
xi = P*X;
w = xi(3);
T = X(end,:);
m3n = sqrt(sum(P(3,1:3).*P(3,1:3)));
Depth = (sign(det(P(:,1:3)))*w)/(T*m3n);