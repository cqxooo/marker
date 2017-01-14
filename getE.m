function E = getE(left, right, K)
p1 = inv(K)*left.point;
p1 = p1./repmat(p1(3,:),3,1);
p2 = inv(K)*right.point;
p2 = p2./repmat(p2(3,:),3,1);
E = eightpointE(p1, p2);
E = refineE(E, p1, p2);