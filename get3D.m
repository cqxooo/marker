function cloud = get3D(data, Proj, K)
for imLoop = 1:length(data)-1
    [left(imLoop) right(imLoop)] = findcorres(data(imLoop), data(imLoop+1));
    x1 = inv(K)*left(imLoop).point;
    x1 = x1./repmat(x1(3,:),3,1);
    x2 = inv(K)*right(imLoop).point;
    x2 = x2./repmat(x2(3,:),3,1);
    E = Tx(Proj(:,end,imLoop+1))*Proj(:,1:3,imLoop+1);
    [u d v] = svd(E);
    dd = (d(1,1)+d(2,2))/2;
    d = diag([dd dd 0]);
    E = u*d*v';
    [R, t] = decomposeE_new(E, x1(:,1), x2(:,1));
    Proj(:,:,imLoop+1) = [R t];
end
cloud = Init3Dpoint(left, right, Proj, K);