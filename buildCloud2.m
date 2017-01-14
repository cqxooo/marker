function cloud = buildCloud2(data, Proj, K)
for imLoop =1:length(data)-1
    [lcor rcor] = findcorres(data(imLoop), data(imLoop+1));
    x1 = inv(K)*lcor.point;
    x1 = x1./repmat(x1(3,:),3,1);
    x2 = inv(K)*rcor.point;
    x2 = x2./repmat(x2(3,:),3,1);
    E = Tx(Proj(:,end,imLoop+1))*Proj(:,1:3,imLoop+1);
    [u d v] = svd(E);
    dd = (d(1,1)+d(2,2))/2;
    d = diag([dd dd 0]);
    E = u*d*v';
    [R, t] = decomposeE_new(E, x1(:,1), x2(:,1));
    Proj(:,:,imLoop+1) = [R t];
    if imLoop == 1
        cloud.id = lcor.id;
        cloud.point = triangulate(x1, x2, Proj(:,:,imLoop), Proj(:,:,imLoop+1));  
    else
        tmp.id = lcor.id;
        tmp.point = triangulate(x1, x2, Proj(:,:,imLoop), Proj(:,:,imLoop+1));  
        cloud = updateCloud(cloud, tmp);
    end
    
end