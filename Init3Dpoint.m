function cloud = Init3Dpoint(left, right, Proj, K)
id = left(1).id;
for imLoop=2:length(left)
    for i=1:length(left(imLoop).id)
       idx = find(id==left(imLoop).id(i));
       if isempty(idx)
          id = [id left(imLoop).id(i)]; 
       end
    end    
end
id = sort(id);
for i=1:length(id)
    tmp = [];
    for imLoop=1:length(left)
        idx = find(left(imLoop).id==id(i));
       if ~isempty(idx)
          tmp = [tmp;imLoop idx]; 
       end
    end
    if ~isempty(tmp)
        for j=1:size(tmp,1)
            x1(:,j) = left(tmp(j,1)).point(:,tmp(j,2)); 
            x2(:,j) = right(tmp(j,1)).point(:,tmp(j,2));
            P1(:,:,j) = Proj(:,:,tmp(j,1));
            P2(:,:,j) = Proj(:,:,tmp(j,1)+1);
        end
        x1 = inv(K)*x1;
        x1 = x1./repmat(x1(3,:),3,1);
        x2 = inv(K)*x2;
        x2 = x2./repmat(x2(3,:),3,1);
        Xw(:,i) = triangulate2(x1, x2, P1, P2);    
    end
end
cloud.id = id;
cloud.point = Xw;