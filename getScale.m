function scale = getScale(imf, Proj, K)
[pathstr,imNameC,~] = fileparts(imf{1}); 
load ([pathstr '/basic.mat'])
for i=1:2
    for imLoop=1:length(basic)-1
        x1(:,imLoop) = basic(imLoop).point(:,i); 
        x2(:,imLoop) = basic(imLoop+1).point(:,i);
        P1(:,:,imLoop) = Proj(:,:,imLoop);
        P2(:,:,imLoop) = Proj(:,:,imLoop+1);
    end
    x1 = inv(K)*x1;
    x1 = x1./repmat(x1(3,:),3,1);
    x2 = inv(K)*x2;
    x2 = x2./repmat(x2(3,:),3,1);
    Xw(:,i) = triangulate2(x1, x2, P1, P2);    
end
scale = 105.02/distPnts(Xw(:,1),Xw(:,2));
