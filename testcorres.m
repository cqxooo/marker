function testcorres
load corres.mat
imgFilename = {};
% for i=1:7
%     eval(['imgFilename{end+1} =''measureImages/high/raw' int2str(i)  '.jpg'';']);
% end
for i=1:3
    eval(['imgFilename{end+1} =''test/raw' int2str(i)  '.jpg'';']);
end
for imLoop = 1:length(left)
    fig1 = figure(imLoop);
   imshow(imread(imgFilename{imLoop}))
   hold on
   for i=1:length(left(imLoop).id)
       x = inv(K)*right(imLoop).point(:,i);
       plot(left(imLoop).point(1,i),left(imLoop).point(2,i),'MarkerEdgeColor', MyPalette(i), 'Marker', '.', 'MarkerSize', 14); 
%        drawLine(inv(K)'*E(:,:,imLoop)'*x,'r-');
   end
   fig2 = figure(imLoop+1);
   imshow(imread(imgFilename{imLoop+1}))
   hold on
   for i=1:length(right(imLoop).id)
       x = inv(K)*left(imLoop).point(:,i);
       plot(right(imLoop).point(1,i),right(imLoop).point(2,i),'MarkerEdgeColor', MyPalette(i), 'Marker', '.', 'MarkerSize', 14); 
%        drawLine(inv(K)'*E(:,:,imLoop)*x,'r-');
   end
   close(fig1)
   close(fig2)
end
% for imLoop = 1:length(left)
%         x1 = inv(K)*left(imLoop).point;
%         x1 = x1./repmat(x1(3,:),3,1);
%         x2 = inv(K)*right(imLoop).point;
%         x2 = x2./repmat(x2(3,:),3,1);
%         [R, t] = decomposeE_new(E(:,:,imLoop), x1(:,1), x2(:,1));
%         P1 = [eye(3) zeros(3,1)];
%         P2 = [R t];
%         X = triangulate(x1(:,1), x2(:,1), P1, P2);  
%         for i=1:size(X,2)
%             xx1 = P1*[X(:,i);1];
%             xx2 = P2*[X(:,i);1];
%             xx1 = xx1/xx1(3);
%             xx2 = xx2/xx2(3);
%         end
% end