function cloud = refine3D(data, K)
for imLoop = 1:length(data)-1
   [lcor rcor] = findcorres(data(imLoop), data(imLoop+1));
%    [E, lcor, rcor] = ransacE(lcor, rcor, K); 
   if imLoop == 1
       [X Proj(:,:,imLoop) Proj(:,:,imLoop+1)] = getpoint3D_real(lcor, rcor, K);
       tmp.id = X.id;
       tmp.point = X.point;
   else
%        [x Xw] = findcorres(data(imLoop+1), tmp);
%        P = compute_p( x.point, Xw.point, K );
       [X Proj(:,:,imLoop) Proj(:,:,imLoop+1)] = getpoint3D_real(lcor, rcor, K, Proj(:,:,imLoop));
       tmp = updateCloud(tmp, X);
   end
%    Xcor(imLoop) = X;
end
% cloud = buildCloud(Xcor);
cloud = tmp;
cloud = globalRefine(cloud, data, Proj, K);
