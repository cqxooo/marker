function cloud = buildCloud(data, E, Proj, K)
for imLoop = 1:length(data)-1
    [data(imLoop) data(imLoop+1)] = findcorres_E(data(imLoop), data(imLoop+1), E(:,:,imLoop), K);  
end
for imLoop = 1:length(data)-1
    [left(imLoop) right(imLoop)] = findcorres(data(imLoop), data(imLoop+1));
end
save corres.mat
cloud = Init3Dpoint(left, right, Proj, K);
cloud.point
cloud = globalRefine1(cloud, data, Proj, K);
cloud.point