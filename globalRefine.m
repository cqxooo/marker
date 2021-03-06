function cloud = globalRefine(Xw, data, Proj, K)
mot = [];
for i=2:size(Proj,3)
    tmp = vrrotmat2vec(Proj(:,1:3,i))';
    mot = [mot tmp(1:3)*tmp(4) Proj(:,end,i)];
end
pnum = size(mot,2);
    % set the options for lsqnonlin
%     options = optimset('lsqnonlin');
    options = optimset('Algorithm',{'levenberg-marquardt', 0.01});
    options = optimset(options, 'display', 'iter');
    options = optimset(options,'TolFun',1e-20);
    options = optimset(options,'TolX',1e-20);
    options = optimset(options,'MaxIter',100);
    options = optimset(options,'MaxFunEvals',10000);
% set the options for lsqnonlin
[param,resnorm,residual] = lsqnonlin(@(x)costfun(Xw.id, data, K, pnum, x),...
    [mot(:);Xw.point(:)],[],[],options);
%reconstruct from refined 3D point
% cloud.id = Xw.id;
% cloud.point = reshape(param(pnum*3+1:end),3,length(Xw.id)); 
%============================================================
%reconstruct from refined projection matrix
Proj(:,:,1) = [eye(3) zeros(3,1)];
mot = reshape(param(1:pnum*3),3,pnum);
for i=1:pnum/2
   v = mot(:,2*i-1);   
   Proj(:,:,i+1) = [vrrotvec2mat([v/norm(v);norm(v)]) mot(:,2*i)];    
end
cloud = get3D(data, Proj, K);
%============================================================
end




function err = costfun(id, data, K, pnum, x)
Proj(:,:,1) = [eye(3) zeros(3,1)];
mot = reshape(x(1:pnum*3),3,pnum);
Xw.id = id;
Xw.point = reshape(x(pnum*3+1:end),3,length(id));
for i=1:pnum/2
   v = mot(:,2*i-1);   
   Proj(:,:,i+1) = [vrrotvec2mat([v/norm(v);norm(v)]) mot(:,2*i)];    
end
err = 0;
for imLoop=1:length(data)
    [xcor Xcor] = findcorres(data(imLoop), Xw);
    ex = K*Proj(:,:,imLoop)*[Xcor.point;ones(1,size(Xcor.point,2))];
    ex = ex./repmat(ex(3,:),3,1);
    cost = sum((xcor.point(1:2,:)-ex(1:2,:)).^2);
    cost = sqrt(mean(cost));
    err = err + cost;
end
end