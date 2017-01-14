function result = getpoint3D(imdata)
for imLoop = 1:length(imdata)-1
    Kl = imdata(imLoop).K;
    Kr = imdata(imLoop+1).K;
    x1 = imdata(imLoop).cens;
    x2 = imdata(imLoop+1).cens;
    [R, t] = getRelative(Kl, Kr, x1, x2);
    data.x1 = x1;
    data.x2 = x2;
    param0 = [Kl(1,1) Kl(2,2) Kl(1:2,3)';
              Kr(1,1) Kr(2,2) Kr(1:2,3)';
              R t];


    % set the options for lsqnonlin
%     options = optimset('lsqnonlin');
    options = optimset('Algorithm',{'levenberg-marquardt', 0.01});
    options = optimset(options, 'display', 'iter');
    options = optimset(options,'TolFun',1e-20);
    options = optimset(options,'TolX',1e-20);
    options = optimset(options,'MaxIter',30);
    options = optimset(options,'MaxFunEvals',800);
% set the options for lsqnonlin
[param,resnorm,residual] = lsqnonlin(@costGold, param0,...
              [],[],options,data);
Kl = [diag([param(1,1) param(1,2)]) param(1,3:4)';0 0 1];
Kr = [diag([param(2,1) param(2,2)]) param(2,3:4)';0 0 1];
R = param(3:end,1:3);
t = param(3:end,4);
P1 = [eye(3) zeros(3,1)];
P2 = [R t];
x1 = inv(Kl)*x1;
x1 = x1./repmat(x1(3,:),3,1);
x2 = inv(Kr)*x2;
x2 = x2./repmat(x2(3,:),3,1);
X = triangulate(x1, x2, P2);    
result(imLoop).id = imdata(imLoop).id;
result(imLoop).point = X;
end

function cost = costGold(param, data)
Kl = [diag([param(1,1) param(1,2)]) param(1,3:4)';0 0 1];
Kr = [diag([param(2,1) param(2,2)]) param(2,3:4)';0 0 1];
R = param(3:end,1:3);
t = param(3:end,4);
P1 = [eye(3) zeros(3,1)];
P2 = [R t];
x1 = data.x1;
x2 = data.x2;

x1 = inv(Kl)*x1;
x1 = x1./repmat(x1(3,:),3,1);
x2 = inv(Kr)*x2;
x2 = x2./repmat(x2(3,:),3,1);
X = triangulate(x1, x2, P2);

ex1 = P1*[X;ones(1,size(X,2))];
ex1 = ex1./repmat(ex1(3,:),3,1);
ex2 = P2*[X;ones(1,size(X,2))];
ex2 = ex2./repmat(ex2(3,:),3,1);
cost = sum((x1(:)-ex1(:)).^2 + (x2(:)-ex2(:)).^2);
cost = sqrt(cost)/size(x1,2);
return;



