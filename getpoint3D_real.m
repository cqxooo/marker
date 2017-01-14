function [X P1 P2] = getpoint3D_real(left, right, K, P1)
x1 = left.point;
x2 = right.point;
[R, t] = getRelative(K, K, x1, x2);

if nargin == 3
    P2 = [R t];
    P1 = [eye(3) zeros(3,1)];
else
    P2 = [R*P1(:,1:3) R*P1(:,end)+t];
end




    % set the options for lsqnonlin
%     options = optimset('lsqnonlin');
    options = optimset('Algorithm',{'levenberg-marquardt', 0.01});
    options = optimset(options, 'display', 'iter');
    options = optimset(options,'TolFun',1e-20);
    options = optimset(options,'TolX',1e-20);
    options = optimset(options,'MaxIter',100);
    options = optimset(options,'MaxFunEvals',800);
% set the options for lsqnonlin
[P2,resnorm,residual] = lsqnonlin(@(P2)costGold(x1, x2, P1, P2, K), P2,...
              [],[],options);
x1 = inv(K)*x1;
x1 = x1./repmat(x1(3,:),3,1);
x2 = inv(K)*x2;
x2 = x2./repmat(x2(3,:),3,1);
E = Tx(P2(:,end))*P2(:,1:3);
[u d v] = svd(E);
dd = (d(1,1)+d(2,2))/2;
d = diag([dd dd 0]);
E = u*d*v';
[R, t] = decomposeE_new(E, x1(:,1), x2(:,1));
P2 = [R t];
X.id = left.id;
X.point = triangulate(x1, x2, P1, P2);    



function cost = costGold(x1, x2, P1, P2, K)
x1 = inv(K)*x1;
x1 = x1./repmat(x1(3,:),3,1);
x2 = inv(K)*x2;
x2 = x2./repmat(x2(3,:),3,1);
X = triangulate(x1, x2, P1, P2);

ex1 = P1*[X;ones(1,size(X,2))];
ex1 = ex1./repmat(ex1(3,:),3,1);
ex2 = P2*[X;ones(1,size(X,2))];
ex2 = ex2./repmat(ex2(3,:),3,1);
cost = sum((x1(1:2,:)-ex1(1:2,:)).^2 + (x2(1:2,:)-ex2(1:2,:)).^2);
cost = sqrt(mean(cost));
return;



