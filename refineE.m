function E = refineE(estE, p1, p2)
%extract R and t
[R, t] = decomposeE_new(estE, p1(:,1), p2(:,1));
tmp = vrrotmat2vec(R)';
mot = [tmp(1:3)*tmp(4);t];

%     options = optimset('Algorithm',{'levenberg-marquardt', 0.01});
    options = optimset('lsqnonlin');
    options = optimset(options, 'display', 'iter');
    options = optimset(options,'TolFun',1e-20);
    options = optimset(options,'TolX',1e-20);
    options = optimset(options,'MaxIter',500);
    options = optimset(options,'MaxFunEvals',10000);
    [param,resnorm,residual] = lsqnonlin(@(mot)refineAlgErr(mot, p1, p2),...
    mot,[],[],options);
v = param(1:3);
R = vrrotvec2mat([v/norm(v);norm(v)]);
t = param(4:6);
E = Tx(t)*R;
% [u d v] = svd(E);
% dd = (d(1,1)+d(2,2))/2;
% d = diag([dd dd 0]);
% E = u*d*v';

end

function err = refineAlgErr(mot, p1, p2)
v = mot(1:3);
R = vrrotvec2mat([v/norm(v);norm(v)]);
t = mot(4:6);
E = Tx(t)*R;
% [u d v] = svd(E);
% dd = (d(1,1)+d(2,2))/2;
% d = diag([dd dd 0]);
% E = u*d*v';
X2tEX1 = zeros(1,size(p1,2));
for n = 1:size(p1,2)
    X2tEX1(n) = p2(:,n)'*E*p1(:,n);
end
EX1 = E*p1;
EtX2 = E'*p2;   
% Evaluate distances
err =  X2tEX1.^2 ./ (EX1(1,:).^2 + EX1(2,:).^2 + EtX2(1,:).^2 + EtX2(2,:).^2);
% err = 0;
% for i=1:size(EX1,2)
%     l = EX1(:,i);
%     l1 = -sign(l(3))*l/norm(l(1:2));
%     l = EtX2(:,i); 
%     l2 = -sign(l(3))*l/norm(l(1:2));
%     err = err+abs(p2(:,i)'*l1)+abs(p1(:,i)'*l2);    
% end
end