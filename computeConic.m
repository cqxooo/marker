function C = computeConic(p, T)
t = T*[p';ones(1,size(p,1))];
t = t./repmat(t(3,:),3,1);
AE = [t(1,:).*t(1,:);2*t(1,:).*t(2,:);t(2,:).*t(2,:);
      2*t(1,:);2*t(2,:);ones(1,size(p,1))]';
if (size(AE, 1)<5)
    disp('not enough points for ellipse fitting!');
    C = [];
    return;
end

[u, d, v] = svd(AE);
C = [v(1,6) v(2,6) v(4,6);v(2,6) v(3,6) v(5,6);v(4,6) v(5,6) v(6,6)];
C = C./C(3,3);
[V,D]=eig(C);
D = [D(1,1) D(2,2) D(3,3)];
D = sort(D);
if D(1)<0 & D(2) >0 & D(3)>0
    C = C;
else
    C = -C;
end
C = T'*C*T;
C = C/C(3,3);