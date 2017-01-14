function [left right] = findcorres_E(left, right, E, K)
eps = 10^(-7);
x1 = inv(K)*left.point;
x1 = x1./repmat(x1(3,:),3,1);
x2 = inv(K)*right.point;
x2 = x2./repmat(x2(3,:),3,1);
countId = max(left.id)+1;
[R, t] = decomposeE_new(E, x1(:,1), x2(:,1));
P1 = [eye(3) zeros(3,1)];
P2 = [R t];
for i=1:length(left.id)
    d = [];
    for j=1:length(right.id)
        X = triangulate(x1(:,i), x2(:,j), P1, P2);  
        xx1 = P1*[X;1];
        xx2 = P2*[X;1];
        xx1 = xx1/xx1(3);
        xx2 = xx2/xx2(3);
        dd = computeSampson(x1(:,i), x2(:,j),E)+distPnts(xx1, x1(:,i))+distPnts(xx2, x2(:,j));
        d = [d dd];
    end
    [~,idx] = min(d);
    while left.id(i)~=0 && right.id(idx)~=0
        d(idx) = 1000;
        [~,idx] = min(d);
    end
    if left.id(i)==0 && right.id(idx)==0
        left.id(i) = countId;
        right.id(idx) = countId;
        countId = countId + 1;
    elseif left.id(i)~=0 && right.id(idx)==0
        right.id(idx) = left.id(i);        
    end
end
end

function d = computeSampson(x1,x2,E)
EX1 = E*x1;
EtX2 = E'*x2;  
X2tEX1 = x2'*E*x1;
d =  X2tEX1^2/(EX1(1)^2 + EX1(2)^2 + EtX2(1)^2 + EtX2(2)^2);
end
