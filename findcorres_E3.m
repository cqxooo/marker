function [left right] = findcorres_E3(left, right, E, K)
eps = 10^(-7);
x1 = inv(K)*left.point;
x1 = x1./repmat(x1(3,:),3,1);
x2 = inv(K)*right.point;
x2 = x2./repmat(x2(3,:),3,1);
countId = max(left.id)+1;
[R, t] = decomposeE_new(E, x1(:,1), x2(:,1));
P1 = [eye(3) zeros(3,1)];
P2 = [R t];
[et e] = getEpipole(E);
et = et/et(3);
e = e/e(3);
for i=1:length(left.id)
    d = [];
    for j=1:length(right.id)
        X = triangulate(x1(:,i), x2(:,j), P1, P2);  
        xx1 = P1*[X;1];
        xx2 = P2*[X;1];
        xx1 = xx1/xx1(3);
        xx2 = xx2/xx2(3);
        d = [d distPnts(xx1, x1(:,i))+distPnts(xx2, x2(:,j))];
    end
    [mind,midx] = min(d);
    idx2 = find(d/mind<10);
    if ~isempty(idx2)
        for k=1:length(idx2)
            d = [];
           for m=1:length(left.id)
               X = triangulate(x1(:,m), x2(:,idx2(k)), P1, P2);  
               xx1 = P1*[X;1];
               xx2 = P2*[X;1];
               xx1 = xx1/xx1(3);
               xx2 = xx2/xx2(3);
               d = [d distPnts(xx1, x1(:,m))+distPnts(xx2, x2(:,idx2(k)))];
           end
           [mind,~] = min(d);
           idx1 = find(d/mind<10);
           if length(idx1)==length(idx2)
               break;
           end
       end
    end
    if length(idx1)==length(idx2)
        d1 = [];
        d2 = [];
        for k=1:length(idx1)
            d1 = [d1 distPnts(x1(:,idx1(k)),e)];
            d2 = [d2 distPnts(x2(:,idx2(k)),et)];
        end
        order1 = sort(d1);
        order2 = sort(d2);
        idx = order1==distPnts(x1(:,i),e);
        n = d2==order2(idx);
        midx = idx2(n); 
    end
    if left.id(i)==0 && right.id(midx)==0
        left.id(i) = countId;
        right.id(midx) = countId;
        countId = countId + 1;
    elseif left.id(i)~=0 && right.id(midx)==0
        right.id(midx) = left.id(i);  
    end
end
end

function d = SampsonDist(x1,x2,E)
EX1 = E*x1;
EtX2 = E'*x2;  
X2tEX1 = x2'*E*x1;
d =  X2tEX1^2/(EX1(1)^2 + EX1(2)^2 + EtX2(1)^2 + EtX2(2)^2);
end
