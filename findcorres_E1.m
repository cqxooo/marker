function [left right] = findcorres_E1(left, right, E, K)
t = 10^(-8);
x1 = inv(K)*left.point;
x1 = x1./repmat(x1(3,:),3,1);
x2 = inv(K)*right.point;
x2 = x2./repmat(x2(3,:),3,1);
countId = max(left.id)+1;
[et e] = getEpipole(E);
et = et/et(3);
e = e/e(3);
for i=1:size(x1,2)
    tmp = [];
    for j=1:size(x2,2)
            EX1 = E*x1(:,i);
            EtX2 = E'*x2(:,j);  
            X2tEX1 = x2(:,j)'*E*x1(:,i);
            d =  X2tEX1^2/(EX1(1)^2 + EX1(2)^2 + EtX2(1)^2 + EtX2(2)^2);
            tmp = [tmp;d];
    end    
    idx2 = find(tmp<t);
    if ~isempty(idx2)
        tmp = [];
        for j=1:size(x1,2)
            EX1 = E*x1(:,j);
            EtX2 = E'*x2(:,idx2(1));  
            X2tEX1 = x2(:,idx2(1))'*E*x1(:,j);
            d =  X2tEX1^2/(EX1(1)^2 + EX1(2)^2 + EtX2(1)^2 + EtX2(2)^2);
            tmp = [tmp;d];
        end
        idx1 = find(tmp<t);
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
           if left.id(i)==0 && right.id(idx2(n))==0
               left.id(i) = countId;
               right.id(idx2(n)) = countId;
               countId = countId+1;
           elseif left.id(i)~=0 && right.id(idx2(n))==0
               if isempty(find(right.id==left.id(i), 1))
                   right.id(idx2(n)) = left.id(i);
               end
           elseif left.id(i)==0 && right.id(idx2(n))~=0
               if isempty(find(left.id==right.id(idx2(n)), 1))
                   left.id(i) = right.id(idx2(n));
               end
           end
        end
    end
end

    
%     if length(idx)==1
%         if left.id(i)==0 && right.id(idx)==0
%             left.id(i) = countId;
%             right.id(idx) = countId;
%             countId = countId+1;
%         elseif left.id(i)~=0 && right.id(idx)==0
%             if isempty(find(right.id==left.id(i), 1))
%                 right.id(idx) = left.id(i);
%             end
%         elseif left.id(i)==0 && right.id(idx)~=0
%             if isempty(find(left.id==right.id(idx), 1))
%                 left.id(i) = right.id(idx);
%             end
%         end
%     end
% end