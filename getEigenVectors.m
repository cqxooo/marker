function [cen, v1, v2] = getEigenVectors(C1, C2)

%          C12 = inv(C1)*C2;
         C12 = inv(C2)*C1;

         %C12 = C12./C12(3,3);

            [V,D]=eig(C12,'nobalance');
            %[V,D]=eig(C12);
          V = [V(1,:)./V(3,:); V(2,:)./V(3,:); ones(1,3)];
        d = [D(1,1) D(2,2) D(3, 3)];
        dd = [abs(d(1)-d(2))+abs(d(1)-d(3)) abs(d(2)-d(1))+abs(d(2)-d(3)) abs(d(3)-d(1))+abs(d(3)-d(2))];
        [a, i] = max(dd);
        cen = V(:,i);
        
        k = 1:3;
        idx = find(k~=i);
        v1 = V(:,idx(1));
        v2 = V(:,idx(2));
        if v1(1)>v2(1)
            temp = v1;
            v1 = v2;  
            v2 = temp;
        end
%         hold on;
%         plot(cen(1), cen(2), 'r+');
%         plot(v1(1), v1(2), 'b+');
%         plot(v2(1), v2(2), 'b+');
%         hold off;

        
return; 