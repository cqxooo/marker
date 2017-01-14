function X = triangulate2(x1, x2, P1, P2)
AA = [];
for i=1:size(x1,2)
    A = [x1(1,i)*P1(3,:,i)-x1(3,i)*P1(1,:,i);
        x1(2,i)*P1(3,:,i)-x1(3,i)*P1(2,:,i);
        x2(1,i)*P2(3,:,i)-x2(3,i)*P2(1,:,i);
        x2(2,i)*P2(3,:,i)-x2(3,i)*P2(2,:,i)];
    A1n = sqrt(sum(A(1,:).*A(1,:)));
    A2n = sqrt(sum(A(2,:).*A(2,:)));
    A3n = sqrt(sum(A(3,:).*A(3,:)));
    A4n = sqrt(sum(A(4,:).*A(4,:))); 
    A = [A(1,:)/A1n;
        A(2,:)/A2n;
        A(3,:)/A3n;
        A(4,:)/A4n];
    AA = [AA;A];
end
    [u d v] = svd(AA);
    w = v(:,end); 
    X = w(1:3)/w(4);

    
