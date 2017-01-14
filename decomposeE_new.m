function [R, T] = decomposeE_new(E, x1, x2)

% [R, t] = decomposeE(E);
%
% This function decomposes an essential matrix 
% in a rotation and a translation.
%
% E - essential matrix;
% R - rotation;
% t- translation.

[U D V] = svd(E);

W = [0 -1 0
     1  0 0
     0  0 1];

t = U(:,3)./norm(U(:,3));

% Two possibilities:
rot1 = U * W  * V';
rot2 = U * W' * V';

if det(rot1)<0
    rot1 = -rot1;
end
if det(rot2)<0
    rot2 = -rot2;
end    

% 4 possible choices of the camera matrix P2
% choices of R and 2 possible signs of t.
p1 = [eye(3) zeros(3,1)];
P2(:,:,1) = [rot1 t]; 
P2(:,:,2) = [rot1 -t]; 
P2(:,:,3) = [rot2 t]; 
P2(:,:,4) = [rot2 -t]; 

for i=1:4
    p2 = P2(:,:,i);
    A = [x1(1)*p1(3,:)-x1(3)*p1(1,:);
         x1(2)*p1(3,:)-x1(3)*p1(2,:);
         x2(1)*p2(3,:)-x2(3)*p2(1,:);
         x2(2)*p2(3,:)-x2(3)*p2(2,:)];
    [u d v] = svd(A);
    d1 = getDepth(p1, v(:,end));
    d2 = getDepth(p2, v(:,end));
    if d1 > 0 && d2 > 0
        R = p2(:,1:3);
        T = p2(:,end);   
    end
end

