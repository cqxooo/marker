function [R, t]  = rectify1(C1, C2, K, ls, lh)
x0 = K(1:2,3);
vy = K*K'*lh;
vx = K*K'*ls;
vz = cross(ls, lh);
c = inv(K)*vz;
t = c*10/norm(c);
r1 = inv(K)*vx;
r2 = inv(K)*vy;
r3 = cross(r1,r2);
R = [r1/norm(r1) r2/norm(r2) r3/norm(r3)];

% c2 = K'*C2*K;
% [V,D] = eig(c2);
% d = diag(D);
% [~,idx]=sort(d);
% d = d(idx(3:-1:1));
% R = V(:,[3 2 1]);
% n = sqrt((d(1)-d(2))/(d(1)-d(3)));
% m = 0;
% l = sqrt((d(2)-d(3))/(d(1)-d(3)));
% lh1 = inv(K*R)'*[l m n]';
% lh1 = -sign(lh1(3))*lh1/norm(lh1(1:2));

H = K*inv(R)*inv(K); 
ls0 = inv(H)'*ls;
ls0 = -sign(ls0(3))*ls0/norm(ls0(1:2))
lh0 = inv(H)'*lh;
lh0 = -sign(lh0(3))*lh0/norm(lh0(1:2))

fig = figure(99);
drawConic(H'*C1*H, 'b-');
hold on; axis ij; axis([0 x0(1)*2 0 x0(2)*2]);
drawConic(H'*C2*H, 'b-');
plot(x0(1), x0(2),'r*','markersize',12);
drawLine(ls0,'b-');
drawLine(lh0,'g-');   
close(fig)