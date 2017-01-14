function RR  = rectify(C1, C2, K, ls, lh)
%===========================================
%rectify
%===========================================
x0 = K(1:2,3);
xfl = K(1,1);
yfl = K(2,2);
fig = 1;

% x_ = getPerPtLine(x0, ls);
% v1 = inv(K)*[x0;1]; v1=v1/norm(v1);
% v2 = inv(K)*[x_;1]; v2=v2/norm(v2);
% ct = v1'*v2;
% theta = acos(ct);
% v = cross(v2,v1);
% v = v/norm(v);
% R0 = rot(theta,v,0);
x_ = (-ls(2)*x0(2)-ls(3))/ls(1);
d = x0(1) - x_;
theta = atan2(d,xfl);
R0 = rot(theta,'y',0);
H0 = K*R0*inv(K);

ls0 = inv(H0')*ls;
ls0 = -sign(ls0(3))*ls0/norm(ls0(1:2));

if fig
    fig1=figure(99);
    H = H0;    
    drawConic(H'*C1*H, 'b-');
    hold on; axis ij; axis([0 x0(1)*2 0 x0(2)*2]);
    drawConic(H'*C2*H, 'b-');
    plot(x0(1), x0(2),'r*','markersize',12);
    plot(x_, x0(2),'b*','markersize',12);
    drawLine(ls0,'b-');
    drawLine(inv(H)'*lh,'g-');    
end    
   
%     axis([0 640 0 480]);
%     set(gca, 'xtick', [], 'ytick', []);


%===========================================
% rotate the camera so that the axis is vertical
R1 = [ls0(1) ls0(2) 0
     -ls0(2) ls0(1) 0
        0      0    1];

H1 = K*R1*R0*inv(K);
%===========================================
% rectified contur in image

if fig
    fig2=figure(98);
    H = H1;
    drawConic(H'*C1*H, 'b-');
    hold on; axis ij; axis([0 x0(1)*2 0 x0(2)*2]);
    drawConic(H'*C2*H, 'b-');
    plot(x0(1), x0(2),'r*','markersize',12);
    drawLine(inv(H)'*ls,'b-');
    drawLine(inv(H)'*lh,'g-');       
end

%     axis([0 640 0 480]);
%     set(gca, 'xtick', [], 'ytick', []);



%===========================================
% rotate the camera so that the intersection of ls and lh align with
% principle point.

% x_ = getPerPtLine(x0, inv(H1)'*lh);
% v1 = inv(K)*[x0;1]; v1=v1/norm(v1);
% v2 = inv(K)*[x_;1]; v2=v2/norm(v2);
% ct = v1'*v2;
% theta = acos(ct);
% v = cross(v2,v1);
% v = v/norm(v);
% R2 = rot(theta,v,0);
lh0 = inv(H1')*lh;
lh0 = -sign(lh0(3))*lh0/norm(lh0(1:2));
y_ = (-lh0(1)*x0(1)-lh0(3))/lh0(2);
d = y_ - x0(2);
theta = atan2(d,yfl);
R2 = rot(theta,'x',0);
H2 = K*R2*R1*R0*inv(K);

ls0 = inv(H2')*ls;
ls0 = -sign(ls0(3))*ls0/norm(ls0(1:2))

lh0 = inv(H2')*lh;
lh0 = -sign(lh0(3))*lh0/norm(lh0(1:2))

if fig
    fig3=figure(97);
    H = H2;
    drawConic(H'*C1*H, 'b-');
    hold on; axis ij; axis([0 x0(1)*2 0 x0(2)*2]);
    drawConic(H'*C2*H, 'b-');
    plot(x0(1), x0(2),'r*','markersize',12);
    plot(x0(1), y_,'b*','markersize',12);
    drawLine(inv(H)'*ls,'b-');
    drawLine(inv(H)'*lh,'g-'); 
end

RR = inv(R2*R1*R0);
if fig
    close(fig1);
    close(fig2);
    close(fig3);
end
% v = H2*vp;
% v = [v(1,:)./v(3,:); v(2,:)./v(3,:); v(3,:)./v(3,:)];
% plot(v(1,1),v(2,1),'r*','MarkerSize',12);
% plot(v(1,2),v(2,2),'b*','MarkerSize',12);
% cp = H2*cp;
% cp = [cp(1,:)./cp(3,:); cp(2,:)./cp(3,:); cp(3,:)./cp(3,:)];
% [alpha beta] = getAngle(v,inv(H2)'*lh, inv(H2)'*ls, cp);

return;