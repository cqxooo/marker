function test
BASE = ([0 1 0 0 0 0;0 0 0 1 0 0;0 0 0 0 0 1]);
figure(4);
plot3(BASE(1,:),BASE(3,:),-BASE(2,:),'b-','linewidth',2);
text(1,0,0,'X_c');
text(-1,1,0,'Z_c');
text(0,0,-1,'Y_c');
text(-1,-1,1,'O_c');