function [pl pr] = getIntersections(C1, C2, l1, l2)
tl1 = [];
tr1 = [];
[p1, p2] = getIntersectLineConic(l1, C1);
tl1 = [tl1 p1];
tr1 = [tr1 p2];

[p1, p2] = getIntersectLineConic(l2, C1);
tl1 = [tl1 p1];
tr1 = [tr1 p2];
tl1 = sortrows(tl1',2)';
tr1 = sortrows(tr1',2)';

tl2 = [];
tr2 = [];
[p1, p2] = getIntersectLineConic(l1, C2);
tl2 = [tl2 p1];
tr2 = [tr2 p2];

[p1, p2] = getIntersectLineConic(l2, C2);
tl2 = [tl2 p1];
tr2 = [tr2 p2];
tl2 = sortrows(tl2',2)';
tr2 = sortrows(tr2',2)';
pl = [tl1 tl2];
pr = [tr1 tr2];
% fig = figure(2);
% axis ij;
% axis([0,640,0,480])
% hold on;
% drawConic(C1, 'b-');
% drawConic(C2, 'b-');
% drawLine(l1, 'r-');
% drawLine(l2, 'r-');
% for i=1:size(pl,2)
%     plot(pl(1,i),pl(2,i),'k*');
%     plot(pr(1,i),pr(2,i),'k*');
% end
% close(fig) 
return;