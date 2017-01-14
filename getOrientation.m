function concen = getOrientation(conic, K, lhs, cens)
for i=1:size(lhs,2)
    lh = lhs(:,i);
    cen = cens(:,i);
    vy = K*K'*lh;
    ls = cross(vy, cen);
    ls = -sign(ls(3))*ls/norm(ls(1:2));
    [concen(i).RR, concen(i).t]  = rectify1(conic(i).C1, conic(i).C2, K, ls, lh);
    concen(i).cen = cen;    
end