function cps = getCp(cmat)
cps = [];
for i=1:length(cmat)
    A = [];
    [~, v1, v2] = getEigenVectors(cmat(i).C1, cmat(i).C2);
    A = [A v1 v2];
    [~, v1, v2] = getEigenVectors(cmat(i).C2, cmat(i).C3);
    A = [A v1 v2];
    [~, v1, v2] = getEigenVectors(cmat(i).C1, cmat(i).C3);
    A = [A v1 v2];
    [~, ~,v] = svd(A');
    l = v(:,end);
    lh = -sign(l(3))*l/norm(l(1:2));
    [cp1, cp2] = getIntersectLineConic(lh, cmat(i).C1);
    cps = [cps cp1 cp2];
    [cp1, cp2] = getIntersectLineConic(lh, cmat(i).C2);
    cps = [cps cp1 cp2];
    [cp1, cp2] = getIntersectLineConic(lh, cmat(i).C3);
    cps = [cps cp1 cp2];
end

