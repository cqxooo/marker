function cen = getTricen(cmat)
[c1, ~, ~] = getEigenVectors(cmat.C1, cmat.C2);
[c2, ~, ~] = getEigenVectors(cmat.C2, cmat.C3);
[c3, ~, ~] = getEigenVectors(cmat.C1, cmat.C3);
cen = mean([c1 c2 c3]')';