function [cps data] = getCpCenter(left, right, imLoop)
point = [];
idx = [];
for i = 1:length(left.id)
    left.cmat(i) = conicRobustFit(left.conic(i),3);
%     conic(i) = conicFit_real(data(i),3);
    if isempty(left.cmat(i).C1) || isempty(left.cmat(i).C2) || isempty(left.cmat(i).C3) 
        idx = [idx i];
    else
        cen = getTricen(left.cmat(i));
        point = [point cen];
    end
end
if ~isempty(idx)
    left.id(idx) = [];
    left.cmat(idx) = [];
end
left.point = point;
point = [];
idx = [];
for i = 1:length(right.id)
    right.cmat(i) = conicRobustFit(right.conic(i),3);
%     conic(i) = conicFit_real(data(i),3);
    if isempty(right.cmat(i).C1) || isempty(right.cmat(i).C2) || isempty(right.cmat(i).C3) 
        idx = [idx i];
    else
        cen = getTricen(right.cmat(i));
        point = [point cen];
    end
end
if ~isempty(idx)
    right.id(idx) = [];
    right.cmat(idx) = [];
end
right.point = point;
[left right] = findcorres2(left, right);
[F, bestidx] = ransacF(left, right);
cps = [];
if imLoop ==1
    cp = getCp(left.cmat(bestidx));
    cps = [cps cp];
end
    cp = getCp(right.cmat(bestidx));
    cps = [cps cp];


