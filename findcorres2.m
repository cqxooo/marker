function [lcor rcor] = findcorres2(left, right)
lcor.id = [];
rcor.id = [];
lcor.point = [];
rcor.point = [];
lcor.cmat = [];
rcor.cmat = [];
for i = 1:length(right.id)
    idx = find(left.id==right.id(i));
    if ~isempty(idx)
        lcor.id = [lcor.id left.id(idx)];
        rcor.id = [rcor.id right.id(i)];
        lcor.point = [lcor.point left.point(:,idx)]; 
        rcor.point = [rcor.point right.point(:,i)];
        lcor.cmat = [lcor.cmat left.cmat(idx)]; 
        rcor.cmat = [rcor.cmat right.cmat(i)];
    end
end