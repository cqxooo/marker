function [lcor rcor] = findcorres(left, right)
lcor.id = [];
rcor.id = [];
lcor.point = [];
rcor.point = [];
for i = 1:length(left.id)
    if left.id(i)>0
        idx = find(left.id(i)==right.id);
        if ~isempty(idx)
            lcor.id = [lcor.id left.id(i)];
            rcor.id = [rcor.id right.id(idx)];
            lcor.point = [lcor.point left.point(:,i)]; 
            rcor.point = [rcor.point right.point(:,idx)];
        end
    end
end