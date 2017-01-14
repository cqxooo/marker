function [conic data] = getCenter2(data, K)
IAC = inv(K*K');
point = [];
idx = [];
for i = 1:length(data.id)
    conic(i) = conicRobustFit(data.conic(i),3);
%     conic(i) = conicFit_real(data(i),3);
    if isempty(conic(i).C1) || isempty(conic(i).C2) || isempty(conic(i).C3) 
        idx = [idx i];
    else
        cen = getTricen(conic(i));
        point = [point cen];
    end
end
if ~isempty(idx)
    data.id(idx) = [];
    data.conic(idx) = [];
    conic(idx) = [];
end
data.point = point;