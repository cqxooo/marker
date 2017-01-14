function [conic data] = getCenter(data, K)
IAC = inv(K*K');
point = [];
idx = [];
for i = 1:length(data.c)
    conic(:,:,i) = conicFit(data.c(i));
%     conic(i) = conicFit_real(data(i),3);
    if isempty(conic(:,:,i)) 
        idx = [idx i];
    else
        [cen, ~, ~] = getEigenVectors(IAC, conic(:,:,i));
        point = [point cen];
    end
end
if ~isempty(idx)
    if isfield(data,'id')
        data.id(idx) = [];
    end    
    data.c(idx) = [];
    conic(:,:,idx) = [];
end
data.point = point;
if ~isfield(data,'id')
    data.id = zeros(1,length(data.c));
end