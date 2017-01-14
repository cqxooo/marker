function result = updateCloud(result, X)
for i = 1:length(X.id)
    idx = find(result.id==X.id(i));
    if isempty(idx)
        result.id = [result.id X.id(i)]; 
        result.point = [result.point X.point(:,i)];
%     else
%         result.point(:,idx) 
%         X.point(:,i)
%         result.point(:,idx) = (result.point(:,idx)+X.point(:,i))/2;
    end
end