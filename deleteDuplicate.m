function data1 = deleteDuplicate(data)
data1.id = [];
data1.c = [];
% data1.corner = [];
for i=1:length(data.id)
    idx = find(data.id(i)==data.id);
    if length(idx)==1
       data1.id = [data1.id data.id(i)];
       data1.c = [data1.c data.c(i)];
%        data1.corner = [data1.corner data.corner(i)];
    end    
end