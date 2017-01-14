function newM = fillMat(M, rect, imsize)
xmin = rect(1);
ymin = rect(2);
xmax = xmin+rect(3);
ymax = ymin+rect(4);
row = imsize(1);
col = imsize(2);
if (xmin<0 && ymin>0) || (xmin<0 && ymax<row)
    num = rect(3)-size(M,2)+1;
    newM = [zeros(rect(4)+1,num) M];
elseif (xmin>0 && ymin<0) || (xmax<col && ymin<0)
    num = rect(3)-size(M,1)+1;
    newM = [zeros(num,rect(3)+1);M];
elseif xmin<0 && ymin<0
    r = rect(3)-size(M,1)+1;
    c = rect(3)-size(M,2)+1;
    newM = [zeros(r,rect(3)+1);zeros(rect(4)+1-r,c) M];
elseif (xmin>0 && ymax>row) || (xmax<col && ymax>row)
    num = rect(3)-size(M,1)+1;
    newM = [M; zeros(num,rect(3)+1)];
elseif xmin<0 && ymax>col
    r = rect(3)-size(M,1)+1;
    c = rect(3)-size(M,2)+1;
    newM = [zeros(r,rect(3)+1);M zeros(rect(4)+1-r,c)];
elseif xmax>col && ymin>0 || (xmax>col && ymax<row)
    num = rect(3)-size(M,2)+1;
    newM = [M zeros(rect(3)+1,num)];
elseif xmax>row && ymin<0
    r = rect(3)-size(M,1)+1;
    c = rect(3)-size(M,2)+1;
    newM = [zeros(rect(4)+1-r,c) M;zeros(r,rect(3)+1)];
elseif xmax>row && ymax>col
    r = rect(3)-size(M,1)+1;
    c = rect(3)-size(M,2)+1;
    newM = [M zeros(rect(4)+1-r,c);zeros(r,rect(3)+1)];    
end