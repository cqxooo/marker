function [left right] = findcorres3(left, right, limf, rimf, winsize)
leftI = imread(limf);
rightI = imread(rimf);
[row, col, ~] = size(rightI);
imsize = [row, col];
countId = max(left.id)+1;
for i=1:length(left.point)
    i
    tmp = [];
    rect1 = round([left.point(1:2,i)'-winsize, 2*winsize, 2*winsize]);
    template1 = im2bw(imcrop(leftI, rect1));
    if std(double(template1(:)))==0
        continue;
    end
    if size(template1,1)<2*winsize+1 || size(template1,2)<2*winsize+1
        template1 = fillMat(template1, rect1, imsize);
    end
%     figure(1)
%     imshow(template1)
    for j=1:length(right.point)
        rect2 = round([right.point(1:2,j)'-winsize, 2*winsize, 2*winsize]);
        template2 = im2bw(imcrop(rightI, rect2));
        if size(template2,1)<2*winsize+1 || size(template2,2)<2*winsize+1
            template2 = fillMat(template2, rect2, imsize);
        end
%         figure(2)
%         imshow(template2)
        c = normxcorr2(template1, template2);
        tmp = [tmp sum(c(:))];
    end
    [~,idx] = max(tmp);
    if left.id(i)==0 && right.id(idx)==0
        left.id(i) = countId;
        right.id(idx) = countId;
        countId = countId+1;
    elseif left.id(i)~=0 && right.id(idx)==0
        right.id(idx) = left.id(i);
    end
end