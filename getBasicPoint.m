function getBasicPoint
clc;
close;
clear;
imgFilename = {};
% for i=1:4
%     eval(['imgFilename{end+1} =''test/test' int2str(i)  '.jpg'';']);
% end
% for i=1:17
%     eval(['imgFilename{end+1} =''measureImages/low/raw' int2str(i)  '.jpg'';']);
% end
% for i=1:7
%     eval(['imgFilename{end+1} =''measureImages/high/raw' int2str(i)  '.jpg'';']);
% end
for i=1:12
    eval(['imgFilename{end+1} =''measureImages/light/raw' int2str(i)  '.jpg'';']);
end
for imLoop=1:length(imgFilename)
    disp(imgFilename{imLoop})
    fig = figure(imLoop)
    imshow(imread(imgFilename{imLoop}))
    hold on 
    for i=1:2
        [point(1,i), point(2,i)] = ginput(1);
        plot(point(1,i), point(2,i), 'r.');
    end
    close(fig)
    basic(imLoop).point = [point;ones(1,2)];
end
[pathstr,imNameC,~] = fileparts(imgFilename{1}); 
save ([pathstr '/basic.mat'])

