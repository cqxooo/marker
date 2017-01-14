function RealConcentric1
% Real concentric circle experiment for 3D measurement

clc;
close;
clear;
imgFilename = {};
%--------------------------------------------------------
for i=1:5
    eval(['imgFilename{end+1} =''test\test' int2str(i)  '.jpg'';']);
end
%ellipse detection

% for imLoop=1:length(imgFilename)
%     disp(imgFilename{imLoop})
%     imdata(imLoop) = img_process(imgFilename{imLoop});
% end
% save data.mat
load data.mat
    
    %=============================================
    %camera calibration
    K = [3657.56395 0 2211.45963;0 3690.11165 1852.34732;0 0 1];

for imLoop=1:length(imdata)
    fig = figure(imLoop);
    im = imread(imgFilename{imLoop});
    imshow(im);
    hold on; axis ij; 
    [conic data(imLoop)] = getCenter(imdata(imLoop), K);
    for i=1:length(conic)
        drawConic(conic(i).C1, 'r-');    
        drawConic(conic(i).C2, 'r-'); 
        drawConic(conic(i).C3, 'r-');   
        plot(data(imLoop).cens(1,i),data(imLoop).cens(2,i),'k.');   
        plot(data(imLoop).conic(i).c1(:,1),data(imLoop).conic(i).c1(:,2),'y.');
        plot(data(imLoop).conic(i).c2(:,1),data(imLoop).conic(i).c2(:,2),'g.');
        plot(data(imLoop).conic(i).c3(:,1),data(imLoop).conic(i).c3(:,2),'b.');
    end
    close(fig)
end
save test.mat
load test.mat
cloud = refine3D(cor, K);
for i=1:length(cloud)
    n = length(cloud(i).id);
    for j=1:n-1
        for k=j+1:n
            if j==1 && k==2
               ratio = 16.7/distPnts(cloud(i).point(:,1),cloud(i).point(:,2));
            end
            d = distPnts(cloud(i).point(:,j),cloud(i).point(:,k))*ratio; %compute distance
            disp(['distance between ' num2str(cloud(i).id(j)) ' to ' num2str(cloud(i).id(k)) ' is ' num2str(d) 'cm']);       %show computed distance
        end        
    end
end
