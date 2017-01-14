function RealConcentric_autocalib
% Real concentric circle experiment for 3D measurement

clc;
close;
clear;
imgFilename = {};
%--------------------------------------------------------
for i=1:4
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
%     camera calibration
      K0 = [4495.762756580495989 0 2663.141192930780107;0 4573.667054720602209 1850.660882632053244;0 0 1];
      K = calib_conic(imdata, K0);
%=============================================
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
        plot(data(imLoop).point(1,i),data(imLoop).point(2,i),'k.');   
        plot(data(imLoop).conic(i).c1(:,1),data(imLoop).conic(i).c1(:,2),'y.');
        plot(data(imLoop).conic(i).c2(:,1),data(imLoop).conic(i).c2(:,2),'g.');
        plot(data(imLoop).conic(i).c3(:,1),data(imLoop).conic(i).c3(:,2),'b.');
    end
    close(fig)
end
save test.mat
load test.mat
cloud = refine3D(data, K);
for i=1:length(cloud)
    n = length(cloud(i).id);
    for j=1:n-1
        for k=j+1:n
            if j==1 && k==2
               ratio = 14.4/distPnts(cloud(i).point(:,1),cloud(i).point(:,2));
            end
            d = distPnts(cloud(i).point(:,j),cloud(i).point(:,k))*ratio; %compute distance
            disp(['distance between ' num2str(cloud(i).id(j)) ' to ' num2str(cloud(i).id(k)) ' is ' num2str(d) 'cm']);       %show computed distance
        end        
    end
end
