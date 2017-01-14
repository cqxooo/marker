function RealConcentric
% Real concentric circle experiment for 3D measurement

clc;
close;
clear;
imgFilename = {};
%--------------------------------------------------------
for i=1:3
    eval(['imgFilename{end+1} =''test/raw' int2str(i)  '.jpg'';']);
end
% for i=1:17
%     eval(['imgFilename{end+1} =''measureImages/low/raw' int2str(i)  '.jpg'';']);
% end
% for i=1:7
%     eval(['imgFilename{end+1} =''measureImages/high/raw' int2str(i)  '.jpg'';']);
% end
% for i=1:13
%     eval(['imgFilename{end+1} =''measureImages/light/raw' int2str(i)  '.jpg'';']);
% end
%ellipse detection
% for imLoop=1:length(imgFilename)
%     disp(imgFilename{imLoop})
%     [BigM(imLoop),SmallM(imLoop)] = imageProcess(imgFilename{imLoop});
%     BigM(imLoop) = deleteDuplicate(BigM(imLoop));
% %     [data(imLoop)] = img_process(imgFilename{imLoop});
% end
[pathstr,imNameC,~] = fileparts(imgFilename{1}); 
% save ([pathstr '/' imNameC '.mat'])
load ([pathstr '/' imNameC '.mat'])
%     
% %     =============================================
%     camera calibration
    K = [4529.355657679932800 0 2242.423458884737600;0 4530.395670821238100 1522.627084047865100 ;0 0 1];
%     K = calib_conic(imdata, K0);
for imLoop=1:length(BigM)
    fig = figure(imLoop);
    im = imread(imgFilename{imLoop});
    imshow(im);
    hold on; axis ij; 
    [conic BMdata(imLoop)] = getCenter(BigM(imLoop), K);
    for i=1:length(BMdata(imLoop).id)
        drawConic(conic(:,:,i), 'g-');    
        plot(BMdata(imLoop).point(1,i),BMdata(imLoop).point(2,i), 'r.');   
    end
    close(fig)
end
for imLoop=1:length(SmallM)
    fig = figure(imLoop);
    im = imread(imgFilename{imLoop});
    imshow(im);
    hold on; axis ij; 
    [conic SMdata(imLoop)] = getCenter(SmallM(imLoop), K);
    for i=1:length(SMdata(imLoop).c)
        drawConic(conic(:,:,i), 'g-');    
        plot(SMdata(imLoop).point(1,i),SMdata(imLoop).point(2,i),'r.');   
    end
    close(fig)
end
[E Proj] = getE_Proj(BMdata, K);
save test.mat
load test.mat
cloud = buildCloud(SMdata, E, Proj, K);
% scale = getScale(imgFilename, Proj,K);
% cloud.point = scale*cloud.point;
figure(1)
scatter3(cloud.point(1,:),cloud.point(2,:),cloud.point(3,:))
% for i=1:length(cloud)
%     n = length(cloud(i).id);
%     for j=1:n-1
%         for k=j+1:n
%             if j==1 && k==2
%                ratio = 19.2/distPnts(cloud(i).point(:,1),cloud(i).point(:,2));
%             end
%             d = distPnts(cloud(i).point(:,j),cloud(i).point(:,k))*ratio; %compute distance
%             disp(['distance between ' num2str(cloud(i).id(j)) ' to ' num2str(cloud(i).id(k)) ' is ' num2str(d) 'cm']);       %show computed distance
%         end        
%     end
% end
