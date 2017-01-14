%
%
%   Date        : 2016-6-3
%   Objective   : process the input image, get id & conics
%   Author      : lampson
%
%   Output      : BigMarker.id ; BigMarker.c ; BigMarker.corner;
%                 SmallMarker.c
%


function [BigMarker,SmallMarker] = imageProcess(file_name)

im = imread(file_name);
% figure;
% imshow(im)

im = rgb2gray(im);
[h,w] = size(im);
%==========revise==============
BigMarker.id = [];
BigMarker.c = [];
SmallMarker.c = [];
%==============================
% Step 1:
% correct the ununiform illumination of images, dilate and erode
fprintf('\n Illumination correction... \n');
se = strel('disk',100);
se1 = strel('disk',100);
fse = imdilate(im,se);
background = imerode(fse,se1);
I2 = background - im;
im = imadjust(I2);


train_data = load('train.mat');
blank_angles = train_data.train_data.blank_angles;
id_index = train_data.train_data.id_index;


% Step 2:
% two kinds of binary images, cir_bw for circle detection, rec_bw for
% rectangle
rec_bw = im2bw(im);
cir_bw = ~rec_bw;

rec_im = zeros(h,w);

% Step 4: Rectangle detection
fprintf('\n Rectangle detection... \n');
[B,L] = bwboundaries(rec_bw, 'noholes');
STATS = regionprops(L, 'all');

mark_num = 0;
big_num = 0;
small_num = 0;

for i = 1 : length(STATS)  
        %properties of regions
        min_axis = STATS(i).MinorAxisLength;
        convex_area = STATS(i).ConvexArea;
        filled_area = STATS(i).FilledArea;
        
        % rectangle features
        diff4 = min_axis/max(h,w);
        diff5 = convex_area/filled_area - 1;
        
        % search rectangles
        if( (diff5 < 0.3) && (diff4 > 0.01) )
            
            %%%% show rectangle boundaries
%             boundary = B{i};
%             hold on;
%             plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)

            rec_center = STATS(i).Centroid;
            
            
            %save convex image index
            convex_hull = STATS(i).ConvexHull;
            hull_bw = poly2mask(convex_hull(:,1),convex_hull(:,2),h,w);
            
            rec_inds = {find(hull_bw == 1)};
               
            
%             figure;
%             imshow(im(rec_inds));
            ind = rec_inds{1,1};
            rec_im(ind) = rec_bw(ind);
            
            cir_im = zeros(h,w);
            cir_im(ind) = cir_bw(ind);
            
%             figure;
%             imshow(cir_im)
            
            [inner_B,inner_L] = bwboundaries(cir_im, 'noholes');
            REGIONS = regionprops(inner_L,'all');
            
            region_length = length(REGIONS);
            areas = zeros(region_length,1);
            all_area = zeros(region_length,1);
            centers = zeros(region_length,2);
            
            for j = 1:region_length
                areas(j) = REGIONS(j).Area;
                all_area(j) = REGIONS(j).FilledArea;
                centers(j,:) = REGIONS(j).Centroid;
            end
            
            dist_centers = centers - repmat(rec_center,[region_length,1]);
            dist_centers = sum(dist_centers.^2,2);
            
            [~,order] = sort(dist_centers);
            
                j = order(1);
                %properties of regions
                inner_area = REGIONS(j).Area;
                inner_maj = REGIONS(j).MajorAxisLength;
                inner_min_a = REGIONS(j).MinorAxisLength;

                % circle features
                diff1 = abs((inner_area/(pi*(inner_maj/2)*(inner_min_a/2))-1));
                diff2 = inner_maj/inner_min_a;
                diff3 = inner_maj/max(h,w);
                
                if( (diff1<0.2) && (diff2 < 3) && (diff3 > 0.005) ) 
                    mark_num = mark_num + 1;
                    fprintf('\n %d marker \n',mark_num);
                    
                    sector_boundary = [];
                    % order pixel numbers in sectors, find the difference
                    % between big and small markers
                    [~,sector_order] = sort(areas,'descend');
                    for ite = 1:region_length
                       c_order = sector_order(ite);
                       if( c_order ~= j )
                          if( areas(c_order) > inner_area*0.3 && ...
                                  areas(c_order)/all_area(c_order)>0.8 &&...
                                  areas(c_order)/all_area(c_order)<1.2)
                              sector_boundary = [sector_boundary;inner_B{c_order}];
                          end
                       else
                           circle_boundary = inner_B{c_order};
                       end
                    end
                    
                    if(~isempty(sector_boundary))
                        big_num = big_num + 1;
                        
                        
%                         hold on;
%                         plot(corner(:,1),corner(:,2),'r*');
% 
%                         hold on;
%                         plot(circle_boundary(:,2), circle_boundary(:,1), 'g', 'LineWidth', 2)

%                         hold on;
%                         plot(sector_boundary(:,2), sector_boundary(:,1), 'g', 'LineWidth', 2)


                        id = getID(sector_boundary,circle_boundary,blank_angles,id_index);
%                         corner = k_means(convex_hull);
                        
                        
                        if(big_num > 0 && big_num == 1)
                            BigMarker.id = id;
                            BigMarker.c = {circle_boundary};
%                             BigMarker.corner = {corner};
                        else
                            BigMarker.id = [BigMarker.id ; id];
                            BigMarker.c = [BigMarker.c ; {circle_boundary}];
%                             BigMarker.corner = [BigMarker.corner ; {corner}];
                        end

                    else
                        small_num = small_num + 1;
                        circle_boundary = inner_B{j};
%                         hold on;
%                         plot(circle_boundary(:,2), circle_boundary(:,1), 'r', 'LineWidth', 2)
                        
                        if(small_num == 1)
                           SmallMarker.c = {circle_boundary};
                        else
                           SmallMarker.c = [SmallMarker.c ; {circle_boundary}];
                        end
                    end
                end
            
    
        end
end