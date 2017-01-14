%
%
%   Date        : 2016-4-5
%   Objective   : Detect the position of particular markers
%   Input       : Images from mobile
%   Output      : Results( ID(ellipse_id) and conic(structure is conic.c1,
%                   conic.c2) )
%
%                 --Results->
%                      -- ellipse_id
%                      -- conic->
%                          --  c1
%                          --  c2
%
%   Author      : lampson
%
%

clc;
clear;
clc;

train_data = load('train.mat');
blank_angles = train_data.blank_angles;

% im = imread('006.png');
im = imread('test\test1_t.jpg');
figure; imshow(im)

im = rgb2gray(im);
[h,w] = size(im);

% Step 1:
% correct the ununiform illumination of images, dilate and erode
fprintf('\n Illumination correction... \n');
se = strel('disk',200);
se1 = strel('disk',50);
fse = imdilate(im,se);
background = imerode(fse,se1);
I2 = background - im;
im = imadjust(I2);


% Step 2:
% two kinds of binary images, cir_bw for circle detection, rec_bw for
% rectangle
rec_bw = im2bw(im);
cir_bw = ~rec_bw;

% imshow(cir_bw)

% Step 3: Circle detection
fprintf('\n Circles detection... \n');
[B,L] = bwboundaries(cir_bw, 'noholes');
STATS = regionprops(L, 'all');

cir_suspect_num = 0;

for i = 1 : length(STATS)  
        %properties of regions
        area = STATS(i).Area;
        maj = STATS(i).MajorAxisLength;
        min = STATS(i).MinorAxisLength;
        convex_area = STATS(i).ConvexArea;
        filled_area = STATS(i).FilledArea;
        
        % circle features
        diff1 = abs((area/(pi*(maj/2)*(min/2))-1));
        diff2 = maj/min;
        diff3 = maj/max(h,w);
        diff4 = convex_area/filled_area - 1;
        
        % search circles
        if((diff1<0.8) && (diff2 < 3) && (diff3 > 0.01) && (diff4 < 0.8))
            cir_suspect_num = cir_suspect_num + 1;
            if(cir_suspect_num == 1)
               cir_inds = {find(L==i)};
               [coord_x,coord_y] = find(L == i);
               cir_coords = {[coord_x,coord_y]};
            else
                cir_inds = [cir_inds; find(L==i)];
                [coord_x,coord_y] = find( L == i );
                cir_coords = [cir_coords;[coord_x,coord_y]];
            end
            %%%% show rectangle boundaries
%             boundary = B{i};
%             hold on;
%             plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
        end
end


% Step 4: Rectangle detection
fprintf('\n Rectangle detection... \n');
[B,L] = bwboundaries(rec_bw, 'noholes');
STATS = regionprops(L, 'all');

rec_suspect_num = 0;

for i = 1 : length(STATS)  
        %properties of regions
        maj = STATS(i).MajorAxisLength;
        convex_area = STATS(i).ConvexArea;
        filled_area = STATS(i).FilledArea;
        
        % rectangle features
        diff3 = maj/max(h,w);
        diff4 = convex_area/filled_area - 1;
        
        % search rectangles
        if( (diff4 < 0.1) && (diff3 > 0.01) )
            rec_suspect_num = rec_suspect_num + 1;
            
            %%%% show rectangle boundaries
%             boundary = B{i};
%             hold on;
%             plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
            
            %save convex image index
            convex_hull = STATS(i).ConvexHull;
            hull_bw = poly2mask(convex_hull(:,1),convex_hull(:,2),h,w);
            if(rec_suspect_num == 1)
               rec_inds = {find(hull_bw == 1)};
               [coord_x,coord_y] = find(hull_bw == 1);
               rec_coords = {[coord_x,coord_y]};
            else
               rec_inds = [rec_inds;find(hull_bw == 1)];
               [coord_x,coord_y] = find( hull_bw == 1 );
               rec_coords = [rec_coords;[coord_x,coord_y]];
            end
        end
end


% Step 5: Figure out the rectangles containing circles
fprintf('\n Region of interest selection... \n');
% count rectangle numbers in image
rect_num = 0;
for i = 1:rec_suspect_num
    
    rec_ite = rec_inds{i,1};
   
    for j = 1:cir_suspect_num
       
        cir_ite = cir_inds{j,1};
        % a circle is totally in the rectangle
        overlap_num = sum(ismember(cir_ite,rec_ite));
        if(overlap_num == length(cir_ite))

            
           % Find ROI, two kinds of images, one for rectangle the other for
           % circles
                   out = zeros(h,w);
                   out(rec_ite) = im(rec_ite);

                   ROI = uint8(out);
                   ROI = ~im2bw(ROI);

                   %rectangle
                   [B,L] = bwboundaries(im2bw(out), 'noholes');
                   STATS = regionprops(L, 'all');
                   %circle
                   [cir_B,cir_L] = bwboundaries(ROI, 'noholes');
                   cir_STATS = regionprops(cir_L, 'all');
           
                    
           
           %
           % Save all rentangle (x,y) and circle boundaries
           
           for bound_i = 1 : length(STATS)
               if (STATS(bound_i).Area / (h*w) > 0.0001)
                       %rectangle central coordinate
                       central_point = STATS(bound_i).Centroid;

%                         boundary = B{bound_i};
%                         hold on;
%                         plot(boundary(:,2), boundary(:,1), 'r-', 'LineWidth', 5);
%                         hold on;
%                         plot(central_point(2), central_point(1), 'r', 'LineWidth', 5);
                        %
                        % %%%  central circle
                        %
                        % show all shapes inside a rectangle
                        shape_num = 0;
                        for c_i = 1 : length(cir_STATS)
                            subcir_ind = find(cir_L == c_i);
                            subcir_sum = sum(ismember(subcir_ind,rec_ite));

                            maj = cir_STATS(c_i).MajorAxisLength;
                            diff3 = maj/max(h,w);
                            % region inside rectangle and big enough
                            if( (subcir_sum == length(subcir_ind)) && (diff3>0.005) )
                                temp_dist = norm( cir_STATS(c_i).Centroid - central_point);

                                shape_num = shape_num + 1;
                                if (shape_num == 1)
                                    marker_boundary = { cir_B{c_i} };
                                    centroid_dist = {temp_dist};
                                else
                                    marker_boundary = [marker_boundary;cir_B{c_i}];
                                    centroid_dist = [centroid_dist;temp_dist];
                                end
                            end
                        end
               end

           end
           
           rect_num = rect_num + 1;
           if (rect_num == 1)
               all_boundary = {marker_boundary};
               central_dist_set = {centroid_dist};
           else
               all_boundary = [all_boundary;{marker_boundary}];
               central_dist_set = [central_dist_set;{centroid_dist}];
           end

                
           break;
        end
   end
end








% Correct the markers to frontal
fprintf('\n Marker Correction ... \n');
for i = 1:rect_num
    dist = cell2mat( central_dist_set{i} );
    [~,sort_ind] = sort(dist);
    
    boundaries = all_boundary{i};
    subcir_boundary = boundaries{sort_ind(1)};
    
    trans_boundary = [];
    len_boundary = length(boundaries);
    if len_boundary < 2
        fprintf('Wrong rectangles');
        continue;
    end
    
    for j = 2:length(boundaries)
        trans_boundary = [trans_boundary;boundaries{sort_ind(j)}];
    end
    

        % ellipse fitting
        params = fitellipse(subcir_boundary(:,1),subcir_boundary(:,2));
        
        % rotation matrix
        rota = [cos(params(5)) -sin(params(5));sin(params(5)) cos(params(5))];
                            
        
        region_y = trans_boundary(:,1);
        region_x = trans_boundary(:,2);

        % rotate and affine transformation
        test_data = trans_boundary - repmat( [params(1) params(2)], length(trans_boundary), 1);
        test_data = test_data * rota;

        affi = [params(4)/params(3) 0; 0 1];
        test_data = test_data * affi;

        out_ind = test_data + repmat( [params(1) params(2)], length(trans_boundary), 1);
        out_ind = ceil(out_ind);

        % warp
        min_x = sort(out_ind(:,1));
        min_y = sort(out_ind(:,2));

           out_ind(:,1) = out_ind(:,1) - repmat( min_x(1), length(min_x),1 ) + 1; 
           out_ind(:,2) = out_ind(:,2) - repmat( min_y(1), length(min_y),1 ) + 1;
           
           cir_point = [params(1)-min_x(1),params(2)-min_y(2)];

        out = zeros(max(out_ind(:,2)),max(out_ind(:,1)));
        [out_h,out_w] = size(out);

        % original indexs and corrected index
        region_ind = (region_x-1)*h + region_y;
        out_index = out_ind(:,2) + (out_ind(:,1) - 1) * out_h;

        white_im = ones(size(im));
        out(out_index) = white_im(region_ind);
%         out = im2bw(out);
%         figure;
%         imshow( out );

        


        

        % 
        % Step 6: get blank anlges
        %
        x_0= cir_point(1);
        y_0 = cir_point(2);

        coord_x = out_ind(:,1) - x_0;
        coord_y = out_ind(:,2) - y_0;


        comple = zeros(size(coord_y));


        %%% x>0; y>0
        % nothing changes

        %%% x<0; y>0 && x<0,y<0
            clear left_ind;
            left_ind = find(coord_x<0);
        comple(left_ind) = pi;

        %%% x>0,y<0
        ind_1 = zeros(size(coord_y));
        ind_2 = zeros(size(coord_y));

            up_ind = find(coord_y<0);
            left_ind = find(coord_x>0);

            ind_1(up_ind) = 1;
            ind_2(left_ind) = 1;

        ind_roi = find( (ind_1 & ind_2) == 1);

        comple(ind_roi) = 2*pi;


        region_b_num = 5;
        region_num = 5;
        
        % % % dist = sqrt( coord_x.^2 + coord_y.^2 );
        % % 
        theta = atan( (coord_y) ./ (coord_x) ) + comple;
        distance = sqrt( (coord_x).^2 + (coord_y).^2 );
        
        
                    angle_num = 180;
                    bins = 0:2*pi/angle_num:2*pi;

                    angle_data = hist(theta,bins);


                    blank_num = zeros(region_num,1);
                    blank_block = 0;

                    %
                    % select angles
                    %
                    angle_i = 1;
                    while angle_i<=angle_num
                        if( angle_data(angle_i)+angle_data(angle_i+1) ~= 0 )                        
                            blank_block = blank_block + 1;
                            
                            if(blank_block>5)
                               fprintf('\n This rectangle is deleted \n');
                               break;
                            end
                            
                            while(angle_i<angle_num && ( angle_data(angle_i)+angle_data(angle_i+1) ~= 0) )
                               blank_num(blank_block) =  blank_num(blank_block) + 1;
                               angle_i = angle_i+1;
                            end
                            
                            if(angle_i == angle_num)
                               if( (angle_data(angle_i)~=0) && (angle_data(1)~=0) && (blank_block>1) )
                                  blank_num(1) = blank_num(1) + blank_num(blank_block);
                                  blank_num(blank_block) = 0;
                               end
                            end  
                        end
                        angle_i = angle_i + 1;
                    end
                            % When the image quality is too low
                            if(blank_block>5)
                               fprintf('\n Next rectangle \n');
                               continue;
                            end
                    
                    region_b_num = region_b_num - length(find(blank_num==0));
                    region_area = blank_num(1:region_b_num);                    
                    
                    [~,first_blank] = sort(region_area,'descend');
                    
                    
                    result = zeros(1,region_num);
                    
                    for s_i = 1:region_b_num
                        s_ind = first_blank(1) + s_i - 1;
                        if(s_ind>region_b_num)
                           s_ind =  mod(s_ind,region_b_num);
                        end
                        result(1,s_i) = region_area(s_ind);
                    end
        %
        % Output ID
        diff = blank_angles - repmat(result,length(blank_angles),1);
        
        dist = sum(diff.^2,2);
        [~,id] = sort(dist);
        
        ellipse_id = id(1);
        
        
                    %
                    %
                    % select points
                    %
                    sector = 1;
                    angle_i = 1;
                    select_points_num = ceil( 40/region_b_num );
                    
                    while angle_i<=angle_num
                        if( angle_data(angle_i)+angle_data(angle_i+1) ~= 0 )
                            
                            % if the first block is splited
                            if(angle_i == 1)
                                ite_i = angle_num;
                                while(ite_i>0 && ( angle_data(ite_i) + angle_data(ite_i-1) ~= 0) )
                                   ite_i = ite_i-1;
                                end
                                start_ind = ite_i;
                                
                                while(angle_i<angle_num && ( angle_data(angle_i)+angle_data(angle_i+1) ~=0 ))
                                   angle_i = angle_i+1;
                                end
                                end_ind = angle_i;
                                
                                %blank_block = blank_block - 1;
                                
                                gap = end_ind + ( angle_num - start_ind);
                                
                                select_angle_start = (start_ind + 1/4 * gap) * 2 / 180 * pi;
                                select_angle_end = (end_ind - 1/4 * gap) * 2 / 180 * pi;
                                
                                select_ind = find(theta>select_angle_start | theta<select_angle_end);
                                [~,order] = sort(distance(select_ind),'ascend');
                                final_ind = select_ind(order);
                                hold on; plot(region_x(final_ind(1:select_points_num)), region_y(final_ind(1:select_points_num)), 'r.'  );
                                
                                conic.c1 = [region_x(final_ind(1:select_points_num)) region_y(final_ind(1:select_points_num))];
                                
                                [~,order] = sort(distance(select_ind),'descend');
                                final_ind = select_ind(order);
                                hold on; plot(region_x(final_ind(1:select_points_num)), region_y(final_ind(1:select_points_num)), 'g.'  );
                                conic.c2 = [region_x(final_ind(1:select_points_num)) region_y(final_ind(1:select_points_num))];
                                
                            elseif(sector < blank_block)
                                    start_ind = angle_i;
                                    while(angle_i<angle_num && ( angle_data(angle_i)+angle_data(angle_i+1) ~=0 ) )
                                       angle_i = angle_i+1;
                                    end

                                    end_ind = angle_i;
                                    
                                    gap = end_ind - start_ind;
                                    
                                    select_angle_start = (start_ind + 1/4 * gap) * 2 / 180 * pi;
                                    select_angle_end = (end_ind - 1/4 * gap) * 2 / 180 * pi;

                                    select_ind = find(theta>select_angle_start & theta<select_angle_end);
                                    
                                    [~,order] = sort(distance(select_ind),'ascend');
                                    final_ind = select_ind(order);
                                    hold on; plot(region_x(final_ind(1:select_points_num)), region_y(final_ind(1:select_points_num)), 'r.'  );
                                    if(sector == 1)
                                        conic.c1 = [region_x(final_ind(1:select_points_num)) region_y(final_ind(1:select_points_num))];
                                    else
                                        conic.c1 = [conic.c1 ; region_x(final_ind(1:select_points_num)) region_y(final_ind(1:select_points_num))];
                                    end
                                    [~,order] = sort(distance(select_ind),'descend');
                                    final_ind = select_ind(order);
                                    hold on; plot(region_x(final_ind(1:select_points_num)), region_y(final_ind(1:select_points_num)), 'g.'  );
                                    if(sector == 1)
                                        conic.c2 = [region_x(final_ind(1:select_points_num)) region_y(final_ind(1:select_points_num))];
                                    else
                                        conic.c2 = [conic.c2 ; region_x(final_ind(1:select_points_num)) region_y(final_ind(1:select_points_num))];
                                    end
                            end
                            sector = sector + 1;
                            
                        end
                        
                        angle_i = angle_i + 1;
                    end
                    
                    if (i == 1)
                        Results.id = ellipse_id;
                        Results.conic = conic;
                    else
                        Results.id = [Results.id;ellipse_id];
                        Results.conic = [Results.conic;conic];
                    end
end

Results