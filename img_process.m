function Results = img_process(imf)
%
%   Date        : 2016-5-1
%   Objective   : Detect the position of particular markers
%   Input       : Images from mobile
%   Output      : Results( ID(ellipse_id) and conic(structure is conic.c1,
%                   conic.c2,conic.c3) )
%
%                 --Results->
%                      -- ellipse_id
%                      -- conic->
%                          --  c1
%                          --  c2
%                          --  c3
%
%   Author      : lampson
%
%

train_data = load('train.mat');
blank_angles = train_data.blank_angles;

% im = imread('006.png');
im = imread(imf);
% figure;
% imshow(im)
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
        min_a = STATS(i).MinorAxisLength;
        convex_area = STATS(i).ConvexArea;
        filled_area = STATS(i).FilledArea;
        
        % circle features
        diff1 = abs((area/(pi*(maj/2)*(min_a/2))-1));
        diff2 = maj/min_a;
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

                        boundary = B{bound_i};
%                         hold on; plot(boundary(:,2), boundary(:,1), 'r-', 'LineWidth', 5);

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
           fprintf('\n      --> Detected %d markers... \n',rect_num);
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
start_result = 0;
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
    
    for j = 2:len_boundary
        trans_boundary = [trans_boundary;boundaries{sort_ind(j)}];
    end
    
    %check whether the image is separated correctly
    sectorTwice = 2*(length(boundaries)-1);
                        
    % add the inner circle to Result
    conic.c1 = [subcir_boundary(:,2) subcir_boundary(:,1)];

%     hold on; plot(subcir_boundary(:,2),subcir_boundary(:,1),'bo')
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


        region_num = 20;
        
        % % % dist = sqrt( coord_x.^2 + coord_y.^2 );
        % % 
        theta = atan( (coord_y) ./ (coord_x) ) + comple;
        distance = sqrt( (coord_x).^2 + (coord_y).^2 );
        
        
                    angle_num = 90;
                    bins = 0:2*pi/angle_num:2*pi;
                    bins_num = angle_num + 1;
                    
                    angle_data = hist(theta,bins);


                    section_num = zeros(region_num,2);
                    blank_block = 0;

                    %
                    % select angles
                    %
                    angle_i = 1;
                    while angle_i<=bins_num
                        
                        blank_block = blank_block + 1;
                        
                        if(angle_i <= bins_num && angle_data(angle_i) ~= 0)
                            %label 1-> not blank
                            section_num(blank_block,2) = 1;
                            while(angle_i <= bins_num && angle_data(angle_i) ~= 0)
                                section_num(blank_block,1) =  section_num(blank_block,1) + 1;
                                
                                if(angle_i == bins_num && blank_block > sectorTwice)
                                    section_num(1,1) = section_num(1,1) + section_num(blank_block,1);
                                    section_num(blank_block,1) = 0;
                                    section_num(blank_block,2) = 0;
                                end
                                
                                angle_i = angle_i + 1;
                            end
                        else
                            %label 2-> blank
                            section_num(blank_block,2) = 2;
                            
                            while(angle_i <= bins_num && angle_data(angle_i) == 0)
                                section_num(blank_block,1) =  section_num(blank_block,1) + 1;
                                
                                if(angle_i == bins_num && blank_block > sectorTwice)
                                    section_num(1,1) = section_num(1,1) + section_num(blank_block,1);
                                    section_num(blank_block,1) = 0;
                                    section_num(blank_block,2) = 0;
                                end
                                
                                angle_i = angle_i + 1;
                            end
                        end
                    end

                    
                     region_b_num = region_num - length(find(section_num(:,1)==0));
                    region_area = section_num(1:region_b_num,:);
                    
                    
                    if(region_b_num ~= sectorTwice)
                        fprintf('Image process fail...');
                        if(i==1)
                           start_result = start_result + 1; 
                        end
                        continue;
                    end                    
                    
                    sector_ind = find(section_num(:,2)==1);
                    sector_area = section_num( sector_ind,1 );
                    [~,first_blank] = sort(sector_area,'descend');
                    
                    
                    result = zeros(1,region_num);
                    
                    s_ind = sector_ind( first_blank(1) );
                    for s_i = 1:region_b_num
                        
                        if(s_ind>region_b_num)
                           s_ind =  mod(s_ind,region_b_num);
                        end
                        result(1,s_i) = region_area(s_ind);
                        s_ind = s_ind + 1;
                    end
        %
        % Output ID
        diff = blank_angles - repmat(result,length(blank_angles),1);
        
        dist = sum(diff.^2,2);
        [~,id] = sort(dist);
        
        ellipse_id = id(1);
        
        
                    %
                    %
                    % Step 7: select points
                    %
                    sector = 1;
                    angle_i = 1;
                    sector_num = sectorTwice/2;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%
                    %
                    %
                    %   alert the size of selected points
                    %
                    %%%%%%%%%%%%%%%%%%%%%%%%
                    select_middle_points_num = ceil( 20/(region_b_num/2) );
                    select_out_points_num = ceil( 50/(region_b_num/2) );
                    
                    while (angle_i<=bins_num)
                        if( angle_data(angle_i) ~= 0 )
                            
                            % if the first block is splited
                            if(angle_i == 1)
                               
                                ite_i = bins_num;
                                while(ite_i>0 && ( angle_data(ite_i) ~= 0) )
                                   ite_i = ite_i-1;
                                end
                                start_ind = ite_i;
                                
                                while(angle_i<bins_num && ( angle_data(angle_i) ~=0 ))
                                   angle_i = angle_i+1;
                                end
                                end_ind = angle_i;
                                
                                
                                gap = end_ind + ( bins_num - start_ind);
                                
                                select_angle_start = (start_ind + 1/4 * gap) / bins_num * (2*pi);
                                select_angle_end = (end_ind - 1/4 * gap) / bins_num * (2*pi);
%                                 select_angle_start = (start_ind +3 ) / bins_num * (2*pi);
%                                 select_angle_end = (end_ind -3 ) / bins_num * (2*pi);
                                
                                select_ind_1 = find( theta >= select_angle_start & theta <= bins_num );
                                select_ind_2 = find( theta >= 0 & theta <= select_angle_end);
                                select_ind = [select_ind_1;select_ind_2];

                                                          
                            else
                                if(sector <= sector_num)
                                    start_ind = angle_i-1;
                                    while(angle_i<=bins_num && ( angle_data(angle_i) ~=0 ) )
                                       angle_i = angle_i+1;
                                    end

                                    end_ind = angle_i;
                                    
                                    gap = end_ind - start_ind;
                                    
                                    select_angle_start = (start_ind + 1/4 * gap) / bins_num * (2*pi);
                                    select_angle_end = (end_ind - 1/4 * gap) / bins_num * (2*pi);
%                                     select_angle_start = (start_ind+3 ) / bins_num * (2*pi);
%                                     select_angle_end = (end_ind-3 ) / bins_num * (2*pi);

                                    select_ind = find(theta>=select_angle_start & theta<=select_angle_end);
                                    
                                    
                                end
                            end
                            
                            if (sector <= sector_num)
                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    %
                                    % Select 4 corners of the sector 
                                    %
                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                            structure = [region_x(select_ind(:)) region_y(select_ind(:))];
                                            %plot(structure(:,1), structure(:,2), 'g+');

                                            dist_set = distance(select_ind(:));
                                            select_num = length(dist_set);
                                            %corners: record the corner
                                            %         index
                                            corners = zeros(4,1);
                                            
                                            %%% get the break points
                                            outboudary_ind = select_ind + 1;
                                            temp = find(ismember(outboudary_ind,select_ind)==0);
                                            clear outboudanry_ind;
                                            
                                            outboudary_ind = select_ind - 1;
                                            temp = [temp;find(ismember(outboudary_ind,select_ind) == 0)];
                                            clear outboudanry_ind;
                                            
                                            temp = unique(temp);
                                            len_temp = length(temp);
                                            %%% get corners
                                            if (len_temp > 4)
                                                %%% index difference
                                                
                                                ordered_temp = sort(temp);
                                                temp_ind = 1:len_temp;
                                                temp_shift_ind = temp_ind - 1;
                                                temp_shift_ind(temp_shift_ind==0) = len_temp;
                                                ordered_temp_diff = abs( ordered_temp(temp_ind) - ordered_temp(temp_shift_ind) );

                                                

                                                if( length(ordered_temp_diff(ordered_temp_diff==1))==2  )
                                                    [~,temp_order] = sort(ordered_temp_diff);
                                                    % clockwise
                                                    it =  temp_order(1);
                                                    corners(1,1) = ordered_temp(it);
                                                    it =  temp_order(2);
                                                    corners(2,1) = ordered_temp(it);

                                                    % anti-clockwise
                                                    temp_shift_ind = temp_ind + 1;
                                                    temp_shift_ind(temp_shift_ind==len_temp+1) = 1;
                                                    ordered_temp_diff = abs( ordered_temp(temp_ind) - ordered_temp(temp_shift_ind) );

                                                    [~,temp_order] = sort(ordered_temp_diff);

                                                    it =  temp_order(1);
                                                    corners(3,1) = ordered_temp(it);
                                                    it =  temp_order(2);
                                                    corners(4,1) = ordered_temp(it);
                                                else
                                                    %%% when index is too close
                                                    corner_dots = structure(temp,:);
                                                    relative_dist = zeros(len_temp,1);

                                                    % scores-> pick corners
                                                    scores = zeros(len_temp,1);
                                                    for corner_i = 1:len_temp
                                                        relative_dist = sum( (corner_dots - repmat(corner_dots(corner_i,:),len_temp,1)).^2,2 );
                                                        scores = scores + relative_dist;
                                                    end

                                                    [~,temp_order] = sort(scores,'descend');

                                                    it =  temp_order(1);
                                                    corners(1,1) = temp(it);
                                                    it =  temp_order(2);
                                                    corners(2,1) = temp(it);
                                                    it =  temp_order(3);
                                                    corners(3,1) = temp(it);
                                                    it =  temp_order(4);
                                                    corners(4,1) = temp(it);
                                                end
                                            else
                                                if(len_temp == 4)
                                                    corners = temp;
                                                else
                                                    dist_set = distance(select_ind(:));
                                                    select_num = length(dist_set);

                                                    %
                                                    % right direction to compute
                                                    % distance diff
                                                    %
                                                    diff_ind = (1:1:select_num) +1;
                                                    diff_ind(select_num) = 1;
                                                    dist_diff = abs( dist_set(1:select_num) - dist_set(diff_ind) );
                                                    [~,dist_order] = sort(dist_diff,'descend');

                                                    it =  dist_order(1);
                                                    corners(1,1) = it;
                                                    it =  dist_order(2);
                                                    corners(2,1) = it;

                                                    %
                                                    % inverse to computer distance
                                                    % diff
                                                    %
                                                    diff_ind = (1:1:select_num) -1;
                                                    diff_ind(1) = select_num;
                                                    dist_diff = abs( dist_set(1:select_num) - dist_set(diff_ind) );
                                                    [~,dist_order] = sort(dist_diff,'descend');

                                                    it =  dist_order(1);
                                                    corners(3,1) = it;
                                                    it =  dist_order(2);
                                                    corners(4,1) = it; 
                                                end
                                            end
                                            % sort points into inner and outer    
                                            %
                                            [~,dots_order] = sort( dist_set(corners) );
                                            %inner dots
                                            t_i = dots_order(1);
                                            inner_dot_1 = select_ind(corners(t_i));
                                            t_i = dots_order(2);
                                            inner_dot_2 = select_ind(corners(t_i));
                                            %outer dots
                                            t_i = dots_order(3);
                                            outer_dot_1 = select_ind(corners(t_i));
                                            t_i = dots_order(4);
                                            outer_dot_2 = select_ind(corners(t_i));
                                            
%                                             plot( region_x(inner_dot_1),region_y(inner_dot_1), 'b+' )
%                                             plot( region_x(inner_dot_2),region_y(inner_dot_2), 'b+' )
%                                             plot( region_x(outer_dot_1),region_y(outer_dot_1), 'ro' )
%                                             plot( region_x(outer_dot_2),region_y(outer_dot_2), 'ro' )
                                            
                                        %
                                        % select inner points index
                                        %
                                        if( max(outer_dot_1,outer_dot_2) < max(inner_dot_1,inner_dot_2) && min(outer_dot_1,outer_dot_2) > min(inner_dot_1,inner_dot_2) )
                                            out_ind = select_ind < min(inner_dot_1,inner_dot_2) | select_ind > max(inner_dot_1,inner_dot_2);
                                        else
                                            out_ind = select_ind > min(inner_dot_1,inner_dot_2) & select_ind < max(inner_dot_1,inner_dot_2);
                                        end
                                        
                                        if(sector == 1)
                                            conic.c2 = [structure(out_ind,1) structure(out_ind,2)];
                                        else
                                            conic.c2 = [conic.c2 ; structure(out_ind,1) structure(out_ind,2)];
                                        end
                                        %hold on; plot(structure(out_ind,1), structure(out_ind,2), 'r+');

                                        %
                                        % select outer points index
                                        %
                                        clear out_ind;
                                        
                                        if(min(outer_dot_1,outer_dot_2) < min(inner_dot_1,inner_dot_2) && max(outer_dot_1,outer_dot_2) > max(inner_dot_1,inner_dot_2) )
                                            out_ind = select_ind < min(outer_dot_1,outer_dot_2) | select_ind > max(outer_dot_1,outer_dot_2);
                                        else
                                            out_ind = select_ind > min(outer_dot_1,outer_dot_2) & select_ind < max(outer_dot_1,outer_dot_2);
                                        end
                                        
                                        if(sector == 1)
                                            conic.c3 = [structure(out_ind,1) structure(out_ind,2)];
                                        else
                                            conic.c3 = [conic.c3 ; structure(out_ind,1) structure(out_ind,2)];
                                        end
                                        %hold on; plot(structure(out_ind,1), structure(out_ind,2), 'go');
                                    
                            end


                            sector = sector + 1;
                            
                        end
                        
                        angle_i = angle_i + 1;
                    end
                    
                    if( (i-start_result)==1)
                        Results.id = ellipse_id;
                        Results.conic = conic;
                    else
                        Results.id = [Results.id;ellipse_id];
                        Results.conic = [Results.conic;conic];
                    end
end

