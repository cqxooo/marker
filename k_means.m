%
%
%   Date        : 2016-6-4
%   Objective   : group the convex hull points into four groups
%   Input       : a set of points
%   Author      : lampson
%
%

function out_dots = k_means(dots)
    k = 4;
    
    % initial group centroid and labels
    dot_length = length(dots);
    group_centroid = zeros(k,2);
    group_centroid(1:k,:) = dots(1:k,:);
%     for i = 1:k
%         group_centroid(i,:) = dots(1+floor( (dot_length/(k-1)-1) * (i-1) ),:);
%     end
    hw_rect = floor( max(dots)+10 );
    out_dots = process(dots,dot_length,group_centroid,k);
    
    mask_a = poly2mask(dots(:,2),dots(:,1),hw_rect(:,1),hw_rect(:,2));
    area_origin = sum( mask_a(:)==1 );
    mask_b = poly2mask(out_dots(:,2),out_dots(:,1),hw_rect(:,1),hw_rect(:,2));
    area_out = sum( mask_b(:)==1 );
    
    area_ratio = area_out/area_origin;
    
    while(area_ratio<0.8 || area_ratio>1.2)
        temp = floor( rand(k,1)*dot_length+1 );
        group_centroid(1:k,:) = dots(temp(1:k),:);
        
        out_dots = process(dots , dot_length , group_centroid , k);
    
        mask_a = poly2mask(dots(:,2),dots(:,1),hw_rect(:,1),hw_rect(:,2));
        area_origin = sum( mask_a(:)==1 );
        mask_b = poly2mask(out_dots(:,2),out_dots(:,1),hw_rect(:,1),hw_rect(:,2));
        area_out = sum( mask_b(:)==1 );
        
        area_ratio = area_out/area_origin;
    end
%     figure;
%     imshow(mask_b);
end


function out_dots = process(dots , dot_length , group_centroid, k )
    % dist for each point to k groups
    dist_group = zeros(dot_length,k);
    
    for i = 1:k
        diff = dots - repmat( group_centroid(i,:), dot_length, 1);
        dist = sum( diff.^2,2 );
        
        dist_group(:,i) = dist;
    end
    
    [~,group_label] = min(dist_group,[],2);
    last_label = zeros(dot_length,1);
    
    % iterate to group all dots
    stop_flag = norm(group_label-last_label);
    while(stop_flag > 0)
        
        for i = 1:k
           group_centroid(i,:) = mean(dots( group_label == i,: ) , 1);
           diff = dots - repmat( group_centroid(i,:), dot_length, 1);
           dist = sum( diff.^2,2 );
        
           dist_group(:,i) = dist;
        end
        
        last_label = group_label;
        [~,group_label] = min(dist_group,[],2);
        
        for i = 1:k
           if(sum(group_label==i) == 0)
               diff = dots - repmat( group_centroid(i,:), dot_length, 1);
               dist = sum( diff.^2,2 );
               %[~,minID] = min(dist);
               [~,maxID] = max(dist);
               group_label(maxID) = i;
           end
        end
        
        stop_flag = norm(group_label-last_label);
    end
    
    
    % remove the most dots in a group util only one left in each group
    for i = 1:k
       while(sum(group_label==i) ~= 1)
           group_ind = find(group_label == i);
           
           group_centroid(i,:) = mean(dots( group_ind,: ) , 1);
           diff = dots( group_ind,: ) - repmat( group_centroid(i,:), length(group_ind), 1);
           dist = sum( diff.^2,2 );
           
%            [~,centerID] = min(dist);
%            group_centroid(i,:) = dots(centerID,:);
%            diff = dots( group_ind,: ) - repmat( group_centroid(i,:), length(group_ind), 1);
%            dist = sum( diff.^2,2 );
           
           [~,maxID] = max(dist);
           group_label( group_ind(maxID) ) = -1;
       end
    end
    
    k=4;
    out_dots = zeros(k,2);
    for i = 1:k
        out_dots(i,:) = dots(group_label == i,:);
    end
end