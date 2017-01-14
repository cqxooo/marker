%
%
%   Date        : 2016-6-3
%   Objective   : Get ID of the target markers
%   Input       : circle boundary, other sectors boundary
%   Output      : ID
%   Author      : lampson
%
%

function ellipse_id = getID(sector_boundary, circle_boundary, blank_angles, id_index)

        params = fitellipse(circle_boundary(:,1),circle_boundary(:,2));
        
        % rotation matrix
        rota = [cos(params(5)) -sin(params(5));sin(params(5)) cos(params(5))];
                            
        

        % rotate and affine transformation
        test_data = sector_boundary - repmat( [params(1) params(2)], length(sector_boundary), 1);
        test_data = test_data * rota;

        affi = [params(4)/params(3) 0; 0 1];
        test_data = test_data * affi;

        out_ind = test_data + repmat( [params(1) params(2)], length(sector_boundary), 1);
        out_ind = ceil(out_ind);

        % warp
        min_x = sort(out_ind(:,1));
        min_y = sort(out_ind(:,2));

           out_ind(:,1) = out_ind(:,1) - repmat( min_x(1), length(min_x),1 ) + 1; 
           out_ind(:,2) = out_ind(:,2) - repmat( min_y(1), length(min_y),1 ) + 1;
           
           cir_point = [params(1)-min_x(1),params(2)-min_y(2)];

   


        

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
        
        
                    angle_num = 180;
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
                            temp = 0;
                            while(angle_i <= bins_num && angle_data(angle_i) ~= 0)
                                temp = temp + 1;
                                angle_i = angle_i + 1;
                            end
                            
                            if(temp > 2)
                                section_num(blank_block,1) = temp;
                                
                                %label 2-> blank
                                section_num(blank_block,2) = 1;
                            
                                % the section is splitted
                                if(angle_i >= bins_num && section_num(blank_block,2) == section_num(1,2) )
                                    section_num(1,1) = section_num(1,1) + section_num(blank_block,1);
                                    section_num(blank_block,1) = 0;
                                    section_num(blank_block,2) = 0;
                                end

                            else
                                
                                % when the blank is only small angles
                                if(blank_block > 1)
                                    blank_block = blank_block - 1;
                                    section_num(blank_block,1) =  section_num(blank_block,1) + temp;
                                else
                                    section_num(blank_block,1) =  temp;
                                    blank_block = blank_block - 1;
                                end
                            end
                        else   
                            temp = 0;
                            while(angle_i <= bins_num && angle_data(angle_i) == 0)
                                temp = temp + 1;
                                angle_i = angle_i + 1;
                            end
                            
                            if(temp > 2)
                                section_num(blank_block,1) = temp;
                                
                                %label 2-> blank
                                section_num(blank_block,2) = 2;
                            
                                % the section is splitted
                                if(angle_i >= bins_num && section_num(blank_block,2) == section_num(1,2) )
                                    section_num(1,1) = section_num(1,1) + section_num(blank_block,1);
                                    section_num(blank_block,1) = 0;
                                    section_num(blank_block,2) = 0;
                                end

                            else
                                
                                % when the blank is only small angles
                                if(blank_block > 1)
                                    blank_block = blank_block - 1;
                                    section_num(blank_block,1) =  section_num(blank_block,1) + temp;
                                else
                                    section_num(blank_block,1) =  temp;
                                    blank_block = blank_block - 1;
                                end
                            end
                        end
                        
                        %combine the two sectors if two continous sectors
                        %exist
                        if(blank_block>1)
                            if( section_num(blank_block,2) == section_num(blank_block-1,2) )
                                temp = section_num(blank_block,1);
                                section_num(blank_block,2) = 0;
                                section_num(blank_block-1,1) =  section_num(blank_block-1,1) + temp;
                                blank_block = blank_block - 1;
                            end
                        end
                    end

                    
                    region_b_num = region_num - length(find(section_num(:,1)==0));
                    region_area = section_num(1:region_b_num,:);
                    
                    
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
        [~,ind] = min(dist);
        
        ellipse_id = id_index(ind);
        
        
%         abs_diff = abs(diff);
%         [r_num,c_num] = size(abs_diff);
%         id_vote = zeros(r_num,1);
%         
%         for i = 1:c_num
%             if(sum(abs_diff(:,i)) > 0)
%                 min_value = min(abs_diff(:,i));
%                 id_vote(abs_diff == min_value) = id_vote(abs_diff == min_value) + 1;
% %               id_vote(min_row) = id_vote(min_row) + 1;
%             end
%         end
%         [~,ellipse_id] = max(id_vote); 
end