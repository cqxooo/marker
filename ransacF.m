function [F, bestidx] = ransacF(left, right) 


% RANSAC algorithm for the estimation of the essential matrix
% Input  - A -> 3xn set of homogeneous points in image A
%        - B -> 3xn set of homogeneous points in image B
% Output - 3x3 fundamental matrix F


%%%%%%%%%% Parameters Set UP %%%%%%%%%%
    % Threshold
    t = 0.0075;
    
    % Number of Maximun Loop
    N = 10000;
    
    % Total Population
    num = size(left.point,2);
    
    %Ratio of Outlier in the total population
    ratio = 0.05;
    
    % Acceptable number of inlier
    Accept_N = round((1-ratio)*num);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Best scores
    bestscore = 0;
    bestidx = [];
    
    % First normalize the points
    T1 = normalize(left.point,1);
    p1 = T1*left.point;
    p1 = p1./repmat(p1(3,:),3,1);
    T2 = normalize(right.point,1);
    p2 = T2*right.point;
    p2 = p2./repmat(p2(3,:),3,1);
    
    for i =1:N
        % Select 8 random points                    
        idx = randperm(num,8);            
        
        % Compute the fundamental matrix
        Ft = eightpoint(p1(:,idx),p2(:,idx));
        
        X2tFX1 = zeros(1,size(p1,2));
            for n = 1:size(p1,2)
                X2tFX1(n) = p2(:,n)'*Ft*p1(:,n);
            end

            FX1 = Ft*p1;
            FtX2 = Ft'*p2;     

            % Evaluate distances
            d =  X2tFX1.^2 ./ (FX1(1,:).^2 + FX1(2,:).^2 + FtX2(1,:).^2 + FtX2(2,:).^2);
             
            inliersL = find(abs(d) < t);     % Indices of inlying points
            
            
            % Check this model (size of the inliers set)
            if(length(inliersL)>bestscore)
                bestscore = length(inliersL);
                bestidx = inliersL;
            end;
            if bestscore >= Accept_N,break,end     
    end;
    
    % Reestimate F based on the inliers of the best model only 
%     id = left.id(bestidx);
    F = eightpoint(p1(:,bestidx),p2(:,bestidx));
    F = T2'*F*T1;
%     left.id = id;
%     right.id = id;
%     left.point = left.point(:,bestidx);
%     right.point = right.point(:,bestidx);
    
