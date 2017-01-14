function [E, left, right] = ransacE(left, right, K) 


% RANSAC algorithm for the estimation of the essential matrix
% Input  - A -> 3xn set of homogeneous points in image A
%        - B -> 3xn set of homogeneous points in image B
% Output - 3x3 essential matrix E


%%%%%%%%%% Parameters Set UP %%%%%%%%%%
    % Threshold
    t = 0.00001;
    
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
    p1 = inv(K)*left.point;
    p1 = p1./repmat(p1(3,:),3,1);
    p2 = inv(K)*right.point;
    p2 = p2./repmat(p2(3,:),3,1);
    
    for i =1:N
        % Select 8 random points                    
        idx = randperm(num,8);            
        
        % Compute the essential matrix
        Et = eightpointE(p1(:,idx),p2(:,idx));
        
        X2tEX1 = zeros(1,size(p1,2));
            for n = 1:size(p1,2)
                X2tEX1(n) = p2(:,n)'*Et*p1(:,n);
            end

            EX1 = Et*p1;
            EtX2 = Et'*p2;     

            % Evaluate distances
            d =  X2tEX1.^2 ./ (EX1(1,:).^2 + EX1(2,:).^2 + EtX2(1,:).^2 + EtX2(2,:).^2);
             
            inliersL = find(abs(d) < t);     % Indices of inlying points
            
            
            % Check this model (size of the inliers set)
            if(length(inliersL)>bestscore)
                bestscore = length(inliersL);
                bestidx = inliersL;
            end;
            if bestscore >= Accept_N,break,end     
    end;
    
    % Reestimate F based on the inliers of the best model only 
    id = left.id(bestidx);
    E = eightpointE(p1(:,bestidx),p2(:,bestidx));
    
    left.id = id;
    right.id = id;
    left.point = left.point(:,bestidx);
    right.point = right.point(:,bestidx);
    
