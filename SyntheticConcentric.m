function SyntheticConcentric
% Synthetic concentric circle experiment for calibration

clc;
close all;
lss = [];
lhs = [];
cps = [];
vxs = [];
%--------------------------------------------------------
% To change the input paramters, only modify this part.
% Initialization of camera intrinsic parameters
imgSize = [640,480];
%imgSize = [6000; 6000];
% f       = 800;

f       = 800;
s = 0;
alpha   = 1.5;
pp = [imgSize(1), imgSize(2)]/2;
K = [alpha*f s pp(1);
    0        f pp(2);
    0        0  1    ];
K0 = K;
% axis([-10000,10000,-10000,10000]);
%Initialization of camera external parameters
rx = 15*pi/180;
ry = -15*pi/180;
rz = -5*pi/180;
t1 = [-10 -5 -70]';
 
Rz = [cos(rz) -sin(rz) 0;
      sin(rz) cos(rz) 0;
      0       0       1];
  
Ry = [cos(ry) 0 sin(ry);
      0       1     0;
      -sin(ry) 0 cos(ry)];
  
Rx = [1      0       0;
      0    cos(rx)  -sin(rx);
      0    sin(rx)  cos(rx)];
  
R1 = Rz*Ry*Rx;
img(1).P = K*R1*[eye(3) -t1];%projection matrix R is the ratation of camera, eye(3)is the ratat between camera coordinate and world coordinate

rx = -5*pi/180;
ry = 5*pi/180;
rz = -5*pi/180;
t2 = [-5 -1 -1]';

Rz = [cos(rz) -sin(rz) 0;
      sin(rz) cos(rz) 0;
      0       0       1];
  
Ry = [cos(ry) 0 sin(ry);
      0       1     0;
      -sin(ry) 0 cos(ry)];
  
Rx = [1      0       0;
      0    cos(rx)  -sin(rx);
      0    sin(rx)  cos(rx)];
  
R2 = Rz*Ry*Rx;
img(2).P = K*R2*[eye(3) -t2]*[R1,-R1*t1;0,0,0,1];%projection matrix R is the ratation of camera, eye(3)is the ratat between camera coordinate and world coordinate

%P = KR[I|T];
% the concentric circle center and radius
r1 = 2;
r2 = 3;
center = [0 0 0;0 8 4;0 15 4;0 22 2;
          8 0 6;8 8 8;8 15 8;8 22 2]';


% cp
cp0 = [];
% cp = P*[1 0 sqrt(-1) 0]';%cp circule points
% cp0 = [cp0 cp./cp(3)]; 
% cp = P*[1 0 -sqrt(-1) 0]';
% cp0 = [cp0 cp./cp(3)]; 
for i=1:length(img)
    cp = img(i).P*[1 sqrt(-1) 0 0]';%cp circule points
    cp0 = [cp0 cp./cp(3)]; 
    cp = img(i).P*[1 -sqrt(-1) 0 0]';
    cp0 = [cp0 cp./cp(3)]; 
    l = cross(cp0(:,1),cp0(:,2));
    lh0 = -sign(l(3))*l/norm(l(1:2));
    img(i).cp = cp0;
    img(i).lh = lh0;
end
for i=1:length(img)
    cen = img(i).P*[center;ones(1,8)];
    cen0 = [cen(1,:)./cen(3,:);cen(2,:)./cen(3,:);cen(3,:)./cen(3,:)];
    img(i).cen = cen0;
end
% generate conics: numConic = 1; there may have multiple concentric circles
numConic = size(center,2);

for imLoop = 1:length(img)
    for i = 1:numConic
        % generate conics
        c = center(:,i);
        Q1 = [eye(3) -c;
            -c'    c'*c-r1*r1];
        
        Q2 = [eye(3) -c;
            -c'    c'*c-r2*r2];
        iC = img(imLoop).P*inv(Q1)*img(imLoop).P';
        C1 = inv(iC);
        conic(imLoop,i).C1 = C1./C1(3,3);
    
        iC = img(imLoop).P*inv(Q2)*img(imLoop).P';
        C2 = inv(iC);
        conic(imLoop,i).C2 = C2./C2(3,3);
        fig = figure(1);
        axis ij;
        axis([0,imgSize(1),0,imgSize(2)])
        hold on;
        drawConic(conic(imLoop,i).C1);
        drawConic(conic(imLoop,i).C2);
        plot(img(imLoop).cen(1,i), img(imLoop).cen(2,i), 'r+');
        drawLine(img(imLoop).lh, 'b-');
    end
    close(fig)
end

%=============================================
N = 100;
for imLoop =1:length(img)
    for i=1:numConic
        s = sampleConic(conic(imLoop,i).C1, N)';
        edges(imLoop,i).C1 = s;
        s = sampleConic(conic(imLoop,i).C2, N)';
        edges(imLoop,i).C2 = s;
    end
end
%end of preparing synthesized conic 
for imLoop = 1:length(img)
    cps = [];
    lhs = [];
    cens = [];
    for i = 1:length(edges)
        [C1, C2] = conicFit(edges(imLoop,i));% the image chaged 180 degrees
        % %============================================
        % % Linear approaches
        [cen, v1, v2] = getEigenVectors(C1, C2);
        l = cross(v1,v2);
        lh = -sign(l(3))*l/norm(l(1:2));
        
        [cp1, cp2] = getIntersectLineConic(lh, C1);
        [cp3, cp4] = getIntersectLineConic(lh, C2);
        a = (real(cp1)+real(cp3))/2;
        if sign(imag(cp1(1))) == sign(imag(cp3(1))) 
            b = (imag(cp1)+imag(cp3))/2;
        else
            b = (imag(cp1)-imag(cp3))/2;
        end
        c = a + b*-1j;
        cp = [c conj(c)];
        
        cps = [cps cp];
        lhs = [lhs lh];         
        cens = [cens cen];  
    end
%     save test.mat
%     load test.mat
    [error, IAC, K] = calib_cps(cps, K0);
%     concen(imLoop,:) = getOrientation(conic(imLoop,:), K, lhs, cens);
    img(imLoop).K = K;
    img(imLoop).cens = cens;
    img(imLoop).id = (1:8)';
end
result = getpoint3D(img);
for i=1:length(result)
    n = length(result(i).id);
    for j=1:n-1
        for k=j+1:n
            if j==1 && k==2
               ratio = distPnts(center(:,1),center(:,2))/distPnts(result(i).point(:,1),result(i).point(:,2));
            end
            d = distPnts(result(i).point(:,j),result(i).point(:,k))*ratio; %compute distance
            gt = distPnts(center(:,j),center(:,k));%ground truth
            e = abs(gt-d)/gt*100;
            disp(['distance between ' num2str(result(i).id(j)) ' to ' num2str(result(i).id(k)) ' is ' num2str(d)]);       %show computed distance
            disp(['distance error between ' num2str(result(i).id(j)) ' to ' num2str(result(i).id(k)) ' is ' num2str(e) '%']); %show error 
        end        
    end
end








