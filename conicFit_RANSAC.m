function Conic = conicFit_RANSAC(p, T, N, magicvalue)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions: conic fitting by RANSAC
%T is the normalize matrix, N is the number of sample, magicvalue is 1.62
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
len = size(p,1);
std_p = 1/T(1,1);
num=1000; % modify
idx = len*rand(N,num);
idx = ceil(idx);
max = -1;
t = T*[p';ones(1,size(p,1))];
t = t./repmat(t(3,:),3,1);
AE = [t(1,:).*t(1,:);2*t(1,:).*t(2,:);t(2,:).*t(2,:);
      2*t(1,:);2*t(2,:);ones(1,size(p,1))]';
if (size(AE, 1)<5)
    disp('not enough points for ellipse fitting!');
    return;
end
for i=1:num
    j = idx(:,i);
    A = AE(j,:);
    [u, d, v] = svd(A);
    v = v(:,6);
    if v(2)^2>v(1)*v(3)
        continue;
    end
    C = [v(1)  v(2) v(4)
        v(2)   v(3) v(5)
        v(4)   v(5) v(6)];
    [V,D]=eig(C);
    D = [D(1,1) D(2,2) D(3,3)];
    D = sort(D);
    if D(1)<0 & D(2) >0 & D(3)>0
        C = C;
    else
        C = -C;
    end
    C = T'*C*T;
    C = C/C(3,3);
    
    d = [];
    for k=1:len
        x = [p(k,1);p(k,2);1];
        d = [d sqrt(abs(x'*C*x))];
    end
    ii = find(d < sqrt(magicvalue)*std_p/4);
    num = size(ii,2);
    if num > max;
        max = num;
        Conic = C;
    end
end
return;



