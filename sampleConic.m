function X = sampleConic(C, N)

% C = V*D*V'
[V, D] = eig(C);
D = diag(D);

% move the -ve eigenvalue/eigenvector to the end
negs = length(find(D < 0));
if negs > 1
  D = -D;
end
if D(3) > 0
  neg = find(D < 0);
  aux = D(3);
  D(3) = D(neg);
  D(neg) = aux;
  aux = V(:,3);
  V(:,3) = V(:,neg);
  V(:,neg) = aux;
end

% lengths of the axes of the conic
rx = sqrt(-D(3)/D(1));
ry = sqrt(-D(3)/D(2));

X = [];
for i = 1:N
  theta = 2*pi*(i-1)/N;
  X = [X [rx*cos(theta)
          ry*sin(theta)
          1            ]];     
end

x = V*X;

px = x(1,:)./x(3,:);
py = x(2,:)./x(3,:);
X = [px' py' ones(length(px),1)];
if x(3,1) < 0
    X = X(length(px):-1:1,:);
end
return


