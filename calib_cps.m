function [error, IAC, K] = calib_cps(cps, K0)
% Calibrate with the constraints given by circular points and the rotaion
% axis, the vanishing point along the x-axis.


% w = [w1 0 w3;
%      0 w2 w4;
%      w3 w4 w5];
coef = [];
for i = 1:size(cps, 2)
    cp = cps(:,i);
    coef = [coef;     cp(1)*cp(1) cp(2)*cp(2)   2*cp(1)  2*cp(2)        1];     % iT w i = 0
end

% normalize each row
[row, column] = size(coef);

for i = 1:row
    coef(i,:)=coef(i,:)./sqrt(sum(coef(i,:).*coef(i,:)));
end

[u, d, v0] = svd(coef);

v = sign(real(v0)).*abs(v0);
C = [v(1,end)     0       v(3,end)
        0      v(2,end)   v(4,end)
     v(3,end)  v(4,end)   v(5,end)];
C = C./ abs(C(3,3));
invC = inv(C);
[tempK, p] = chol(C);
if p ~= 0
    C = -1.*C;
    [tempK, p] = chol(C);
end

K = inv(tempK);
K = K ./K(3,3)
% C is the imaged absolute conic
 IAC = C;

eaf = abs(K(1,1)-K0(1,1))/K0(1,1);
ef = abs(K(2,2)-K0(2,2))/K0(2,2);
eu = abs(K(1,3)-K0(1,3))/K0(1,3);
ev = abs(K(2,3)-K0(2,3))/K0(2,3);

error= 100*[eaf ef eu ev]


return;