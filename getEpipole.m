function [et e] = getEpipole(E)
[u d v] = svd(E);
et = u(:,end);
e = v(:,end);