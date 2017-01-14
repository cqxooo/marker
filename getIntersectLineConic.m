function [xl, xr] = getIntersectLineConic(l, C)
% get the intersection between line l and coinc C. 
% the result should be two points in general.

conic = sym(C);
pt = sym('[x y 1]');
p = sym('[x; y; 1]');
eqn1 = pt*conic*p;
% l = sym(l');
eqn2 = l'*p;
r = solve(eqn1, eqn2, 'x','y');
num = 2;
xx = [subs(r.x) subs(r.y) ones(num,1)]';
if xx(1,1)<xx(1,2)
    xl = xx(:,1);
    xr = xx(:,2);
else
    xr = xx(:,1);
    xl = xx(:,2);
end

% fid = fopen('c:\Program Files\Wolfram Research\Mathematica\5.0\lineConicIntersection.txt','w');
% for idx = 1:3
%     fprintf(fid,'%.30f ',l(idx));
% end
% fprintf(fid,'\n');
% for idxi = 1:3
%     for idxj = 1:3
%         fprintf(fid,'%.30f ',C(idxi,idxj));
%     end
%     fprintf(fid,'\n');
% end
% fclose(fid);

% x = [124.02009907565227     233.7374842438434  1];
% xp = [166.43284480623697    392.3283036029343  1];

% x = [98.41973771129268   138.0119608256585  1];
% xp = [103.27266979615649  156.15816811253023 1];
    
    
return;