function cp = getIntersectConic(C1, C2)
% get the intersection between line l and coinc C. 
% the result should be two points in general.

conic1 = sym(C1);
conic2 = sym(C2);

pt = sym('[x y 1]');
p = sym('[x; y; 1]');
eqn1 = pt*conic1*p;
eqn2 = pt*conic2*p;
r = solve(eqn1, eqn2, 'x','y');
num = size(r.x,1);
xx = [subs(r.x) subs(r.y) ones(num,1)]';
a = real(xx(:,1));
b = imag(xx(:,1));
c = a + b*(-1j);
cp = [c conj(c)];
if cp(1,1)>=cp(1,2)
    cp = cp(:,[2 1]);
end