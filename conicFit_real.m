function conic = conicFit_real(xn, num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions: conic fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalize the conics before fit to a 3x3 metrix
T = normalize( [xn.c1' xn.c2' xn.c3'], 0);
for i=1:num
    eval(['conic.C' int2str(i) '=conicFit_RANSAC(xn.c' int2str(i) ', T, 16, 1.62);']);
end





