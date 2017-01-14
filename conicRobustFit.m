function conic = conicRobustFit(xn,num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions: conic fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalize the conics before fit to a 3x3 metrix
T = normalize( [xn.c1' xn.c2' xn.c3'], 0);
for i=1:num
    eval(['conic.C' int2str(i) '=computeConic(xn.c' int2str(i) ', T);']);
end
