function conic = conicFit(xn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions: conic fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalize the conics before fit to a 3x3 metrix
c = cell2mat(xn);
c = c(:,[2 1]);
T = normalize( c', 0);
conic = computeConic(c, T);

