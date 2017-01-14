function T = normalize(x,homo)

%    T = normalize(x)
%    T = transformation induced by data normalization
if homo
    dim = size(x,1)-1;
else
    dim = size(x,1);
    x = [x;ones(1,size(x,2))];
end
  % Transform taking X's centroid to the origin

  Ttrans = [ eye( dim ), -mean( x(1:dim,:)' )' ; zeros( 1, dim ), 1 ];

  % Calculate appropriate scaling factor

  x = Ttrans * x;
  lengths = sqrt( sum( x(1:dim,:).^2 ));
  s = sqrt(dim) / mean(lengths);

  % Transform scaling x to an average length of sqrt(2)

  Tscale = [ s * eye(dim), zeros(dim,1) ; zeros(1,dim), 1 ];

  % Compose the transforms

  T = Tscale * Ttrans;