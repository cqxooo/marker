function P = compute_p( x, Xw, K )

  % Estimate the camera projection P best fitting
  % 3D points X and observed image points x
if size(Xw,1) == 3
    Xw = [Xw;ones(1,size(Xw,2))];
end
% if nargin == 3
%     x = inv(K)*x;
%     x = x./repmat(x(3,:),3,1);
% end
      

  % Compute isotropic scaling transforms to normalize x and Xw

  T = normalize( x,1 );
  U = normalize( Xw,1 );

  xn = T * x;
  Xwn = U * Xw;

  % Get the linear estimate of P

  Plin = dlt( xn, Xwn );
 %Refine P by non-linear
 %     options = optimset('lsqnonlin');
    options = optimset('Algorithm',{'levenberg-marquardt', 0.01});
    options = optimset(options, 'display', 'iter');
    options = optimset(options,'TolFun',1e-20);
    options = optimset(options,'TolX',1e-20);
    options = optimset(options,'MaxIter',100);
    options = optimset(options,'MaxFunEvals',800);
% set the options for lsqnonlin
  Pnonlin = lsqnonlin(@(P)reprojection_error(xn,Xwn,P),Plin,[],[],options);


%   Denormalize

  P = inv( T ) * Pnonlin * U;  
  if nargin == 3
    P = inv(K)*P;    
  end
P = P / norm(P(3,1:3));
end
%------------------------------------------------------------------------

function P = dlt( x, Xw )

  n = size( Xw, 2 );
  if n < 6
    error( 'DLT for P requires at least 6 points' );
  end;

  A = [];

  for i = 1:n

    xi = x( 1, i );
    yi = x( 2, i );
    wi = x( 3, i );
    Xwi = Xw( :, i );

    Ai = [ 0, 0, 0, 0,   -wi * Xwi',    yi * Xwi' ;
            wi * Xwi',   0, 0, 0, 0,   -xi * Xwi' ];

    A = [ A ; Ai ];
  end;

  [U,D,V] = svd( A );

  % The solution P minimizing |Ap| subject to |p|=1 is the last column of V

  P = reshape( V(:,12), 4, 3 )';
  P = P / norm(P(3,1:3));
end
%------------------------------------------------------------------------

function err = reprojection_error( x, Xw, P )

  xtilde = P * Xw;
  xtilde = xtilde ./ repmat( xtilde(3,:), 3, 1 );
  err = sqrt( mean( sum( ( x(1:2,:) - xtilde(1:2,:) ).^2 )));
end
%------------------------------------------------------------------------