
function srPyr = superResC13( protonImg, pyrImg, varargin )
  % srPyr = superResC13( protonImg, pyrImg )
  %
  % Outputs:
  % srPyr - the super resolved pyruvate image
  %
  % Written by Nicholas Dwork, Copyright 2019

  ds = 9;  % ratio of proton resolution to pyr resolution
  sigma = 3;
  verbose = 1;

  p = inputParser;
  p.addParameter( 'lambda', 0, @isnumeric );  % Tikhonov regularization parameter
  p.parse( varargin{:} );
  lambda = p.Results.lambda;

  lambda = lambda / numel( pyrImg );

  sProton = size( protonImg );
  newPyrSize = floor( sProton / ds );
  newProtonSize = newPyrSize * ds;
  protonImg = protonImg( 1 : newProtonSize(1), 1 : newProtonSize(2) );
  pyrImg = pyrImg( 1 : newPyrSize(1), 1 : newPyrSize(2) );

  sProtonImg = size( protonImg );
  srPyr0 = imresize( pyrImg, sProtonImg, 'bilinear' );

  function out = f( x )
    smoothX = smoothImg( x, ds, 'gaussian', sigma );
    out = smoothX( ceil(ds/2) : ds : end, ceil(ds/2) : ds : end );  %downsample
  end

  D = @(x) computeGradient( x );

  function [fx, wDx] = applyA( x )
    fx = f( x );
    Dx = D( x );
    wDx = Dx;  % weighted Derivatives
    wDx(:,:,1) = wDx(:,:,1) .* srPyr0;
    wDx(:,:,2) = wDx(:,:,2) .* srPyr0;
    wDx = lambda * wDx;
  end

  function out = fAdj( x )
    % Adjoint of g
    us = zeros( size(pyrImg) * ds );  % upsample
    us( ceil(ds/2) : ds : end, ceil(ds/2) : ds : end ) = x;
    out = smoothImg( us, ds, 'gaussian', sigma, 'op', 'transp' );
  end

  DT = @(x) computeGradient( x, 'op', 'transp' );

  function out = applyAdjA( x1, x2 )
    out1 = fAdj( x1 );
    wx2 = x2;
    wx2(:,:,1) = wx2(:,:,1) .* srPyr0;
    wx2(:,:,2) = wx2(:,:,2) .* srPyr0;
    out2 = lambda * DT( wx2 );
    out = out1 + out2;
  end

  sPyrImg = size(pyrImg);
  function out = A( x, op )
    if nargin < 2, op = []; end
    if strcmp( op, 'transp' )
      x1 = reshape( x(1:numel(pyrImg)), sPyrImg );
      x2 = reshape( x(numel(pyrImg)+1:end), [ sProtonImg 2 ] );
      out = applyAdjA( x1, x2 );
      out = out(:);
    else
      x = reshape( x, sProtonImg );
      [ fx, wDx ] = applyA( x );
      out = [ fx(:); wDx(:); ];
    end
  end


  Dw = D( protonImg );
  Dw(:,:,1) = Dw(:,:,1) .* srPyr0;
  Dw(:,:,2) = Dw(:,:,2) .* srPyr0;
  b = [ pyrImg(:); lambda * Dw(:); ];


  function out = g( x )
    out = 0.5 * norm( A(x) - b, 2 ).^2;
  end

  function out = gGrad( x )
    out = A( A( x ), 'transp' ) - A( b, 'transp' );
  end

  function out = h( x )
    out = 0;
    if min( x ) < 0, out=Inf; end
  end

  proxth = @(x,t) max( x, 0 );

  %[check,err] = checkAdjoint( srPyr0, @f, 'fAdj', @fAdj );
  %[check,err] = checkAdjoint( srPyr0, D, 'fAdj', DT );
  %[check,err] = checkAdjoint( srPyr0, @A );


  %t = 1.0;
  %proj = @(x) max( x, 0 );
  %srPyr = projSubgrad( srPyr0(:), @gGrad, proj, 't', t, 'verbose', verbose );



  [srPyr,objectiveValues] = fista( srPyr0(:), @g, @gGrad, proxth, ...
    'h', @h, 'N', 100, 'verbose', verbose );
  %[srPyr,lsqrFlag] = lsqr( @A, b, [], [], [], [], srPyr0(:) );
  srPyr = reshape( srPyr, sProtonImg );

  function out = objectiveValue( x )
    [fx, Dx] = applyA( x );
    Dx(:,:,1) = Dx(:,:,1) .* srPyr0;
    Dx(:,:,2) = Dx(:,:,2) .* srPyr0;
    out1 = norm( fx(:) - pyrImg(:), 2 ).^2;
    out2 = lambda * norm( Dx(:) - Dw(:), 2 ).^2;
    out = 0.5 * ( out1 + out2 );
  end

%   function out = g( x )
%     tmp = f( x );
%     diff = pyrImg(:) - tmp(:);
%     out = 0.5 * norm( diff(:), 2 ).^2;
%   end
% 
%   function out = gGrad( x )
%     out = fAdj( f( x ) ) - fAdj( pyrImg );
%   end
% 
%   function out = h( x )
%     Ddiff = computeGradient( x - protonImg );
%     Ddiff1 = squeeze( srPyr0 .* Ddiff(:,:,1) );
%     Ddiff2 = squeeze( srPyr0 .* Ddiff(:,:,2) );
%     out = 0.5 * lambda * ( norm( Ddiff1(:), 2 ).^2 + norm( Ddiff2(:), 2 ).^2 );           
%   end
% 
%   function out = hGrad( x )
%     Ddiff = computeGradient( x - protonImg );
%     DTDdiff = computeGradient( Ddiff, 'op', 'transp' );
%     out = 0.5 * lambda * DTDdiff;
%   end

  disp([ 'Final objectiveValue: ', num2str( objectiveValue( srPyr ) ) ]);
  if exist( 'lsqrFlag', 'var' )
    disp([ 'lsqr Flag: ', num2str(lsqrFlag) ]);
  end
  srPyrNearest = imresize( pyrImg, sProtonImg, 'nearest' );
  figure; imshowscale( srPyrNearest, 5 );  titlenice( 'Initial Guess' );
  figure; imshowscale( srPyr, 5 );  titlenice( 'Super Res Pyruvate' );
  figure; imshowscale( protonImg, 5 );  titlenice( 'Proton Image' );
  if exist( 'objectiveValues', 'var' )
    figure; plotnice( objectiveValues ); titlenice( 'FISTA Objective Values' );
  end
end

