
function srImg = diChromSuperRes( highResImg, lowResImg, sigma, varargin )
  % srImg = superResC13( highResImg, imgLowRes )
  %
  % Inputs:
  % highResImg - a 2D array representing a grayscale of the high resolution proton image
  %
  % Optional Inputs:
  % N - number of iterations for fista
  %
  % Outputs:
  % srPyr - the super resolved pyruvate image
  %
  % Written by Nicholas Dwork, Copyright 2019

  p = inputParser;
  p.addRequired( 'highResImg', @isnumeric );
  p.addRequired( 'lowResImg', @isnumeric );
  p.addParameter( 'lambda', 1d3, @isnumeric );  % Tikhonov regularization parameter
  p.addParameter( 'N', 30, @isnumeric );
  p.addParameter( 'ws', [], @isnumeric );
  p.addParameter( 'optAlg', 'fista_wLS', @(x) true );
  p.addParameter( 'verbose', 1, @(x) isnumeric(x) || islogical(x) );
  p.parse( highResImg, lowResImg, varargin{:} );
  lambda = p.Results.lambda;
  N = p.Results.N;
  ws = p.Results.ws;
  optAlg = p.Results.optAlg;
  verbose = p.Results.verbose;

  lambda = lambda / numel( lowResImg );

  if numel( ws ) > 0, ws = repmat( ws, [1 1 2] ); end

  sHighRes = size( highResImg );
  srImg0 = imresize( lowResImg, sHighRes, 'bilinear' );

  sLowRes = size( lowResImg );
  xs = ones(sHighRes(1),1) * (1:sHighRes(2));
  ys = (1:sHighRes(1))' * ones(1,sHighRes(2));
  xqs = imresize( xs, sLowRes, 'nearest' );
  yqs = imresize( ys, sLowRes, 'nearest' );

  function out = f( x )
    smoothX = smoothImg( x, 'gaussian', sigma );
    out = imresize( smoothX, sLowRes, 'nearest' );
  end

  D = @(x) computeGradient( x );

  function [fx, wDx] = applyA( x )
    fx = f( x );
    Dx = D( x );
    if numel( ws ) > 0 
      wDx = lambda * Dx .* ws;
    else
      wDx = lambda * Dx;
    end
  end

  function out = fAdj( x )
    upSamp = zeros( sHighRes );
    for i=1:sLowRes(2)
      for j=1:sLowRes(1)
        upSamp( yqs(j,i), xqs(j,i) ) = x( j, i );
      end
    end
    out = smoothImg( upSamp, 'gaussian', sigma, 'op', 'transp' );
  end

  DT = @(x) computeGradient( x, 'op', 'transp' );

  function out = applyAdjA( x1, x2 )
    out1 = fAdj( x1 );
    wx2 = x2;
    if numel( ws ) > 0, wx2 = wx2 .* ws; end
    out2 = lambda * DT( wx2 );
    out = out1 + out2;
  end

  sPyrImg = size(lowResImg);
  nLowRes = numel(lowResImg);
  function out = A( x, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      x = reshape( x, sHighRes );
      [ fx, wDx ] = applyA( x );
      out = [ fx(:); wDx(:); ];
    else
      x1 = reshape( x(1:nLowRes), sPyrImg );
      x2 = reshape( x(nLowRes+1:end), [ sHighRes 2 ] );
      out = applyAdjA( x1, x2 );
      out = out(:);
    end
  end

  Dw = D( highResImg );
  if numel( ws ) > 0, Dw = Dw .* ws; end
  b = [ lowResImg(:); lambda * Dw(:); ];

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

  %proxth = @(x,t) min( max( x, 0 ), 1 );
  proxth = @(x,t) max( x, 0 );

  %[check,err] = checkAdjoint( srImg0, @f, 'fAdj', @fAdj );
  %[check,err] = checkAdjoint( srImg0, D, 'fAdj', DT );
  %[check,err] = checkAdjoint( srImg0, @A );


  if strcmp( optAlg, 'projSubgrad' )
    t = 1.0;
    proj = @(x) min( max( x, 0 ), 1 );
    srImg = projSubgrad( srImg0(:), @gGrad, proj, 't', t, 'verbose', verbose );
  elseif strcmp( optAlg, 'fista' )
    [srImg,objectiveValues] = fista( srImg0(:), @g, @gGrad, proxth, ...
      'h', @h, 'N', N, 'verbose', verbose, 'printEvery', 10 );
  else
    [srImg,objectiveValues] = fista_wLS( srImg0(:), @g, @gGrad, proxth, ...
      'h', @h, 'N', N, 'verbose', verbose, 'printEvery', 10 );
  end
  srImg = reshape( srImg, sHighRes );


  function out = objectiveValue( x )
    [fx, Dx] = applyA( x );
    Dx(:,:,1) = Dx(:,:,1) .* srImg0;
    Dx(:,:,2) = Dx(:,:,2) .* srImg0;
    out1 = norm( fx(:) - lowResImg(:), 2 ).^2;
    out2 = lambda * norm( Dx(:) - Dw(:), 2 ).^2;
    out = 0.5 * ( out1 + out2 );
  end

  if verbose ~= 0 && numel( gcp( 'nocreate' ) ) == 0
    disp([ 'Final objectiveValue: ', num2str( objectiveValue( srImg ) ) ]);
    if exist( 'lsqrFlag', 'var' )
      disp([ 'lsqr Flag: ', num2str(lsqrFlag) ]);
    end
    srPyrNearest = imresize( lowResImg, sHighRes, 'nearest' );
    figure; imshowscale( srPyrNearest, 5 );  titlenice( 'Initial Guess' );
    figure; imshowscale( srImg, 5 );  titlenice( 'Super Res Pyruvate' );
    figure; imshowscale( highResImg, 5 );  titlenice( 'Proton Image' );

    if exist( 'objectiveValues', 'var' )
      figure; plotnice( objectiveValues ); titlenice( 'FISTA Objective Values' );
    end
  end
end

