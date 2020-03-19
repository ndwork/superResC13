
function w = findColorChannelWeights( M, C, varargin )
  % w = findColorChannelWeights( M, C [, 'objective', objective ] )
  %
  % Find w such that M = sum w(i) * C(:,:,i)
  %
  % Inputs:
  % M - monochrome image
  % C - color image
  %
  % Optional Inputs:
  % objective - either 'L2' or 'L1'
  %
  % Written by Nicholas Dwork, Copyright 2019

  p = inputParser;
  p.addParameter( 'objective', 'L2', @(x) true );
  p.parse( varargin{:} );
  objective = p.Results.objective;

  sC = size(C);

  A = reshape( C, [ sC(1)*sC(2) sC(3) ] );
  b = M(:);

  w = A \ b;

  f = @(w) norm( A*w - b, 1 );
  if strcmp( objective, 'L1' )
    wLB = w*0;
    w = fmincon( f, w, A, b, [], [], wLB, [] );
  end
end
