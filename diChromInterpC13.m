
function srPyr = diChromInterpC13( protonImg, metaboliteImg, sigma, varargin )
  % srPyr = superResC13( protonImg, pyrImg )
  %
  % Inputs:
  % protonImg - a 2D array representing a grayscale of the high resolution proton image
  %
  % Optional Inputs:
  % N - number of iterations for fista
  %
  % Outputs:
  % srPyr - the super resolved pyruvate image
  %
  % Written by Nicholas Dwork, Copyright 2019

  sProtonImg = size( protonImg );
  srMetabolite0 = imresize( metaboliteImg, sProtonImg, 'bilinear' );
  %srMetabolite0 = imresize( metaboliteImg, sProtonImg, 'bicubic' );

  srPyr = diChromInterp( protonImg, metaboliteImg, sigma, 'ws', srMetabolite0, varargin{:} );
end

