
function run_codisData
  close all;  clear;

  showScale = 5;
  nPts = 10;
  sigma = 1.0;
  lambda = 1d3;

  a=[];  b=[];
  load( '/Users/nicholasdwork/.DellEMS/Vols/disk2s1/feromoxytol/high_res.mat' );
  a = double( a );  a = a / max( a(:) );
  b = double( b );  b = b / max( b(:) );

  highResImg = rot90( squeeze( a(128,11:248,:) ) );
  lowResImg = rot90( squeeze( b(48,:,:) ) );

  dsHighResImg = imresize( highResImg, size( lowResImg ) );
  figure; imshowscale( dsHighResImg, showScale*2 );
  titlenice( 'Downsampled high res img' );
  dsHighResPts = getFeaturesFromImg( nPts, showScale*2 );
  labelImgPts( dsHighResPts * showScale*2 );


  figure; imshowscale( lowResImg, showScale*2 );
  titlenice( 'Low res img' );
  lowResPts = getFeaturesFromImg( nPts, showScale*2 );
  labelImgPts( lowResPts * showScale*2 );

  H21 = homographyFromPts2D( lowResPts, dsHighResPts );
  projLowRes = projectImage( lowResImg, H21 );
  figure; imshowscale( projLowRes, showScale*2 );
  titlenice( 'Projected Image' );

  figure; imshowscale( imresize( lowResImg, size(highResImg), 'nearest' ), showScale );
  titlenice( 'Low Res Img' );


  %srPyr = superResC13( highResImg, projLowRes, sigma, 'lambda', lambda );
  srPyr = diChromSuperRes( highResImg, projLowRes, sigma, 'lambda', lambda );

  cube = cat( 3, imresize(projLowRes,size(srPyr)*3,'nearest'), imresize(srPyr,3,'nearest') );
  cube2Gif( cube, 'forCodi.gif', 1.0 );
end
