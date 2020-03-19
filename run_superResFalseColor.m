
function run_superResFalseColor
  close all; clear; rng(1);

  datacase = 1;

  [M,C,sigma,lambda] = loadPanSharpenData( datacase );

  outDir = './output/falseColor';
  
  if ~exist( outDir, 'dir' ), mkdir( outDir ); end
  

  nChannels = size( C, 3 );
  superRes = cell( 1, 1, nChannels );
  parfor ch = 1:nChannels
    %sM = size( M );
    %sr0 = imresize( C(:,:,ch), sM, 'bilinear' );
    superRes{ ch } = diChromSuperRes( M, C(:,:,ch), sigma, 'lambda', lambda, ... %'nWs', sr0, ...
      'verbose', true );
  end
  colorSuperRes = cell2mat( superRes );


  showScale = 5;
  figure; imshowscale( M, showScale );  titlenice( 'Monochrome Image' );
  colorImg = imColorResize( C, [size(M) 3], 'nearest' );
  figure; imshowscale( colorImg, showScale );  titlenice( 'Original Color' );  
  figure; imshowscale( colorSuperRes, showScale );  titlenice( 'Sharpened Color' );

  if exist( 'outDir', 'var' ) && numel( outDir ) > 0
    imwrite( M, [ outDir, '/monochrome.jpg'] );
    imwrite( colorImg, [ outDir, '/origColor.jpg'] );
    imwrite( colorSuperRes, [ outDir, '/origColor.jpg'] );
  end


  switch nChannels
    case 3
      chWeights = [ 0.3 0.6 0.1 ];
      chTitles = { 'red', 'green', 'blue' };
    case 4
      %chWeights = [ 0.23 0.24 0.11 0.42 ];  % quickbird
      chWeights = [ 0.25 0.25 0.25 0.25 ];  % ikonos
      chTitles = { 'red', 'green', 'blue', 'ir' };
  end

  synMono = zeros( size(M) );
  for ch=1:nChannels
    synMono = synMono + chWeights(ch) * colorSuperRes(:,:,ch);
  end
  figure;  imshowscale( synMono, 5 );  titlenice('Syn Monochrome')

  for falseCh = 1:3
    %falseColor = synMono;
    falseColor = M;
    for otherCh = 1:nChannels
      if otherCh == falseCh, continue; end
      falseColor = falseColor - chWeights(otherCh) * colorSuperRes(:,:,otherCh);
    end
    falseColorSR = colorSuperRes;
    falseColorSR(:,:,falseCh) = 0;
    falseColorSR(:,:,falseCh) = falseColor / chWeights( falseCh );
    figure; imshowscale( falseColorSR, showScale );
    titlenice([ 'Sharpened False ', chTitles{falseCh}, ' Color' ]);
  end

end

