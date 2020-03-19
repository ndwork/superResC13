
function run_validateSuperRes( datacases )
  showImages = false;

  if nargin < 1
    close all; rng(1); clear;
    datacases = 1;
    showImages = true;
  end
  subOutDir = './output/output_validate';

  for datacase = datacases
    outDir = [ subOutDir, '_', num2str(datacase) ];
    mkdir( outDir );

    im_proton=0;  im_pyr=0;   %#ok<NASGU>
    [im_proton, im_pyr, ~, ~, ~, ~, ~, ds] = loadC13SuperResData( datacase );
    switch datacase
      case 1
        lambda = 8d1;
      case 2
        lambda = 8d1;
      otherwise
        error( 'You need to set lambda' );
    end

    im_proton = smoothImg( im_proton, 5 );
    im_pyr = im_pyr( ceil(ds/2) : ds : end, ceil(ds/2) : ds : end, : );

    % decimate thisProton to size of thisPyr
    falseDs = 2;
    dProton = imColorResize( im_proton, size(im_pyr), 'bilinear' );
    dPyr = imColorResize( im_pyr, round( 1/falseDs * size(im_pyr) ), 'bilinear' );

    maxProton = max( dProton(:) );
    dProton = dProton / maxProton;

    maxPyr = max( dPyr(:) );
    im_pyr = im_pyr / maxPyr;
    dPyr = dPyr / maxPyr;

    if ismatrix( dPyr )
      nChannels = 1;
    else
      nChannels = size( dPyr, 3 );
    end    
    restoredPyrs = cell( 1, 1, nChannels );
    parfor ch = 1 : nChannels
      restoredPyrs{ch} = superResC13( dProton(:,:,ch), dPyr(:,:,ch), 0.5*falseDs, ...
        'lambda', lambda, 'verbose', false );
    end
    restoredPyrs = cell2mat( restoredPyrs );

    figH = figure; showImageCube( dProton, 10, 'nImgsPerRow', 5, 'border', 1, 'borderValue', 'max' );
    titlenice( 'Restorative Proton' );
    saveas( gcf, [ outDir, '/dProton.jpg' ] );
    if showImages ~= true, close( figH ); end

    figH = figure; showImageCube( im_pyr, 10, 'nImgsPerRow', 5, ...
      'border', 1, 'borderValue', 'max' );
    titlenice( 'Original Pyruvate' );
    saveas( gcf, [ outDir, '/im_pyr.jpg' ] );
    if showImages ~= true, close( figH ); end

    figH = figure; showImageCube( restoredPyrs, 10, 'nImgsPerRow', 5, ...
      'border', 1, 'borderValue', 'max' );
    titlenice( 'Restored Pyruvate' );
    saveas( gcf, [ outDir, '/restoredPyr.jpg' ] );
    if showImages ~= true, close( figH ); end

    errSR = abs( im_pyr - restoredPyrs );
    sr0 = imColorResize( dPyr, size(dProton), 'nearest' );
    errSR0 = abs( im_pyr - sr0 );
    maxErr = max( [ errSR(:); errSR0(:); ] );

    figH = figure;  showImageCube( errSR, 10, 'nImgsPerRow', 5, 'range', [0 maxErr], ...
      'border', 1, 'borderValue', 'max' );
    titlenice( 'SR Abs Errors' );  colormap( gca, 'jet' );  colorbarnice;
    saveas( gcf, [ outDir, '/errSR.jpg' ] );
    if showImages ~= true, close( figH ); end

    figH = figure;  showImageCube( errSR0, 10, 'nImgsPerRow', 5, 'range', [0 maxErr], ...
      'border', 1, 'borderValue', 'max' );
    titlenice( 'Nearest Interp Abs Errors' );  colormap( gca, 'jet' );  colorbarnice;
    saveas( gcf, [ outDir, '/errSR0.jpg' ] );
    if showImages ~= true, close( figH ); end
  end

end

