
function run_c13superRes( datacases, varargin )
  showImages = false;
  if nargin < 1, close all; rng(1); clear; end

  blurFraction = 0.5;
  gamma = 4;
  scaleIndividually = true;

  if nargin < 1
    datacases = 6;
    showImages = false;
  end

  for datacase = datacases
    if ~exist( 'showImages', 'var' ) || showImages == false, close all; end

    outDir = [ './output/output_c13superRes_', num2str(datacase) ];
    mkdir( outDir );

    im_proton=0;  im_pyr=0;  im_lac=0;  im_bic=0;  im_lpRatio=0;   %#ok<NASGU>
    [im_proton, im_pyr, im_lac, im_bic, im_lpRatio, lambda, subRegion, ds] = ...
      loadC13SuperResData( datacase );

    maxProton = max( im_proton(:) );
    im_proton = im_proton / maxProton;
    maxPyr = max( im_pyr(:) );
    im_pyr = im_pyr / maxPyr;

    if numel( im_lac ) > 1
      im_lpRatio = im_pyr ./ im_lac;
      mask = im_pyr > 0.1 * max( im_pyr(:) );
      im_lpRatio( mask==0 ) = 0;
      maxLpRatio = max( im_lpRatio(:) );
      im_lpRatio = im_lpRatio / maxLpRatio;

      maxLac = max( im_lac(:) );
      im_lac = im_lac / maxLac;      
    end

    if numel( im_bic ) > 1
      maxBic = max( im_bic(:) );
      im_bic = im_bic / maxBic;
    end

    nSlices = size( im_pyr, 3 );
    origPyrs = cell( 1, 1, nSlices );
    superPyrs = cell( 1, 1, nSlices );
    origLacs = cell( 1, 1, nSlices );
    superLacs = cell( 1, 1, nSlices );
    origLpRatios = cell( 1, 1, nSlices );
    superLpRatios = cell( 1, 1, nSlices );
    origBics = cell( 1, 1, nSlices );
    superBics = cell( 1, 1, nSlices );
    resizedProtons = cell( 1, 1, nSlices );

    parfor sliceIndx = 1 : nSlices
      disp([ 'Working on slice ', num2str( sliceIndx ) ]);
      thisProton = smoothImg( im_proton(:,:,sliceIndx), 5 );
      if numel( subRegion ) > 0
        xL = subRegion(1);  xR = subRegion(3);
        yL = subRegion(2);  yR = subRegion(4);
      else
        sProton = size( thisProton );
        xL = 1;  xR = sProton(2);
        yL = 1;  yR = sProton(1);
      end
      thisProton2Save = thisProton( yL : yR, xL : xR );
      if scaleIndividually == true, thisProton2Save = thisProton2Save / max( thisProton2Save(:) ); end
      imwrite( thisProton2Save, [ outDir, '/thisProton_', num2str(sliceIndx,'%3.3i'), '.jpg'] );

      thisPyr = im_pyr( ceil(ds/2) : ds : end, ceil(ds/2) : ds : end, sliceIndx );    %#ok<PFBNS>
      thisPyrResized = imresize( thisPyr, size(thisProton), 'nearest' );
      origPyrs{ sliceIndx } = thisPyrResized( yL : yR, xL : xR );
      thisPyr2Save = thisPyrResized( yL : yR, xL : xR );
      if scaleIndividually == true, thisPyr2Save = thisPyr2Save / max( thisPyr2Save(:) ); end
      imwrite( thisPyr2Save, [ outDir, '/thisPyr_', num2str(sliceIndx,'%3.3i'), '.jpg'] );

      thisSuperPyr = superResC13( thisProton, thisPyr, blurFraction*ds, 'lambda', lambda );
      thisSuperPyr2Save = thisSuperPyr( yL : yR, xL : xR );
      if scaleIndividually == true, thisSuperPyr2Save = thisSuperPyr2Save / max( thisSuperPyr2Save(:) ); end
      imwrite( thisSuperPyr2Save, [ outDir, '/thisSuperPyr_', num2str(sliceIndx,'%3.3i'), '.jpg'] );
      superPyrs{ sliceIndx } = thisSuperPyr( yL : yR, xL : xR );

      if numel( im_lac ) > 1
        thisLac = im_lac( ceil(ds/2) : ds : end, ceil(ds/2) : ds : end, sliceIndx );
        thisLacResized = imresize( thisLac, size( thisProton ), 'nearest' );
        origLacs{ sliceIndx } = thisLacResized( yL : yR, xL : xR );
        thisLac2Save = thisLacResized( yL : yR, xL : xR );
        if scaleIndividually == true, thisLac2Save = thisLac2Save / max( thisLac2Save(:) ); end
        imwrite( thisLac2Save, [ outDir, '/thisLac_', num2str(sliceIndx,'%3.3i'), '.jpg'] );

        thisSuperLac = superResC13( thisProton, thisLac, blurFraction*ds, 'lambda', lambda );
        thisSuperLac2Save = thisSuperLac( yL : yR, xL : xR );
        if scaleIndividually == true, thisSuperLac2Save = thisSuperLac2Save / max( thisSuperLac2Save(:) ); end
        imwrite( thisSuperLac2Save, [ outDir, '/thisSuperLac_', num2str(sliceIndx,'%3.3i'), '.jpg'] );
        superLacs{ sliceIndx } = thisSuperLac( yL : yR, xL : xR );

        thisRatio = im_lpRatio( ceil(ds/2) : ds : end, ceil(ds/2) : ds : end, sliceIndx );    %#ok<PFBNS>
        thisRatioResized = imresize( thisRatio, size( thisProton ), 'nearest' );
        origLpRatios{ sliceIndx } = thisRatioResized( yL : yR, xL : xR );
        imwrite( thisRatioResized( yL : yR, xL : xR ), ...
          [ outDir, '/thisLpRatio_', num2str(sliceIndx,'%3.3i'), '.jpg'] );

        thisSuperLpRatio = superResC13( thisProton, thisRatio, blurFraction*ds, 'lambda', lambda );
        imwrite( thisSuperLpRatio( yL : yR, xL : xR ), ...
          [ outDir, '/thisSuperRatio_', num2str(sliceIndx,'%3.3i'), '.jpg'] );
        superLpRatios{ sliceIndx } = thisSuperLpRatio( yL : yR, xL : xR );
      end

      if numel( im_bic ) > 1
        thisBic = im_bic( ceil(ds/2) : ds : end, ceil(ds/2) : ds : end, sliceIndx );
        thisBicResized = imresize( thisBic, size(thisProton), 'nearest' );
        origBics{ sliceIndx } = thisBicResized( yL : yR, xL : xR );
        thisBic2Save = thisBicResized( yL : yR, xL : xR );
        if scaleIndividually == true, thisBic2Save = thisBic2Save / max( thisBic2Save(:) ); end
        imwrite( thisBic2Save, [ outDir, '/thisBic_', num2str(sliceIndx,'%3.3i'), '.jpg'] );

        thisSuperBic = superResC13( thisProton, thisBic, 0.25*ds, 'lambda', lambda );
        thisSuperBic2Save = thisSuperBic( yL : yR, xL : xR );
        if scaleIndividually == true, thisSuperBic2Save = thisSuperBic2Save / max( thisSuperBic2Save(:) ); end
        imwrite( thisSuperBic2Save, [ outDir, '/thisSuperBic_', num2str(sliceIndx,'%3.3i'), '.jpg'] );
        superBics{ sliceIndx } = thisSuperBic( yL : yR, xL : xR );
      end

      thisResizedProton = thisProton( 1 : size(thisSuperPyr,1), 1 : size(thisSuperPyr,2) );
      resizedProtons{ sliceIndx } = thisResizedProton( yL : yR, xL : xR );
      thisResizedProton2Save = thisResizedProton( yL : yR, xL : xR );
      if scaleIndividually == true
        thisResizedProton2Save = thisResizedProton2Save / max( thisResizedProton2Save(:) );
      end
      imwrite( thisResizedProton2Save, [ outDir, '/thisResizedProton_', num2str(sliceIndx,'%3.3i'), '.jpg'] );
      
      if showImages == false, close all; end
    end

    superPyrs = cell2mat( superPyrs );
    origPyrs = cell2mat( origPyrs );
    if numel( im_lac ) > 1
      origLacs = cell2mat( origLacs );
      superLacs = cell2mat( superLacs );
      origLpRatios = cell2mat( origLpRatios );
      superLpRatios = cell2mat( superLpRatios );
    end
    if numel( im_bic ) > 1
      origBics = cell2mat( origBics );
      superBics = cell2mat( superBics );
    end
    resizedProtons = cell2mat( resizedProtons );

    nImgsPerRow = 5;

    f = figure;
    showImageCube( resizedProtons/max(resizedProtons(:))*255, 3, ...
      'nImgsPerRow', nImgsPerRow, 'border', 3, 'borderValue', 'max', 'range', 'nice' );
    pFrame=getframe(gcf);  proton4Fusion = double( pFrame.cdata ) / 255.;
    if showImages == false, close(f); end
    proton4Fusion = mean( proton4Fusion, 3 );
    imwrite( proton4Fusion, [ outDir, '/protons.jpg' ] );

    f = figure;
    showImageCube( origPyrs, 3, 'nImgsPerRow', nImgsPerRow, 'border', 3, 'borderValue', 'max' );
    titlenice( 'Original Pyruvate Images' );
    pFrame=getframe(gcf);  origPyrs = double( pFrame.cdata ) / 255.;
    if showImages ~= true, close(f); end
    origPyrs = mean( origPyrs, 3 );
    imwrite( origPyrs, [ outDir, '/pyruvatesOrig.jpg' ] );

    f = figure;
    showImageCube( superPyrs, 3, 'nImgsPerRow', nImgsPerRow, 'border', 3, 'borderValue', 'max' );
    titlenice( 'Pyruvate Images' );
    pFrame=getframe(gcf);  pyruvates = double( pFrame.cdata ) / 255.;
    if showImages ~= true, close(f); end
    pyruvates = mean( pyruvates, 3 );
    imwrite( pyruvates, [ outDir, '/pyruvatesSR.jpg' ] );
    
    f = figure;
    showImageCube( superPyrs, 3, 'nImgsPerRow', nImgsPerRow, 'border', 3, 'borderValue', 'max' );
    titlenice( 'Pyruvate Images' );  colormap( gca, 'hot' );
    pFrame=getframe(gcf);  pHot = double( pFrame.cdata ) / 255.;
    if showImages == false, close(f); end
    imwrite( pHot, [ outDir, '/pyruvatesSR_hot.jpg' ] );

    pFused = alphaBlend( pHot, proton4Fusion, 0.3 );
    if showImages == true
      figure; imshowscale( pFused );  titlenice( 'Pyruvate Fused - Alpha Blend' );
    end
    imwrite( pFused, [ outDir, '/fusedPyrAlpha.jpg' ] );

    pFused = clsFusion( pHot, proton4Fusion, gamma );
    if showImages == true
      figure; imshowscale( pFused );  titlenice( 'Pyruvate Fused - CLS' );
    end
    imwrite( pFused, [ outDir, '/fusedPyrCLS.jpg' ] );

    if numel( im_lac ) > 1
      
      f = figure;
      showImageCube( origLacs, 3, 'nImgsPerRow', nImgsPerRow, 'border', 3, 'borderValue', 'max' );
      titlenice( 'Original Lactate Images' );
      pFrame=getframe(gcf);  lactates = double( pFrame.cdata ) / 255.;
      if showImages ~= true, close(f); end
      imwrite( lactates, [ outDir, '/lactatesOrig.jpg' ] );

      f = figure;
      showImageCube( origLacs, 3, 'nImgsPerRow', nImgsPerRow, 'border', 3, 'borderValue', 'max' );
      titlenice( 'Lactate Images' );  colormap( gca, 'hot' );
      pFrame=getframe(gcf);  lacHot = double( pFrame.cdata ) / 255.;
      if showImages ~= true, close(f); end
      imwrite( lacHot, [ outDir, '/lactates_hot.jpg' ] );

      f = figure;
      showImageCube( superLacs, 3, 'nImgsPerRow', nImgsPerRow, 'border', 3, 'borderValue', 'max' );
      titlenice( 'Lactate Images' );
      pFrame=getframe(gcf);  lacHot = double( pFrame.cdata ) / 255.;
      if showImages ~= true, close(f); end
      imwrite( lacHot, [ outDir, '/lactatesSR.jpg' ] );

      f = figure;
      showImageCube( superLacs, 3, 'nImgsPerRow', nImgsPerRow, 'border', 3, 'borderValue', 'max' );
      titlenice( 'Lactate Images' );  colormap( gca, 'hot' );
      pFrame=getframe(gcf);  lacHot = double( pFrame.cdata ) / 255.;
      if showImages ~= true, close(f); end
      imwrite( lacHot, [ outDir, '/lactatesSR_hot.jpg' ] );
      
      lacFused = alphaBlend( lacHot, proton4Fusion, 0.3 );
      if showImages == true
        figure; imshowscale( lacFused );  titlenice( 'Lactate Fused - Alpha' );
      end
      imwrite( lacFused, [ outDir, '/fusedLacAlpha.jpg' ] );

      lacFused = clsFusion( lacHot, proton4Fusion, gamma );
      if showImages == true
        figure; imshowscale( lacFused );  titlenice( 'Lactate Fused - CLS' );
      end
      imwrite( lacFused, [ outDir, '/fusedLacCLS.jpg' ] );


      f = figure;
      showImageCube( origLpRatios, 3, 'nImgsPerRow', nImgsPerRow, 'border', 3, 'borderValue', 'max' );
      titlenice( 'Original L/P Images' );
      pFrame=getframe(gcf);  lpRatiosOrig = mean( double( pFrame.cdata ) / 255., 3 );
      if showImages ~= true, close(f); end
      imwrite( lpRatiosOrig, [ outDir, '/lpRatiosOrig.jpg' ] );

      f = figure;
      showImageCube( superLpRatios, 3, 'nImgsPerRow', nImgsPerRow, 'border', 3, 'borderValue', 'max' );
      titlenice( 'L/P Images' );
      pFrame=getframe(gcf);  lpRatiosSR = mean( double( pFrame.cdata ) / 255., 3 );
      if showImages ~= true, close(f); end
      imwrite( lpRatiosSR, [ outDir, '/lpRatiosSR.jpg' ] );

      f = figure;
      showImageCube( superLpRatios, 3, 'nImgsPerRow', nImgsPerRow, 'border', 3, 'borderValue', 'max' );
      titlenice( 'L/P Images' );  colormap( gca, 'hot' );
      lpRatiosFrame=getframe(gcf);  lpRatiosColor = double( lpRatiosFrame.cdata ) / 255.;
      if showImages == false, close(f); end
      imwrite( lpRatiosColor, [ outDir, '/lpRatiosColor.jpg' ] );

      lpRatiosAlphaFused = alphaBlend( lpRatiosColor, proton4Fusion, 0.3 );
      if showImages == true
        figure;
        imshowscale( lpRatiosAlphaFused );  titlenice( 'L/P Fused - Alpha' );
      end
      imwrite( lpRatiosAlphaFused, [ outDir, '/fusedLpRatioAlpha.jpg' ] );

      lpRatiosFused = clsFusion( lpRatiosColor, proton4Fusion, gamma );
      if showImages == true
        figure; imshowscale( lpRatiosFused );  titlenice( 'L/P Fused - CLS' );
      end
      imwrite( lpRatiosFused, [ outDir, '/fusedLpRatioCLS.jpg' ] );
    end

    if numel( im_bic ) > 1
      f = figure;
      showImageCube( origBics, 3, 'nImgsPerRow', nImgsPerRow, 'border', 3, 'borderValue', 'max' );
      titlenice( 'Original Bicarbonate Images' );
      pFrame=getframe(gcf);  bicarbsOrig = mean( double( pFrame.cdata ) / 255., 3 );
      if showImages ~= true, close(f); end
      imwrite( bicarbsOrig, [ outDir, '/bicarbsOrig.jpg' ] );

      f = figure;
      showImageCube( superBics, 3, 'nImgsPerRow', nImgsPerRow, 'border', 3, 'borderValue', 'max' );
      titlenice( 'Bicarbonate Images SuperResolved' );
      pFrame=getframe(gcf);  bicarbsSR = mean( double( pFrame.cdata ) / 255., 3 );
      if showImages ~= true, close(f); end
      imwrite( bicarbsSR, [ outDir, '/bicarbsSR.jpg' ] );

      f = figure;
      showImageCube( superBics, 3, 'nImgsPerRow', nImgsPerRow, 'border', 3, 'borderValue', 'max' );
      titlenice( 'Bicarbonate Images' );  colormap( gca, 'hot' );
      bFrame=getframe(gcf);  bColor = double( bFrame.cdata ) / 255.;
      if showImages == false, close(f); end
      imwrite( bColor, [ outDir, '/bicarbsSR_hot.jpg' ] );

      bicFused = alphaBlend( bColor, proton4Fusion, 0.3 );
      if showImages == true
        figure; imshowscale( bicFused );  titlenice( 'Bicarbonate Fused - Alpha' );
      end
      imwrite( bicFused, [ outDir, '/fusedBicAlpha.jpg' ] );

      bicFused = clsFusion( bColor, proton4Fusion, gamma );
      if showImages == true
        figure; imshowscale( bicFused );  titlenice( 'Bicarbonate Fused - CLS' );
      end
      imwrite( bicFused, [ outDir, '/fusedBicCLS.jpg' ] );
    end

  end

end

