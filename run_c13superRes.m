
function run_c13superRes
  close all;  clear;  rng(1);
  
  showImages = false;

  for datacase = 3:-1:1
    if ~exist( 'showImages', 'var' ) || showImages == false
      close all;
    end

    outDir = [ './output_', num2str(datacase) ];
    mkdir( outDir );

    im_proton=0;  im_pyr=0;  im_lac=0;  im_bic=0;  im_lpRatio=0;
    switch datacase
      case 1
        lambda = 1d4;
        load( '../data/oneSliceHeart/data.mat' );
        ds = 5;

      case 2
        lambda = 1d4;
        load( '../data/manySlicesHeart/data.mat' );
        ds = 5;

      case 3
        %lambda = 1d1;
        lambda = 1d4;
        areas=0;  % stored within mat file
        load( '../data/prostate/pc9154_data.mat' );

        pngImg = imread( '../data/prostate/pc9154.png' );
        pngImg = double( pngImg(:,:,1) ) / 255.;

        ds = size( pngImg ) ./ [ size(areas,1) size(areas,2) ];
        ds = min( floor( ds ) );

        sAs = size( areas );
        rotAreas = zeros( [ sAs(2), sAs(1), sAs(3:end) ] );
        for k = 1:sAs(5)
          for j = 1:sAs(4)
            for i = 1:sAs(3)
              rotAreas(:,:,i,j,k) = circshift( rot90( areas(:,:,i,j,k), -1 ), [-1 0] );
            end
          end
        end
        im_pyr = abs( squeeze( rotAreas( :, :, 9, 1, 7:end ) ) );
        im_lac = abs( squeeze( rotAreas( :, :, 9, 2, 7:end ) ) );

        pngImg = cropData( pngImg, ds*[ size(areas,1) size(areas,2) ] );
        im_proton = repmat( pngImg, [1 1 size( im_pyr, 3)] );
        clear rotAreas areas;

        sP = size( im_proton );
        resizedPyr = zeros( sP );
        resizedLac = zeros( sP );
        for sliceIndx = 1 : sP(3)
          resizedPyr(:,:,sliceIndx) = imresize( im_pyr(:,:,sliceIndx), [sP(1), sP(2)], 'nearest' );
          resizedLac(:,:,sliceIndx) = imresize( im_lac(:,:,sliceIndx), [sP(1), sP(2)], 'nearest' );
        end
        im_pyr = resizedPyr;  clear resizedPyr;
        im_lac = resizedLac;  clear resizedLac;
    end
    

    if numel( im_lac ) > 1
      im_lpRatio = im_pyr ./ im_lac;
      mask = im_pyr > 0.1 * max( im_pyr(:) );
      im_lpRatio( mask==0 ) = 0;
      maxLpRatio = max( im_lpRatio(:) );
      im_lpRatio = im_lpRatio / maxLpRatio;

      maxLac = max( im_lac(:) );
      im_lac = im_lac / maxLac;      
    end

    maxProton = max( im_proton(:) );
    im_proton = im_proton / maxProton;

    maxPyr = max( im_pyr(:) );
    im_pyr = im_pyr / maxPyr;


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
%for sliceIndx = 1 : nSlices
      disp([ 'Working on slice ', num2str( sliceIndx ) ]);
      thisProton = smoothImg( im_proton(:,:,sliceIndx), 5 );
      imwrite( thisProton, [ outDir, '/thisProton_', num2str(sliceIndx,'%3.3i'), '.jpg'] );

      thisPyr = im_pyr( ceil(ds/2) : ds : end, ceil(ds/2) : ds : end, sliceIndx );    %#ok<PFBNS>
      thisPyrResized = imresize( thisPyr, size(thisProton), 'nearest' );
      origPyrs{ sliceIndx } = thisPyrResized;
      imwrite( thisPyrResized, [ outDir, '/thisPyr_', num2str(sliceIndx,'%3.3i'), '.jpg'] );
      
      thisSuperPyr = superResC13( thisProton, thisPyr, 0.5*ds, 'lambda', lambda );
      imwrite( thisSuperPyr, [ outDir, '/thisSuperPyr_', num2str(sliceIndx,'%3.3i'), '.jpg'] );
      superPyrs{ sliceIndx } = thisSuperPyr;

      if numel( im_lac ) > 1
        thisLac = im_lac( ceil(ds/2) : ds : end, ceil(ds/2) : ds : end, sliceIndx );
        thisLacResized = imresize( thisLac, size( thisProton ), 'nearest' );
        origLacs{ sliceIndx } = thisLacResized;
        imwrite( thisLacResized, [ outDir, '/thisLac_', num2str(sliceIndx,'%3.3i'), '.jpg'] );

        thisSuperLac = superResC13( thisProton, thisLac, 0.5*ds, 'lambda', lambda );
        imwrite( thisSuperLac, [ outDir, '/thisSuperLac_', num2str(sliceIndx,'%3.3i'), '.jpg'] );
        superLacs{ sliceIndx } = thisSuperLac;

        thisRatio = im_lpRatio( ceil(ds/2) : ds : end, ceil(ds/2) : ds : end, sliceIndx );    %#ok<PFBNS>
        thisRatioResized = imresize( thisRatio, size( thisProton ), 'nearest' );
        origLpRatios{ sliceIndx } = thisRatioResized;
        imwrite( thisRatioResized, [ outDir, '/thisLpRatio_', num2str(sliceIndx,'%3.3i'), '.jpg'] );

        thisSuperLpRatio = superResC13( thisProton, thisRatio, 0.5*ds, 'lambda', lambda );
        imwrite( thisSuperLpRatio, [ outDir, '/thisSuperRatio_', num2str(sliceIndx,'%3.3i'), '.jpg'] );
        superLpRatios{ sliceIndx } = thisSuperLpRatio;
      end

      if numel( im_bic ) > 1
        thisBic = im_bic( ceil(ds/2) : ds : end, ceil(ds/2) : ds : end, sliceIndx );
        thisBicResized = imresize( thisBic, size(thisProton), 'nearest' );
        origBics{ sliceIndx } = thisBicResized;
        imwrite( thisBicResized, [ outDir, '/thisBic_', num2str(sliceIndx,'%3.3i'), '.jpg'] );

        thisSuperBic = superResC13( thisProton, thisBic, 0.5*ds, 'lambda', lambda );
        imwrite( thisSuperBic, [ outDir, '/thisSuperBic_', num2str(sliceIndx,'%3.3i'), '.jpg'] );
        superBics{ sliceIndx } = thisSuperBic;
      end

      thisResizedProton = thisProton( 1 : size(thisSuperPyr,1), 1 : size(thisSuperPyr,2) );
      resizedProtons{ sliceIndx } = thisResizedProton;
      imwrite( thisResizedProton, ...
        [ outDir, '/thisResizedProton_', num2str(sliceIndx,'%3.3i'), '.jpg'] );
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

    figure; showImageCube( resizedProtons/max(resizedProtons(:))*255, 3, ...
      'nImgsPerRow', nImgsPerRow, 'border', 3, 'borderValue', 'max', 'range', 'nice' );
    pFrame=getframe(gcf);  proton4Fusion = double( pFrame.cdata ) / 255.;
    proton4Fusion = mean( proton4Fusion, 3 );

    figure;
    showImageCube( origPyrs, 3, 'nImgsPerRow', nImgsPerRow, 'border', 3, 'borderValue', 'max' );
    titlenice( 'Original Pyruvate Images' );
    figure;
    showImageCube( superPyrs, 3, 'nImgsPerRow', nImgsPerRow, 'border', 3, 'borderValue', 'max' );
    titlenice( 'Pyruvate Images' );

    figure; showImageCube( superPyrs, 3, 'nImgsPerRow', nImgsPerRow, 'border', 3, 'borderValue', 'max' );
    titlenice( 'Pyruvate Images' );  colormap( gca, 'hot' );
    pFrame=getframe(gcf);  pColor = double( pFrame.cdata ) / 255.;
    pFused = alphaBlend( pColor, 2*proton4Fusion, 0.3 );
    figure; imshowscale( pFused );  titlenice( 'Pyruvate Fused - Alpha Blend' );
    imwrite( pFused, [ outDir, '/fusedPyrAlpha.jpg' ] );

    pFused = clsFusion( pColor, 4*proton4Fusion );
    figure; imshowscale( pFused );  titlenice( 'Pyruvate Fused - CLS' );
    imwrite( pFused, [ outDir, '/fusedPyrCLS.jpg' ] );

    if numel( im_lac ) > 1
      figure;
      showImageCube( origLacs, 3, 'nImgsPerRow', nImgsPerRow, 'border', 3, 'borderValue', 'max' );
      titlenice( 'Original Lactate Images' );
      figure;
      showImageCube( superLacs, 3, 'nImgsPerRow', nImgsPerRow, 'border', 3, 'borderValue', 'max' );
      titlenice( 'Lactate Images' );

      figure;
      showImageCube( superLacs, 3, 'nImgsPerRow', nImgsPerRow, 'border', 3, 'borderValue', 'max' );
      titlenice( 'Lactate Images' );  colormap( gca, 'hot' );
      lacFrame=getframe(gcf);  lacColor = double( lacFrame.cdata ) / 255.;
      lacFused = alphaBlend( lacColor, 2*proton4Fusion, 0.3 );
      figure; imshowscale( lacFused );  titlenice( 'Lactate Fused - Alpha' );
      imwrite( lacFused, [ outDir, '/fusedLacAlpha.jpg' ] );

      lacFused = clsFusion( lacColor, 4*proton4Fusion );
      figure; imshowscale( lacFused );  titlenice( 'Lactate Fused - CLS' );
      imwrite( lacFused, [ outDir, '/fusedLacCLS.jpg' ] );

      figure;
      showImageCube( origLpRatios, 3, 'nImgsPerRow', nImgsPerRow, 'border', 3, 'borderValue', 'max' );
      titlenice( 'Original L/P Images' );
      figure;
      showImageCube( superLpRatios, 3, 'nImgsPerRow', nImgsPerRow, 'border', 3, 'borderValue', 'max' );
      titlenice( 'L/P Images' );
      
      figure;
      showImageCube( superLpRatios, 3, 'nImgsPerRow', nImgsPerRow, 'border', 3, 'borderValue', 'max' );
      titlenice( 'L/P Images' );  colormap( gca, 'hot' );
      lpRatiosFrame=getframe(gcf);  lpRatiosColor = double( lpRatiosFrame.cdata ) / 255.;
      lpRatiosFused = alphaBlend( lpRatiosColor, 2*proton4Fusion, 0.3 );
      figure; imshowscale( lpRatiosFused );  titlenice( 'L/P Fused - Alpha' );
      imwrite( lpRatiosFused, [ outDir, '/fusedLpRatioAlpha.jpg' ] );

      lpRatiosFused = clsFusion( 0.5*lpRatiosColor, 4*proton4Fusion );
      figure; imshowscale( lpRatiosFused );  titlenice( 'L/P Fused - CLS' );
      imwrite( lpRatiosFused, [ outDir, '/fusedLpRatioCLS.jpg' ] );

    end

    if numel( im_bic ) > 1
      figure; showImageCube( origBics, 3, 'nImgsPerRow', nImgsPerRow, 'border', 3, 'borderValue', 'max' );
      titlenice( 'Original Bicarbonate Images' );
      figure; showImageCube( superBics, 3, 'nImgsPerRow', nImgsPerRow, 'border', 3, 'borderValue', 'max' );
      titlenice( 'Bicarbonate Images' );

      figure; showImageCube( superBics, 3, 'nImgsPerRow', nImgsPerRow, 'border', 3, 'borderValue', 'max' );
      titlenice( 'Bicarbonate Images' );  colormap( gca, 'hot' );
      bFrame=getframe(gcf);  bColor = double( bFrame.cdata ) / 255.;
      bicFused = alphaBlend( bColor, 2*proton4Fusion, 0.3 );
      figure; imshowscale( bicFused );  titlenice( 'Bicarbonate Fused - Alpha' );
      imwrite( bicFused, [ outDir, '/fusedBicAlpha.jpg' ] );

      bicFused = clsFusion( bColor, 4*proton4Fusion );
      figure; imshowscale( bicFused );  titlenice( 'Bicarbonate Fused - CLS' );
      imwrite( bicFused, [ outDir, '/fusedBicCLS.jpg' ] );
    end

  end

end

