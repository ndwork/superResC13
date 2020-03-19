
function run_dicom( datacases )
  % run_dicom( [ datacase ] )
  close all;  rng(1);

  if nargin < 1
    clear;
    datacases=4;
  end

  showImages = false;
  scaleAmount = 0.4;

  for datacase = datacases
    showScale = 5;

    [t1Dir, t2Dir, pDir, lDir, bDir, pds, lds, bds, imgSliceIndxs] = ...
      loadDicomDatacase( datacase );

    for imgSliceIndx = imgSliceIndxs

      outDir = [ './output/output_dicom_', num2str(datacase), '/slice_', ...
        indx2str(imgSliceIndx,max(imgSliceIndxs)), '/' ];
      if ~exist( outDir, 'dir' ), mkdir( outDir ); end
      
      if numel( t2Dir ) > 0
        [t1,t2] = alignDicoms( t1Dir, t2Dir );
        [~,pyr] = alignDicoms( t1Dir, pDir );
      else
        t2 = [];  t2Img = [];
        [t1,pyr] = alignDicoms( t1Dir, pDir );
      end
      [~,lac] = alignDicoms( t1Dir, lDir );
      [~,bic] = alignDicoms( t1Dir, bDir );

      t1Img = t1( :, :, imgSliceIndx );
      if numel( t2 ) > 0, t2Img = t2( :, :, imgSliceIndx ); end
      pImg = pyr( :, :, imgSliceIndx );
      lImg = lac( :, :, imgSliceIndx );
      bImg = bic( :, :, imgSliceIndx );

      if showImages == true
        figure; imshowscale( t1Img, showScale );  titlenice( 't1' );
        if numel( t2 ) > 0, figure; imshowscale( t2Img, showScale );  titlenice( 't2' ); end
      end


      noisyP = pImg(end-50:end,end-50:end);  pImg = max( pImg - median( noisyP(:) ), 0 );
      noisyL = lImg(end-50:end,end-50:end);  lImg = max( lImg - median( noisyL(:) ), 0 );
      noisyB = bImg(end-50:end,end-50:end);  bImg = max( bImg - median( noisyB(:) ), 0 );

      pImg = scaleImg( pImg, [0 5000] ) * 5000;
      lImg = scaleImg( lImg, [0 500] ) * 500;

      maxT1 = max( t1Img(:) );    t1Img = t1Img ./ maxT1;
      if numel( t2 ) > 0
        maxT2 = max( t2Img(:) );  t2Img = t2Img ./ maxT2;
      end
      maxPyr = max( pImg(:) );    pImg = pImg ./ maxPyr;
      maxLac = max( lImg(:) );    lImg = lImg ./ maxLac;
      maxBic = max( bImg(:) );    bImg = bImg ./ maxBic;

      pImg = pImg( ceil(pds/2) : pds : end, ceil(pds/2) : pds : end );
      lImg = lImg( ceil(lds/2) : lds : end, ceil(lds/2) : lds : end );
      bImg = bImg( ceil(bds/2) : bds : end, ceil(bds/2) : bds : end );

      tmp = imresize( pImg, size(t1Img), 'nearest' );
      if showImages == true
        figure; imshowscale( tmp, showScale );
        titlenice( 'pyr' );  %colormap( gca, 'hot' );
      end
      imwrite( tmp / max( tmp(:) ), ...
        [ outDir, '/origPyr_', num2str(datacase,'%3.3i'), '.jpg'] );
      imwrite( scaleImg( tmp / max( tmp(:) ), [0 scaleAmount] ), ...
        [ outDir, '/origPyr_scaled_', num2str(datacase,'%3.3i'), '.jpg'] );

      tmp = imresize( lImg, size(t1Img), 'nearest' );
      if showImages == true
        figure; imshowscale( tmp, showScale );
        titlenice( 'lac' );  %colormap( gca, 'hot' );
      end
      imwrite( tmp / max( tmp(:) ), ...
        [ outDir, '/origLac_', num2str(datacase,'%3.3i'), '.jpg'] );
      imwrite( scaleImg( tmp / max( tmp(:) ), [0 scaleAmount] ), ...
        [ outDir, '/origLac_scaled_', num2str(datacase,'%3.3i'), '.jpg'] );

      tmp = imresize( bImg, size(t1Img), 'nearest' );
      if showImages == true
        figure; imshowscale( tmp, showScale );
        titlenice( 'bic' );  %colormap( gca, 'hot' );
      end
      imwrite( tmp / max( tmp(:) ), ...
        [ outDir, '/origBic_', num2str(datacase,'%3.3i'), '.jpg'] );
      imwrite( scaleImg( tmp / max( tmp(:) ), [0 scaleAmount] ), ...
        [ outDir, '/origBic_scaled_', num2str(datacase,'%3.3i'), '.jpg'] );

      outputs = cell(1,7);
      parfor i = 1 : numel( outputs )

        switch i
          case 1
            lambda = 1d2;
            superPyrT1 = superResC13( t1Img, pImg, 0.3 * pds, 'lambda', lambda );
            superPyrT1 = superPyrT1 * maxPyr;
            outputs{i} = superPyrT1;

          case 2
            lambda = 1d2;
            superPyrT1 = superResC13( t1Img, lImg, 0.3 * lds, 'lambda', lambda );
            superPyrT1 = superPyrT1 * maxPyr;
            outputs{i} = superPyrT1;

          case 3
            if numel( t2 ) > 0
              lambda = 1d2;
              superPyrT2 = superResC13( t2Img, pImg, 0.3 * pds, 'lambda', lambda );
              superPyrT2 = superPyrT2 * maxPyr;
              outputs{i} = superPyrT2;
            end

          case 4
            if numel( t2 ) > 0
              lambda = 1d2;
              superLacT2 = superResC13( t2Img, lImg, 0.3 * lds, 'lambda', lambda );
              superLacT2 = superLacT2 * maxLac;
              outputs{i} = superLacT2;
            end

          case 5
            if numel( t2 ) > 0
              lambda = 1d2;
              superBicT2 = superResC13( t2Img, bImg, 0.3 * bds, 'lambda', lambda );
              superBicT2 = superBicT2 * maxBic;
              outputs{i} = superBicT2;
            end

          case 6
            lambda = 1d2;
            lpds = lds / pds;
            superPL = superResC13( pImg, lImg, 0.3 * lpds, 'lambda', lambda );
            superPL = superPL * maxLac;
            outputs{i} = superPL;

          case 7
            lambda = 1d2;
            bpds = bds / pds;
            superPB = superResC13( pImg, bImg, 0.3 * bpds, 'lambda', lambda );
            superPB = superPB * maxBic;
            outputs{i} = superPB;
        end
      end

      if numel( outputs{1} ) > 0
        superPyrT1 = outputs{1};
        if showImages == true
          figure; imshowscale( superPyrT1, showScale );  titlenice( 'superPyr/T1' );
        end
        imwrite( superPyrT1 / max( superPyrT1(:) ), ...
          [ outDir, '/superPyrT1_', num2str(datacase,'%3.3i'), '.jpg'] );
        imwrite( scaleImg( superPyrT1 / max( superPyrT1(:) ), [0 scaleAmount] ), ...
          [ outDir, '/superPyrT1_scaled_', num2str(datacase,'%3.3i'), '.jpg'] );
      end

      if numel( outputs{2} ) > 0
        superLacT1 = outputs{2};
        if showImages == true
          figure; imshowscale( superLacT1, showScale );  titlenice( 'superPyr/T1' );
        end
        imwrite( superLacT1 / max( superLacT1(:) ), ...
          [ outDir, '/superPyrT1_', num2str(datacase,'%3.3i'), '.jpg'] );
        imwrite( scaleImg( superLacT1 / max( superLacT1(:) ), [0 scaleAmount] ), ...
          [ outDir, '/superPyrT1_scaled_', num2str(datacase,'%3.3i'), '.jpg'] );
      end

      if numel( outputs{3} ) > 0
        superPyrT2 = outputs{3};
        if showImages == true
          figure;  imshowscale( superPyrT2, showScale );  titlenice( 'superPyr/T2' );
        end
        imwrite( superPyrT2 / max( superPyrT2(:) ), ...
          [ outDir, '/superPyrT2_', num2str(datacase,'%3.3i'), '.jpg'] );
        imwrite( scaleImg( superPyrT2 / max( superPyrT2(:) ), [0 scaleAmount] ), ...
          [ outDir, '/superPyrT2_scaled_', num2str(datacase,'%3.3i'), '.jpg'] );
      end

      if numel( outputs{4} ) > 0
        superLacT2 = outputs{4};
        if showImages == true
          figure;  imshowscale( superLacT2, showScale );  titlenice( 'superLac/T2' );
        end
        imwrite( superLacT2 / max( superLacT2(:) ), ...
          [ outDir, '/superLacT2_', num2str(datacase,'%3.3i'), '.jpg'] );
        imwrite( scaleImg( superLacT2 / max( superLacT2(:) ), [0 scaleAmount] ), ...
          [ outDir, '/superLacT2_scaled_', num2str(datacase,'%3.3i'), '.jpg'] );
      end

      if numel( outputs{5} ) > 0
        superBicT2 = outputs{5};
        if showImages == true
          figure;  tmp = imresize( superBicT2, size(t1Img), 'nearest' );
          imshowscale( tmp, showScale );  titlenice( 'superBic/T2' );
        end
        imwrite( superBicT2 / max( superBicT2(:) ), ...
          [ outDir, '/superBicT2_', num2str(datacase,'%3.3i'), '.jpg'] );
        imwrite( scaleImg( superBicT2 / max( superBicT2(:) ), [0 scaleAmount] ), ...
          [ outDir, '/superBicT2_scaled_', num2str(datacase,'%3.3i'), '.jpg'] );
      end

      if numel( outputs{6} ) > 0
        superLacPyr = outputs{6};
        if showImages == true
          figure;  tmp = imresize( superLacPyr, size(t1Img), 'nearest');
          imshowscale( tmp, showScale );  titlenice( 'superLac/Pyr' );
        end
        imwrite( imresize(superLacPyr, size(t1Img), 'nearest' ) / max( superLacPyr(:) ), ...
          [ outDir, '/superPL_', num2str(datacase,'%3.3i'), '.jpg'] );
        imwrite( scaleImg( imresize(superLacPyr, size(t1Img), 'nearest') / max( superLacPyr(:) ), ...
          [0 scaleAmount] ), [ outDir, '/superPL_scaled_', num2str(datacase,'%3.3i'), '.jpg'] );
      end

      if numel( outputs{7} ) > 0
        superBicPyr = outputs{7};
        if showImages == true
          figure;  tmp = imresize( superBicPyr, size(t1Img), 'nearest');
          imshowscale( tmp, showScale );  titlenice( 'superBic/Pyr' );
        end
        imwrite( imresize( superBicPyr, size(t1Img), 'nearest') / max( superBicPyr(:) ), ...
          [ outDir, '/superPB_', num2str(datacase,'%3.3i'), '.jpg'] );
        imwrite( scaleImg( imresize( superBicPyr, size(t1Img), 'nearest') / max( superBicPyr(:) ), ...
          [0 scaleAmount] ), [ outDir, '/superPB_scaled_', num2str(datacase,'%3.3i'), '.jpg'] );
      end

    end
  end

end

