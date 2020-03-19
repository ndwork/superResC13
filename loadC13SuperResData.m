
function [im_proton, im_pyr, im_lac, im_bic, im_lpRatio, lambda, subRegion, ds] = ...
  loadC13SuperResData( datacase )

  im_proton=0;  im_pyr=0;  im_lac=0;  im_bic=0;  im_lpRatio=0;

  switch datacase
    case 1
      %lambda = 1d4;
      lambda = 1d3;
      load( '../data/oneSliceHeart/data.mat' );
      subRegion = [];  % xL yL xR yR
      ds = 5;

    case 2
      %lambda = 1d4;
      lambda = 1d3;
      load( '../data/manySlicesHeart/data.mat' );
      ds = 5;
      subRegion = [ 30 65 150 165 ];  % xL yL xR yR

    case 3
      %lambda = 1d4;
      lambda = 1d3;
      areas=0;  % stored within mat file
      load( '../data/prostate/pc9154_data.mat' );
      subRegion = [];  % xL yL xR yR

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

      pngImg = cropData( pngImg, ds*[ size(rotAreas,1) size(rotAreas,2) ] );
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
      
    case 4
      %lambda = 1d4;
      lambda = 1d3;
      areas=0;  %#ok<*NASGU> % stored within mat file
      load( '../data/prostate/pc9154_data_cc.mat' );
      subRegion = [];  % xL yL xR yR

      pngImg = imread( '../data/prostate/pc9154.png' );
      pngImg = double( pngImg(:,:,1) ) / 255.;

      areas = zeros( 12, 12, 16, 3, 21 );
      areas(:,:,:,1,:) = squeeze( Data_4D_coilcorrect.Pyruvate );
      areas(:,:,:,2,:) = squeeze( Data_4D_coilcorrect.Lactate );

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
      im_pyr(1,:,:,:,:) = 0;  im_pyr(end,:,:,:,:) = 0;

      im_lac = abs( squeeze( rotAreas( :, :, 9, 2, 7:end ) ) );
      im_lac(1,:,:,:,:) = 0;  im_lac(end,:,:,:,:) = 0;

      pngImg = cropData( pngImg, ds*[ size(rotAreas,1) size(rotAreas,2) ] );
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

    case 5
      caseDataDir = '/Users/nicholasdwork/.DellEMS/Vols/disk2s1/c13_heart/';
      %lambda = 1d4;
      lambda = 1d3;
      pyr_images = [];  lac_images = [];  bicarb_images = [];
      load([ caseDataDir, 'c13_data.mat' ]);
      localizers = [];
      load([ caseDataDir, 'proton_localizers.mat' ]);
      ds = 3;
      %subRegion = [ 108 73 183 150 ];
      subRegion = [ 118 83 173 140 ];

      timeIndx = 6;
      im_pyr_small = squeeze( abs( pyr_images( :, :, :, timeIndx ) ) );
      im_lac_small = squeeze( abs( lac_images( :, :, :, timeIndx ) ) );
      im_bic_small = squeeze( abs( bicarb_images( :, :, :, timeIndx ) ) );
      im_proton = squeeze( abs( localizers ) );

      nSlices = size( im_pyr_small, 3 );
      im_pyr = zeros( size( im_proton ) );
      im_lac = zeros( size( im_proton ) );
      im_bic = zeros( size( im_proton ) );
      for sliceIndx = 1 : nSlices
        im_pyr(:,:,sliceIndx) = imresize( im_pyr_small(:,:,sliceIndx), size(im_proton(:,:,1)) );
        im_lac(:,:,sliceIndx) = imresize( im_lac_small(:,:,sliceIndx), size(im_proton(:,:,1)) );
        im_bic(:,:,sliceIndx) = imresize( im_bic_small(:,:,sliceIndx), size(im_proton(:,:,1)) );
      end

    case 6
      caseDataDir = '/Users/nicholasdwork/.DellEMS/Vols/disk2s1/c13_heart/';
      %lambda = 1d4;
      lambda = 1d3;
      pyr_images = [];  lac_images = [];  bicarb_images = [];
      load([ caseDataDir, 'c13_data.mat' ]);
      localizers = [];
      load([ caseDataDir, 'proton_localizers.mat' ]);
      ds = 3;
      subRegion = [ 63 79 197 163 ];

      sliceIndx = 4;
      timeIndx = 7;
      im_pyr_small = squeeze( abs( pyr_images( :, :, sliceIndx, timeIndx ) ) );
      im_lac_small = squeeze( abs( lac_images( :, :, sliceIndx, timeIndx ) ) );
      im_bic_small = squeeze( abs( bicarb_images( :, :, sliceIndx, timeIndx ) ) );
      im_proton = squeeze( abs( localizers(:,:,sliceIndx) ) );

      im_pyr = imresize( im_pyr_small, size( im_proton ) );
      im_lac = imresize( im_lac_small, size( im_proton ) );
      im_bic = imresize( im_bic_small, size( im_proton ) );
  end

end
