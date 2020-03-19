
function run_severalLambdas( datacases )
  showImages = false;

  if nargin < 1
    close all; rng(1); clear;
    datacases = 5;
    showImages = true;
  end


  for datacase = datacases

    blurFraction = 0.5;

    if ~exist( 'showImages', 'var' ) || showImages == false, close all; end

    outDir = [ './output/severalLambdas_', num2str(datacase) ];
    mkdir( outDir );

    im_proton=0;  im_pyr=0;  im_lac=0;  im_bic=0;  im_lpRatio=0;   %#ok<NASGU>
    [im_proton, im_pyr, im_lac, im_bic, im_lpRatio, lambda, subRegion, ds] = ...
      loadC13SuperResData( datacase );   %#ok<ASGLU>

    maxProton = max( im_proton(:) );
    im_proton = im_proton / maxProton;

    maxPyr = max( im_pyr(:) );
    im_pyr = im_pyr / maxPyr;

    lambdaPowers = -2 : 1 : 2;
    nPowers = numel( lambdaPowers );
    superPyrs = cell( 1, 1, nPowers );

    nSlices = size( im_proton, 3 );
    for slice = 1 : nSlices
      thisProton = smoothImg( im_proton(:,:,slice), 5 );
      thisPyr = im_pyr( ceil(ds/2) : ds : end, ceil(ds/2) : ds : end, slice );

      if numel( subRegion ) > 0
        xL = subRegion(1);  xR = subRegion(3);
        yL = subRegion(2);  yR = subRegion(4);
      else
        sProton = size( thisProton );
        xL = 1;  xR = sProton(2);
        yL = 1;  yR = sProton(1);
      end

      parfor lambdaPowerIndx = 1 : nPowers
        disp([ 'Working on lambdaPower ', num2str( lambdaPowerIndx ), ' of ', num2str( nPowers ) ]);
        lambdaPower = lambdaPowers( lambdaPowerIndx );
        thisLambda = lambda * 10^lambdaPower;

        thisSuperPyr = superResC13( thisProton, thisPyr, blurFraction*ds, 'lambda', thisLambda );
        superPyrs{ lambdaPowerIndx } = thisSuperPyr( yL : yR, xL : xR );
      end
      
      superPyrsMat = cell2mat( superPyrs );
      figure; showImageCube( superPyrsMat, 3, 'border', 2, 'borderValue', 'max', 'nImgsPerRow', nPowers );
      if nSlices == 1
        saveas( gcf, [ outDir, '/superPyrs.jpg' ] );
      else
        saveas( gcf, [ outDir, '/superPyrs_', indx2str(slice,nSlices), '.jpg' ] );
      end
      close all;
    end

  end

end

