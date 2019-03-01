
function run_c13superRes
  close all;  clear;  rng(1);

  lambda = 2d2;

  im_pyr = 0;
  load( '../data/data.mat' );

  im_pyr = im_pyr( 4 : 9 : end, 4 : 9 : end );
  im_proton = smoothImg( im_proton, 5 );

  im_proton = im_proton / max( im_proton(:) );
  im_pyr = im_pyr / max( im_pyr(:) );

  superPyr = superResC13( im_proton, im_pyr, 'lambda', lambda );

  %figure; imshowscale( superPyr, 3 );

end

