
function run_panSharpen
  close all;  clear;  rng(1);

  showScale = 5;

  for datacase = 2

    [M,C,sigma] = loadPanSharpenData( datacase );

    tic
    out = panSharpen( M, C, 'sigma', sigma );
    timeTaken = toc;

    disp([ 'Time taken: ', num2str( timeTaken ), ' (s)' ]);
    figure; imshowscale( M, showScale );  titlenice( 'Monochrome Image' );
    figure; imshowscale( imColorResize( C(:,:,1:3), [size(M) 3], 'bilinear' ), showScale );
    titlenice( 'Original Color' );
    figure; imshowscale( out, showScale );  titlenice( 'Sharpened Color' );
  end

end

