
function [M,C,sigma,lambda] = loadPanSharpenData( datacase )

    switch datacase
      case 0
        M = phantom();
        R = imresize( phantom(), ceil(size(M)*0.75), 'bilinear' );
        G = imresize( phantom(), ceil(size(M)*0.75), 'bilinear' );
        B = imresize( phantom(), ceil(size(M)*0.75), 'bilinear' );
        C = cat( 3, R, G, B );
        sigma = 0.5;
        lambda = 1d3;

      case 1
        color = imread( '../data/pics/rocks.png' );
        color = double( color ) / 255;
        M = 0.2989 * color(:,:,1) + 0.5870 * color(:,:,2) + 0.1140 * color(:,:,3);
        sigma = 2.0;
        C = smoothImg( color, 'gaussian', sigma );
        C = cropData( C, [size(M)-9 3] );
        M = cropData( M, size(M)-9 );
        C = C( 2:4:end, 2:4:end, : );
        lambda = 1d3;

      case 2
        load( '../data/satellite/satellite.mat' );  % quickbird
        C = double(color)/255;
        M = double(monochrome)/255;
        sigma = res/2.0;
        lambda = 1d3;

    end
  
end