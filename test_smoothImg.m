
function test_smoothImg
  close all;  clear;  rng(1);

  N = 7;

  ph = phantom();
  %ph = ones(5);  ph(3,3) = 1;

  sph = smoothImg( ph, N );

  sphT = smoothImg( sph, N, 'op', 'transp' );

  figure; imshowscale( [ ph sph sphT ] );  titlenice( 'Smoothed' );

  
  f = @(x) smoothImg( x, N );
  fAdj = @(y) smoothImg( y, N, 'op', 'transp' );

  [checkPassed,err] = checkAdjoint( ph, f, fAdj );
  if checkPassed
    disp( 'Adjoint test passed' );
  else
    disp( [ 'Error: Adjoint test failed with error: ', num2str(err) ] );
  end

end
