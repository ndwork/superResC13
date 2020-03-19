
function run_makePaperImages

  disp( 'run_severalLambdas' );
  run_severalLambdas([ 3 1 ]);
  close all;

  disp( 'run_c13superRes' );
  parfor datacase = 1 : 4  
    run_c13superRes( datacase )
    close all;
  end

  parfor imgType = 1 : 4

    if imgType <= 2
      disp( 'run_validateSuperRes' );
      %datacases = 1:2;
      datacase = imgType;
      run_validateSuperRes( datacase );
      close all;
    end

    if imgType == 3
      disp( 'run_superResFalseColor' );
      run_superResFalseColor
      close all;
    end

    if imgType == 4
      disp( 'run_dicom' );
      run_dicom
      close all;
    end
  end

end

