
function [t1Dir, t2Dir, pDir, lDir, bDir, pds, lds, bds, imgSliceIndxs] = ...
  loadDicomDatacase( datacase )

  dataDir = '/Users/nicholasdwork/.DellEMS/Vols/disk2s1/superResData/';

  switch datacase
    case 1
      imgSliceIndxs = 10;
      t1Dir = [ dataDir, '/brain_15x15mm/t1_irspgr' ];
      t2Dir = [ dataDir, '/brain_15x15mm/t2_fse' ];
      pDir = [ dataDir, '/brain_15x15mm/auc_pyr' ];
      lDir = [ dataDir, '/brain_15x15mm/auc_lac' ];
      bDir = [ dataDir, '/brain_15x15mm/auc_bic' ];
      pds = 30;
      lds = 30;
      bds = 30;

    case 2
      imgSliceIndxs = 15;
      t1Dir = [ dataDir, '/multi_res_brain/t1_spgr' ];
      t2Dir = [ dataDir, '/multi_res_brain/t2_fse' ];
      pDir = [ dataDir, '/multi_res_brain/auc_pyr' ];
      lDir = [ dataDir, '/multi_res_brain/auc_lac' ];
      bDir = [ dataDir, '/multi_res_brain/auc_bic' ];
      pds = 10;
      lds = 10;
      bds = 10;

    case 3
      imgSliceIndxs = 20;
      t1Dir = [ dataDir, '/multi_res_brain/t1_spgr' ];
      t2Dir = [ dataDir, '/multi_res_brain/t2_fse' ];
      pDir = [ dataDir, '/multi_res_brain/auc_pyr' ];
      lDir = [ dataDir, '/multi_res_brain/auc_lac_raw_fov' ];
      bDir = [ dataDir, '/multi_res_brain/auc_bic_raw_fov' ];
      pds = 10;
      lds = 20;
      bds = 30;

    case 4
      imgSliceIndxs = 20:110;
      t1Dir = [ dataDir, '/super_res_for_nick/11_T1W' ];
      t2Dir = [ ];
      pDir = [ dataDir, '/super_res_for_nick/AUC_pyruvate' ];
      lDir = [ dataDir, '/super_res_for_nick/AUC_lactate' ];
      bDir = [ dataDir, '/super_res_for_nick/AUC_bicarbonate' ];
      pds = 10;
      lds = 20;
      bds = 30;
  end

end
