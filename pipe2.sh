cp $PIPE1/data/products/CIRADA_VLASS*_table3_subtile_info*.csv $PIPE2/test_data/test_subtiles.csv
cp $PIPE1/data/products/VLASS*_UOFM_*_Catalogue_*.csv $PIPE2/test_data/test_pybdsf_catenator_out.csv
cd $PIPE2
nohup python3 component_table/vlass_compcat_vlad_stage2_yg.py test_data/test_pybdsf_catenator_out.csv  test_data/test_subtiles.csv  other_survey_data/CATALOG41.FIT other_survey_data/first_14dec17.fits > step1.out 
nohup python3 host_table/vlass_iso_and_cd_finding_v2.py VLASS_components.csv > step2.out 
cd host_table
mkdir LR_output
nohup python3 vlass_uw_lr_v2.py ../VLASS_source_candidates.csv ../other_survey_data/unWISE_coad_directory.csv > step3.out 
nohup python3 stack_matches_v1.py ../VLASS_source_candidates.csv ../VLASS_components.csv > step4.out 
cd ..
mv host_table/VLASS_table* .
nohup python3 finalise_cat/VLASSQL1CIR_catalogue_finalise.py VLASS_table1_components.csv VLASS_table2_hosts.csv test_data/test_subtiles.csv > step5.out 
python3 duplicate_ridder.py
