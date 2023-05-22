EXE="dipole_sq_mod.exe"
rm $EXE
#gfortran -g -fbacktrace  mod_param.F90 ./900_rand_dir/*.F* class_zipped_v3.F90 class_ak_sum_v1.F90 class_pm_sum_v1.F90 mod_g12.F90  prog_ph_12d_int_v1.F90   -o $EXE
#gfortran -Ofast  mod_param.F90 ./900_rand_dir/*.F* class_zipped_v5.F90  sq_mod_dipole.F90  prog_dipole.F90   -o $EXE
gfortran -g -fbacktrace  mod_param.F90 ./900_rand_dir/*.F* class_zipped_v5.F90  sq_mod_dipole.F90  prog_dipole.F90   -o $EXE
