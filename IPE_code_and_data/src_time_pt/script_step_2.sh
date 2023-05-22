EXE="step_2.exe"
rm $EXE
#gfortran -g -fbacktrace  mod_param.F90 ./900_rand_dir/*.F* class_zipped_v3.F90 class_ak_sum_v1.F90 class_pm_sum_v1.F90 mod_g12.F90  prog_ph_12d_int_v1.F90   -o $EXE
#gfortran -Ofast
gfortran -g -fbacktrace -fcheck=all -ffpe-trap=invalid  mod_param.F90 ./900_rand_dir/*.F* ./900_eigen_dir/*.F* ./mod_complex_v1.F90 ./mod_step_2_v1.F90  prog_step_2_v1.F90 -llapack -lblas   -o $EXE
