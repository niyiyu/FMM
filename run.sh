#!/bin/bash
source ~/.profile
cd /Users/niyiyu/Documents/Fortran/FMM/main/
matlab -nodesktop -nosplash -nojvm -r "run /users/niyiyu/documents/Fortran/fmm/main/matlab/increase_v.m;quit;"
gfortran ./main_5.f90 -o FMMprogram_5.out
./FMMprogram_5.out
rm FMMprogram_5.out

matlab -nodesktop -nosplash -nojvm -r "run /users/niyiyu/documents/Fortran/fmm/main/matlab/decrease_travel_time.m;run /users/niyiyu/documents/Fortran/fmm/main/matlab/origin.m;quit;"
gfortran ./main_10.f90 -o FMMprogram_10.out
./FMMprogram_10.out
rm FMMprogram_10.out


matlab -nodesktop -nosplash -nojvm -r "run /users/niyiyu/documents/Fortran/fmm/main/matlab/decrease_travel_time.m;run /users/niyiyu/documents/Fortran/fmm/main/matlab/decrease_v.m;quit;"
gfortran ./main_20.f90 -o FMMprogram_20.out
./FMMprogram_20.out
rm FMMprogram_20.out

#rm FMMprogram_20.out

