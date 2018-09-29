#!/bin/bash
source ~/.profile
cp ./main/source/main_lv1.f90 './main/'
cp ./main/source/main_lv2.f90 './main/'
cp ./main/source/main_lv3.f90 './main/'

read -p "size of the model(n):" n
echo "n=$n"
read -p "size of the model(m):" m
echo "m=$m"
read -p "source position(si):" si
echo "si=$si"
read -p "source position(sj):" sj
echo "sj=$sj"
read -p "radius(lv1):" r1
echo "r1=$r1"
read -p "radius(lv2):" r2
echo "r2=$r2"
read -p "radius(lv3):" r3
echo "r3=$r3"

#parameter(n=#size_of_n#,m=#size_of_m#,si=#pi_of_source#,sj=#pj_of_source#,radius=#radius#)
n=$(($n*2))
m=$(($m*2))
si=$(($si*2))
sj=$(($sj*2))


sed -i '' 's/n=#size_of_n#/n='"$n"'/g' ./main/main_lv1.f90
sed -i '' 's/m=#size_of_m#/m='"$m"'/g' ./main/main_lv1.f90
sed -i '' 's/si=#pi_of_source#/si='"$si"'/g' ./main/main_lv1.f90
sed -i '' 's/sj=#pj_of_source#/sj='"$sj"'/g' ./main/main_lv1.f90
sed -i '' 's/radius=#radius#/radius='"$r1"'/g' ./main/main_lv1.f90



n=$(($n/2))
m=$(($m/2))
si=$(($si/2))
sj=$(($sj/2))

sed -i '' 's/n=#size_of_n#/n='"$n"'/g' ./main/main_lv2.f90
sed -i '' 's/m=#size_of_m#/m='"$m"'/g' ./main/main_lv2.f90
sed -i '' 's/si=#pi_of_source#/si='"$si"'/g' ./main/main_lv2.f90
sed -i '' 's/sj=#pj_of_source#/sj='"$sj"'/g' ./main/main_lv2.f90
sed -i '' 's/radius=#radius#/radius='"$r2"'/g' ./main/main_lv2.f90


n=$(($n/2))
m=$(($m/2))
si=$(($si/2))
sj=$(($sj/2))


sed -i '' 's/n=#size_of_n#/n='"$n"'/g' ./main/main_lv3.f90
sed -i '' 's/m=#size_of_m#/m='"$m"'/g' ./main/main_lv3.f90
sed -i '' 's/si=#pi_of_source#/si='"$si"'/g' ./main/main_lv3.f90
sed -i '' 's/sj=#pj_of_source#/sj='"$sj"'/g' ./main/main_lv3.f90
sed -i '' 's/radius=#radius#/radius='"$r3"'/g' ./main/main_lv3.f90




####
cd /Users/niyiyu/Documents/Fortran/FMM/main/
matlab -nodesktop -nosplash -nojvm -r "run /users/niyiyu/documents/Fortran/fmm/main/matlab/increase_v.m;quit;"
gfortran ./main_lv1.f90 -o FMMprogram_lv1.out
echo "Level 1 run time:" >Time_log.txt
echo -e >>Time_log.txt
{ time -p ./FMMprogram_lv1.out; } 2>>Time_log.txt
rm FMMprogram_lv1.out

matlab -nodesktop -nosplash -nojvm -r "run /users/niyiyu/documents/Fortran/fmm/main/matlab/decrease_travel_time.m;run /users/niyiyu/documents/Fortran/fmm/main/matlab/origin.m;quit;"
gfortran ./main_lv2.f90 -o FMMprogram_lv2.out
echo "Level 2 run time:" >>Time_log.txt
echo -e >>Time_log.txt
{ time -p ./FMMprogram_lv2.out; } 2>>Time_log.txt
rm FMMprogram_lv2.out


matlab -nodesktop -nosplash -nojvm -r "run /users/niyiyu/documents/Fortran/fmm/main/matlab/decrease_travel_time.m;run /users/niyiyu/documents/Fortran/fmm/main/matlab/decrease_v.m;quit;"
gfortran ./main_lv3.f90 -o FMMprogram_lv3.out
echo "Level 3 run time:" >>Time_log.txt
echo -e >>Time_log.txt
{ time -p ./FMMprogram_lv3.out; } 2>>Time_log.txt
rm FMMprogram_lv3.out
rm band.txt
rm main_lv1.f90
rm main_lv2.f90
rm main_lv3.f90

