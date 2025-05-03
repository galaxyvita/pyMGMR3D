#!/bin/bash
# 
# 
#read -rsp $'Press enter to continue...\n'
# gfortran -c ~/NumLib/LSQ/nl2sol.f90 
# cp nl2sol.o ~/NumLib/LSQ/nl2sol.o

export RunFolder=$(pwd)

Home="/home/olaf/"
export ProgDir="${Home}MGMR3D/Program/"
export AntennaFun="/home/olaf/LMA-fit/LMA2019/AntenFunct/"
echo 'Program folder=' ${ProgDir}
echo 'run folder=' ${RunFolder}
echo 'AntennaFunction=' ${AntennaFun}
#
cd ${ProgDir}
make -f "MGMR3D_fit-makefile-v5.mak"
cd ${RunFolder}
${ProgDir}MGMR3D_fit-v5 <TestCase.in
# read -rsp $'Press enter to continue...\n'
 cd plot
#gle -d pdf PulseNuall.gle
#gle -d pdf pulseNud.GLE
gle -d pdf -o ../FitStokes.pdf ${ProgDir}FitStokes.GLE ${RunFolder}/plot/FitResult
gle -d jpg -r 200 -o ../FitStokes-map.jpg ${ProgDir}FitStokes-map.GLE ${RunFolder}/plot/
gle -d pdf -o ../sh-current.pdf ${ProgDir}sh-current.GLE ${RunFolder}/plot/
#cp FitStokes.pdf ../FitStokes.pdf
#cp FitStokes-map-v2.pdf ../FitStokes-map-v2.pdf
#cp sh-current.pdf ../sh-current.pdf
exit

# nohup ./MGMR3D_Gia.sh  >MGMR3D_Gia.log 2>&1  &   
