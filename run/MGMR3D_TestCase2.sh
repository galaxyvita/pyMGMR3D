#!/bin/bash
# 
# 
#read -rsp $'Press enter to continue...\n'
# gfortran -c ~/NumLib/LSQ/nl2sol.f90 
# cp nl2sol.o ~/NumLib/LSQ/nl2sol.o

export RunFolder=$(pwd)
export MG_Base="/Users/users/scholten/MGMR3D"
export LL_Base=/Users/users/scholten/LOFLI
export ProgDir="${MG_Base}/Program"
export AntennaFun="${LL_Base}/AntenFunct/"
# or
#   export AntennaFun="${LL_Base}/AntenFunct/v2-"


echo 'Program folder=' ${ProgDir}
echo 'run folder=' ${RunFolder}
echo 'AntennaFunction=' ${AntennaFun}
#
cd ${ProgDir}
make -f "MGMR3D_fit-makefile-v5.mak"
cd ${RunFolder}

${ProgDir}/MGMR3D_fit-v5 <TestCase2.in

 cd plot
#gle -d pdf PulseNuall.gle
#gle -d pdf pulseNud.GLE
gle -d pdf -o ../FitStokes.pdf ${ProgDir}/FitStokes.gle ${RunFolder}/plot/FitResult
gle -d jpg -r 200 -o ../FitStokes-map.jpg ${ProgDir}/FitStokes-map.gle ${RunFolder}/plot/
gle -d pdf -o ../sh-current.pdf ${ProgDir}/sh-current.gle ${RunFolder}/plot/
#cp FitStokes.pdf ../FitStokes.pdf
#cp FitStokes-map-v2.pdf ../FitStokes-map-v2.pdf
#cp sh-current.pdf ../sh-current.pdf
exit

# nohup ./MGMR3D_Gia.sh  >MGMR3D_Gia.log 2>&1  &   
