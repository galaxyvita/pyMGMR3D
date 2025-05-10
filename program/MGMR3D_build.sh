#!/bin/bash
# 
# 
#read -rsp $'Press enter to continue...\n'
# gfortran -c ~/NumLib/LSQ/nl2sol.f90 
# cp nl2sol.o ~/NumLib/LSQ/nl2sol.o

   export MG_Base="/Users/users/scholten/MGMR3D"
#ProgDir="~/../olaf/"
FFD="${MG_Base}/Program"
cd "${FFD}"
#echo ${FFD}
#
make -f "${FFD}/MGMR3D_fit-makefile-v5.mak"
# rm  .mod
#read -rsp $'Press enter to continue...\n'
