#!/bin/bash
#
#set source = fftpack5.1d  in  ${FFTPackBase}
# Make sure the folder  ${FFTPackBase}/bin/  exists
source  ${LL_Base}/ShortCuts.sh

cd ${LL_src}
make -f ${LL_scripts}/f90split.make

cd ${FFTPackBase}/FFTPACK/
mkdir temp
cd temp
${LL_bin}/f90split ../fftpack5.1d.f90
#
for FILE in `ls -1 *.f90`;
do
  gfortran -fallow-argument-mismatch -c -O3 $FILE
  if [ $? -ne 0 ]; then
    echo "Errors compiling " $FILE
    exit
  fi
done
rm *.f90
#
ar cr libfftpack5.1d.a *.o
rm *.o
rm *.out
#
mv libfftpack5.1d.a ${FFTPackBase}/bin
cd ..
rmdir temp
#
echo "Library installed as ${FFTPackBase}/bin/libfftpack5.1.a."
