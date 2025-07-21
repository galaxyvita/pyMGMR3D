FC = gfortran
FCFLAGS = -ggdb -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid -finit-real=inf

LL_Base = ../LightningPros
Lib_Base = ../Library/bin

#LL_Base=~/Github/pyMGMR3D/LightningPros/
#Lib_Base=~/Github/pyMGMR3D/Library/bin/
F_NLSOL = $(LL_Base)
F_AntFie = $(LL_Base)

EXECUTE = MGMR3D_fit-v5
DEPENDENCIES = MGMR3D-v4.f90 MGMR3D_Fit_RadioFoot-v4.f90 MGMR3D_shower-v5.f90 MGMR3D_analyse-v5.f90 MGMR3D_SetParams-v4.f90 MGMR3D_RFootPars-v5.f90 \
 MGMR3D_BA-v4.f90 MGMR3D_subr.f90 MGMR3D_spline.f90 MGMR3D_FFT.f90  
LIBRARY = $(Lib_Base)/libfftpack5.1d.a  # FFT library, double precision

all: $(EXECUTE)
OBJECTS = $(EXECUTE).o
SOURCES = $(EXECUTE).f90

$(EXECUTE):  nl2sol.o AntFuncCnst.o AntFunct.o $(OBJECTS)
	$(FC) $(FCFLAGS) -o $(EXECUTE) $(OBJECTS) nl2sol.o AntFuncCnst.o AntFunct.o -lm $(LIBRARY)
nl2sol.o: $(F_NLSOL)/nl2sol.f90
	gfortran -c $(F_NLSOL)/nl2sol.f90 -o nl2sol.o
AntFuncCnst.o: $(F_AntFie)/AntFuncCnst.f90 MGMR3D_RFootPars-v5.o
	gfortran -c $(F_AntFie)/AntFuncCnst.f90 -o AntFuncCnst.o
AntFunct.o: $(F_AntFie)/AntFunct.f90 AntFuncCnst.o
	gfortran -c $(F_AntFie)/AntFunct.f90 -o AntFunct.o
MGMR3D_RFootPars-v5.o: MGMR3D_RFootPars-v5.f90
	gfortran -c MGMR3D_RFootPars-v5.f90 -o MGMR3D_RFootPars-v5.o
$(OBJECTS): $(SOURCES) $(DEPENDENCIES)
	$(FC) $(FCFLAGS) -c $(SOURCES)

clean:
	rm -f *.mod *.o $(EXECUTE)
