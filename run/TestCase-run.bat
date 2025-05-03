set "LIBRARY=-lm C:\OlafsUtil\NumLib\bin\libFFTPack-d.a"
set Base="C:\Users\Olaf Scholten\Documents\AstroPhys\Lightning\"
set "AntennaFun=%Base%Imaging\LMA\LMA2019\AntenFunct\"
set ProgDir=%Base%MGMR3D\Program\

set "RunFolder=%cd%"
echo current directory: %RunFolder%
echo Program directory: %ProgDir%

cd %ProgDir%
set "FCFLAGS = -ggdb -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid -finit-real=inf"
set FitStokes="%cd%\FitStokes.GLE"
set map="%cd%\FitStokes-map.GLE"
set sh-current="%cd%\sh-current.GLE"
set t-tret="%cd%\t-tret.GLE"
gfortran  %FCFLAGS% -o MGMR3D_fit MGMR3D_fit-v5.f90 %LIBRARY%

cd %RunFolder%
:: echo on
@..\program\MGMR3D_fit <TestCase.in
:: pause
cd plot
:: echo sh-current: %sh-current%
 call gle /d pdf /o ../sh-current.pdf %sh-current% "%cd%/"
:: echo FitStokes: %FitStokes%
 call gle /d pdf /o ../FitStokes.pdf %FitStokes% "%cd%/FitResult"
:: echo map: %map%
 call gle /d jpg /r 200 /o ../FitStokes-map.jpg %map% "%cd%/"
:: echo (t vs t_retarded): %t_tret%
:: call gle /d pdf /o ../t-tret.pdf %t-tret% "%cd%/"
 pause
exit
