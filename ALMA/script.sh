g95 FMSAVE.f90  -c    -o FMSAVE.o
g95 FMZM90.f90  -c     -o FMZM90.o
g95 FM.f90  -O  -c     -o FM.o
g95 TestFM.f90  -c     -o TestFM.o
g95 TestFM.o    FMSAVE.o    FMZM90.o     FM.o   -o TestFM
g95 SampleFM.f90  -c     -o  SampleFM.o
g95 SampleFM.o    FMSAVE.o     FMZM90.o      FM.o     -o SampleFM
