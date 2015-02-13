ifort   FMSAVE.f90  -c  -o FMSAVE.o
ifort   FMZM90.f90  -c  -o FMZM90.o
ifort   FM.f90  -O  -c  -o FM.o
ifort   TestFM.f90  -c  -o TestFM.o
ifort   TestFM.o    FMSAVE.o    FMZM90.o     FM.o  -o TestFM
ifort   SampleFM.f90  -c   -o  SampleFM.o
ifort   SampleFM.o    FMSAVE.o     FMZM90.o      FM.o -o SampleFM
