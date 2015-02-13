echo  
echo FMLIB_compile.sh purges this directory from previously existing FMLIB modules
echo and compiles the FMLIB programs needed for TABOO ======= GS May 2010 =======
echo 
echo 1- Refreshing the directory... 
echo
#
# --- Revised DM Nov 2012 - the compile command is now passed by the
# --- caller script through the FORTCMD variable 
#
if [ -f ./fmvals.mod ]
 then
 echo Deleting file fmvals.mod
 /bin/rm ./fmvals.mod 
fi
#
if [ -f ./fmzm.mod ]
 then
 echo Deleting file fmzm.mod 
 /bin/rm ./fmzm.mod 
fi
#
if [ -f ./test_vars.mod ]
 then
 echo Deleting file test_vars.mod
 /bin/rm ./test_vars.mod 
fi
#
if [ -f ./fmzm_1.mod ]
 then
 echo Deleting file fmzm_1.mod
 /bin/rm ./fmzm_1.mod
fi
#
if [ -f ./fmzm_2.mod ]
 then
 echo Deleting file fmzm_2.mod
 /bin/rm ./fmzm_2.mod
fi
#
if [ -f ./fmzm_3.mod ]
 then
 echo Deleting file fmzm_3.mod
 /bin/rm ./fmzm_3.mod
fi
#
if [ -f ./fmzm_4.mod ]
 then
 echo Deleting file fmzm_4.mod
 /bin/rm ./fmzm_4.mod
fi
#
if [ -f ./fmzm_5.mod ]
 then
 echo Deleting file fmzm_5.mod
 /bin/rm ./fmzm_5.mod
fi
#
if [ -f ./fmzm_6.mod ]
 then
 echo Deleting file fmzm_6.mod
 /bin/rm ./fmzm_6.mod
fi
#
if [ -f ./fmzm_7.mod ]
 then
 echo Deleting file fmzm_7.mod
 /bin/rm ./fmzm_7.mod
fi
#
if [ -f ./fmzm_8.mod ]
 then
 echo Deleting file fmzm_8.mod
 /bin/rm ./fmzm_8.mod
fi
#
if [ -f ./fmzm_9.mod ]
 then
 echo Deleting file fmzm_9.mod
 /bin/rm ./fmzm_9.mod
fi
echo done!
echo
echo 2- Compiling the FMLIB Fortran programs...
echo
#
echo FMSAVE
$FORTCMD FMSAVE.f90  -c    -o FMSAVE.o
#
echo FMZM90
$FORTCMD FMZM90.f90  -c    -o FMZM90.o
echo FM
$FORTCMD FM.f90      -c    -o FM.o
echo TestFM
$FORTCMD TestFM.f90  -c    -o TestFM.o
$FORTCMD TestFM.o    FMSAVE.o    FMZM90.o     FM.o   -o TestFM
echo SampleFM
$FORTCMD SampleFM.f90   -c -o  SampleFM.o
$FORTCMD SampleFM.o    FMSAVE.o     FMZM90.o      FM.o     -o SampleFM
echo done!
echo






