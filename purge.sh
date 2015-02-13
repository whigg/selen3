echo
echo --------------------------------------------
echo Removing undesired files from home directory  
echo --------------------------------------------
echo 
echo Carefully answer the following y/n questions: 
echo
/bin/rm -i  *.exe
/bin/rm -i *.grd
/bin/rm -i *.jun*  
/bin/rm -i *.log
/bin/rm -i *.tmp 
/bin/rm -i tmp*
/bin/rm -i *.cpt*
/bin/rm -i rect*.dat 
/bin/rm -i tmp*
echo ----------------------------------
echo Tarring the SELEN files for backup  
echo ----------------------------------
echo All F90 and f90  files...
tar cvf a.tar *.F90 
tar rvf a.tar *.f90 


echo
echo ====== Script for compilation of FMLIB programs ======
tar rvf a.tar FMLIB_compile.sh

echo
echo ====== The FMLIB programs ======
tar rvf a.tar *FM*.f90 

echo
echo ====== The ALMA directory ======
tar rvf a.tar ./ALMA 

echo
echo ====== The rheology repository /VSC ======
tar rvf a.tar ./VSC	

echo
#echo ====== The Ice models repository ICE-MODELS ======
#tar rvf a.tar ./ICE-MODELS 

echo
echo ====== The utilities repository /UTILITIES ======
tar rvf a.tar ./UTILITIES 

#echo
#echo ====== The Love numbers repository /LOVE_DEPOSIT ======
#tar rvf a.tar ./LOVE_DEPOSIT

echo
echo ====== The data repository /DATA ======
tar rvf a.tar ./DATA/*.dat
tar rvf a.tar ./DATA/*.txt

echo
echo ====== The config files ======
tar rvf a.tar config*.dat 
tar rvf a.tar config*.f90 

echo
echo ====== The make file makeselen.sh ======
tar rvf a.tar makeselen.sh 

echo
echo ====== This script purge.sh ======
tar rvf a.tar purge.sh 

echo
echo ---
echo Compressing all the data...
echo ---
tar cvfz selen_323_04jun2014.tgz a.tar 
/bin/rm a.tar  
echo
echo ---
echo Bye  
echo ---





