! 
! *****************************
!  This is program "config.f90" 
! *****************************
!
! Last modified GS & FC July 8, 2008  = INTEL PORT = V 2.6 
! Also touched a number of times during July, 2008-
! Re-touched on February 6 2009 for version 2.7 (on the way back from Liverpool) 
! Reviewed GS & FC July 2009 - v3.1-  "Varying coastlines" (reprise!) 
! Reviewed GS & FC July 2009 - v3.1-  "ALMA coupling"
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
! Also revised April 2010 by GS for completing the g95 version in 3.1 
! Updated on April 13, 2010 for degree 1 implementation
! Updated on May 16, 2010 (eb* files move to ALMA love number repository)
! Modified June 2010 - Reference A and GRA for Free air Gravity anomalies (data.inc)
! Modified GS July 2010: Refreshed Stokes coefficient modules & elastic rebound
! GS July 2010: External elastic Love numbers switch included 
! GS July '10 continued: Tidal Love numbers analysis...(i.e. two possible kinds ...)
! ... of TABOO analysys are now possible, one is standard the other is multi-precision). 
! GS July 2010: Polar Motion Transfer Function module added 
! *** Revised GS July 2010 - g95 - Double precision implementation
! *** Revised GS Aug 2010 - polar motion modules 
! *** Revised GS Nov 2010 - Rotational feedback 
! *** Revised GS & GR Nov 2010 - Implementation of MED1
! *** Revised GS Dec 2010 - Implementation Geodetic sites for elastic rebound
! *** Revised GS Dec 2010 - RSL database of type "3" (for data from Dorit Sivan) 
! === Revised GS Jan 2011 - Reference Frame for Love numbers of degree 1 
! === Revised GS & GG Apr 2011 - Including ANT elastic ice models
! === Revised GG Apr 25 2011 - for the ANT component  
! === Revised GG May 03 2011 - for the ANT component  
! === Revised GG May 10 2011 - New format for the tide gauge analysis   
! === Revised GG May 13 2011 - New format for the elareb analysis  
! === Revised GS & FC May 21 2011 - DELTA parameter in ice history     
! *** Revised GS for the rotational part (external Love numbers) <<<<====== IN PROGRESS 
! ::: Revised GS & FC (implementation of "Transient ICE") <<<<====== IN PROGRESS 
! ::: Revised GS (RECTANGULAR ICE for the benchmark SLE II ) <<<<== IN PROGRESS 
! Feb 2012: GS: Implementation of the numerical derivative "in the future" 
! Feb 2012: Multi-step elastic rebound for ice2sea 
! March 3 2012: ESL option is now differently handled for 'ho' and 'pm' ice models
! March 7 2012: Implementation of tidel elastic Love numbers 
! March 8, a bug in GEO.F90 in the data.inc section, now fixed...
! March 12, working to vertical displacement for 'pm' ice models 
! April 12, working to small modifications, ice2sea data are ARRIVED
! April 19 2012, working to GIA-correcting the ice2sea projections (Urbino Hospital)
! June 30, 2012; regular grid for the one-step elastic rebound (EOS paper?)
! July 15, 2012; regular grid for the global maps at presen time...(work with Anny)
! Aug 5, 2012; Equivalent water height (understanding the Chambers et al. 2010 results)
! Nov 6, 2012; Portability options, GMT execution switch and minor optimizations (DM)
! Jul 1, 2013; The "FLO" ice for coupling with her models 
! April 11, 2014. The BENOIT (BEN) files for the collaboration with LEGOS 
! June 21, 2014. The ALPS file for ICELAND for the collaboration with Martin Brader!
! Feb 14, 2015: Implementation of the evolving coastlines WITH rotational feedback
!
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
! Copyright (C) 2008 Giorgio Spada, Florence Colleoni, and Paolo Stocchi 
!
! This file is part of SELEN. 
!  
! SELEN is free software: you can redistribute it and/or modify it under the 
! terms of the GNU General Public License as published by the Free Software 
! Foundation, either version 3 of the License, or at your option any later 
! version. 
!
! SELEN is distributed in the hope that it will be useful, but WITHOUT ANY 
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
! details. 
! 
! You should have received a copy of the GNU General Public License along 
! with SELEN.  If not, see <http://www.gnu.org/licenses/>.
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! --------------------------------------------------------------------------- 
!
!   Program "config.f90" is the interface between the configuration file 
!   "config.dat" and the bash script "selen.sh" which drives the execution 
!   of SELEN. During execution, "config.f90" writes a summary of the SELEN 
!   settings on monitor. The outputs of "config.f90" are: 
!
!   - "data.inc", a Fortran script with general declarations, which serves 
!                 as input for most of the SELEN Fortran programs, 
! 
!   - "selen.sh", a bash shell script that compiles and executes the SELEN 
!                 programs and organizes their outputs according to the 
!		  settings in "config.dat", 
!
!   - several GMT scripts based on the options in "config.dat". 
!
!   Last change GS & FC February 6 2009 
!
! --------------------------------------------------------------------------- 
!
!
!
INCLUDE "harmonics.f90" 
PROGRAM CONFIG
IMPLICIT NONE 
!
!
! =======================================================
!                 --- Declarations ---
! =======================================================
!
!
! -A single precision ZERO 
      REAL, PARAMETER :: ZERO=0.
      REAL ZEROV(201)
!
!
! -Some bounds  
      INTEGER, PARAMETER :: NALFA=52, MAX_LAYERS=100, NIMP=100, & 
  	                    VERY_LARGE_INTEGER=100000, LARGE_INTEGER=1000, & 
			    RES_MIN=12, NRG=10, NMAX_3D_REGIONS=10, RES_MAX=48, & 
			    NP_MAX=2*RES_MAX*(RES_MAX-1)*20+12, & 
   			    LN_MAX=1024 		     	
!			   
! -Do-loop indices & other integers 
      INTEGER I, J, JJ, K, N, II, NP		
!
! -Some counters (mainly counters)
      INTEGER HEADER_LINES, N_HEADER_LINES, ICE_FLAG, IMUL, NLINES, NNV, & 
      	      NDEGREE, NRESOLUTION, NOUT, NRSL, NRSLC, NICE, & 
	      NTIDEGAUGES, NRSL_DATA, LEN_ICE, LEN_RSL, & 
	      LEN_VISCO, LEN_TGAUGES, LEN_RSLC, LEN_GEOD, NGEOD, & 
	      LEN_3D_REGIONS, N_3D_REGIONS, NP_3D_REGIONS(NMAX_3D_REGIONS), & 
	      LEN_TRED_REGIONS_NAME(NMAX_3D_REGIONS), & 
	      LEN_FILE_TOPO, LEN_FILE_PXTOPO, ITER_CI, NNINC, DEN, IREB, & 
	      LEN_FILE_TG_ELA_MS, &
	      NNTHREAD,NNTASK, LEN_TMP
!	    
! -A few real variables 
      REAL*4 DELTA	  
      REAL*8 LATS, LONS
      REAL*4 RJUNK, RSL_DATUM, D_RSLDATUM, & 
             TIME_RSL, TIME_RSL1, TIME_RSL2, TIME_BPC, &
	     LN_H, LN_L, LN_K        	     
!
! -A string useful for compilation commands
!      CHARACTER*8, PARAMETER :: & 
!      Compile = "g95 -O2 "
! -Sequential and parallel compilation commands
      CHARACTER*65 :: CompileSeq, CompileSmp, CompileMpi
!
! -MPI execution prefix
      CHARACTER*100 :: MpiRunCmd      
!
! -The alphabet 
      CHARACTER*1 ALFA(NALFA)
      DATA ALFA/'A','B','C','D','E','F','G','H','I','J','K','L','M', & 
                'N','O','P','Q','R','S','T','U','V','Y','Y','Z','W', & 
	        'a','b','c','d','e','f','g','h','i','j','k','l','m', & 		   	     
		'n','o','p','q','r','s','t','u','v','y','y','z','w'/
!
! -Strings for basic SELEN settings variables
      CHARACTER*10  RESOLUTION, DEGREE, TITLICE, DEGREE_ST_MIN, DEGREE_ST_MAX  
      CHARACTER*10  TIME_BPCC, MIN_RSLC, MAX_RSLC, RSL_INT
      CHARACTER*10  LONMINC, LONMAXC, LATMINC, LATMAXC
      CHARACTER*10  RADIUS_ZOF
      CHARACTER*10  ICE_DELTA 
      CHARACTER*1   MODE
      CHARACTER*2   ITER, ITER_REV             ! New December 2012
      CHARACTER*10  ITER_C  
      CHARACTER*3   NINC, NV, CODE
      CHARACTER*3   MINIMUM_PERIOD  
      CHARACTER*50  NAME_OF_REGION
      CHARACTER*100 WDIR
      CHARACTER*10  LATSC10, LONSC10 
      CHARACTER*3   NTHREAD, NTASK
      CHARACTER*10   RES_REG_GRC, RES_REG_ANC 
!
! -File and GMT scripts names 
      CHARACTER*20 & 
      SH_FILE, 	      ICE_FILE,            ICE_FILENAME,   &
      SHOF_FILE,      FILE_GMT,                            & 
      FILE1_GMT,      FILE2_GMT,           FILE_CAP,       & 
      RSLC_LONLAT_FILE, &       
      FILE_TOPO,      FILE_PXTOPO,         LNAP_FILENAME,  & 
      FILE_GMT_G,     FILE2_GMT_R,         FILE_PXTABLE   
      CHARACTER*40 SHAPE_FILE   
!
      CHARACTER*30 RSL_FILE,     RSL_DATABASE,        & 
      		   FILE_REGION,  FILE_REGION_LONLAT,  & 
                   TGAUGES_FILE, TGAUGES_DATABASE,    & 
		   FILE_3D, FILE_3DD, & 
		   GEO_DATABASE_ELAREB, GEO_DATABASE, &
		   VISCO_FILE, 			      &
      FILE_3D_REGIONS,  TRED_REGIONS_DATABASE,        & 
      NAME_3D_REGIONS(NMAX_3D_REGIONS),               & 
      PIX_3D_FILENAMES(NMAX_3D_REGIONS),              &
      TRED_REGIONS_NAME(NMAX_3D_REGIONS),             & 
      FILE_TG_ELA_MS,                                 &
      DATABASE_TG_ELA_MS
!
      CHARACTER*100 VISCO_FILENAME, ALMA_FILE_MODE_2
      CHARACTER*100 ALMA_FILE_PWA_H, ALMA_FILE_PWA_L, ALMA_FILE_PWA_K
      CHARACTER*100 ALMA_LOGFILE_PWA 
      CHARACTER*100 SHORT_VISCO_FILENAME
      CHARACTER*100 LNAP_FILE 
      CHARACTER*100 TMP_PREFIX      
!                   	
!
! -PAth to ICE-MODELS 
      CHARACTER*100 XFILE
!
! -Rheological parameters 
      CHARACTER*100 VSTRING, LITHO, VSC(MAX_LAYERS)
!
! -Input text strings 
      CHARACTER*100 SS(NIMP)
      CHARACTER*200 LINE, LINEP, ROW, JUNK
      CHARACTER*80  ANOTHERROW
      CHARACTER*10 CJUNK, AJUNK, SS2(2) 
      CHARACTER*1 RSL_DATABASE_FORMAT 
      CHARACTER*1 TGAUGES_DATABASE_FORMAT 
      CHARACTER*1 ELAREB_DATABASE_FORMAT 
      CHARACTER*1 ELAREB_MS_DATABASE_FORMAT 
!
! -Options   
      CHARACTER*1 & 
      OPTION_SH,       OPTION_SF,     OPTION_WI,      OPTION_PX,           OPTION_OF,       & 
      OPTION_RI,       OPTION_LN,     OPTION_GM,      OPTION_TGPLOT,       OPTION_ST,       &
      OPTION_CL,       OPTION_LB,     OPTION_US,      OPTION_OH,           OPTION_OR,       &
      OPTION_RSLDB,    OPTION_RSLP,   OPTION_RSL,     OPTION_RSLZ,         OPTION_RSLMF,    & 
      OPTION_RSLSCA,   OPTION_TG,     OPTION_TGSCA,   OPTION_BIN,          OPTION_ESL,      & 
      OPTION_RSLC,     OPTION_TGA,    OPTION_RSLA,    OPTION_RSLTAB,       OPTION_OFDV,     & 
      OPTION_RM(0:NRG),OPTION_ROF,    OPTION_3D,      OPTION_3D_REGIONS,   OPTION_TOPO,     & 
      OPTION_PXTOPO,   OPTION_NM,     OPTION_PW,      OPTION_LOVE_NUMBERS, OPTION_PTMAP,    & 
      OPTION_PWA,      OPTION_DEG1,   OPTION_TLOVE,                                         &
      OPTION_DOTS,     OPTION_DOTU,   OPTION_DOTN,    OPTION_DOTFA,        OPTION_DOTSS,    &
      OPTION_DOTLOI,   OPTION_DOTLOO, OPTION_DOTLOT,                       OPTION_REB,      &
      OPTION_REB_GR,   OPTION_REB_GG, OPTION_REB_AR,  OPTION_REB_AG,       OPTION_REB_SM,   &
      OPTION_LNAP,     OPTION_TIDAL,  OPTION_PMTF,    OPTION_PMD,          OPTION_PMPLOT,   &
      OPTION_RFB,      OPTION_3D_REB, OPTION_DOTG,    OPTION_TIDAL_EXT,    OPTION_TLOVE_EXT,&
      OPTION_DER,      OPTION_ELA_MS, OPTION_DOTW,      & 
      OPTION_GM_ELA_MS_S_VAR,  OPTION_GM_ELA_MS_S_DOT,  & 
      OPTION_GM_ELA_MS_N_VAR,  OPTION_GM_ELA_MS_N_DOT,  & 
      OPTION_GM_ELA_MS_U_VAR,  OPTION_GM_ELA_MS_U_DOT,  & 
      OPTION_GM_ELA_MS_SE_VAR, OPTION_GM_ELA_MS_SE_DOT, &
      OPTION_TG_ELA_MS, & 
      OPTION_ESL_PM,    & 
      OPTION_GIA_CORR,  & 
      OPTION_XYGRID,    & 
      OPTION_SCALFA, OPTION_GMAPS_XYGRID,    &
      OPTION_SYS, OPTION_OMP, OPTION_MPI, OPTION_NLS, OPTION_GMT, OPTION_NPX
!            
      CHARACTER*2 ICE_TYPE           
!
! -More options   
      CHARACTER*2      OPTION_RFRAME       
!
!  
! -Repository name 
      CHARACTER*12 depot 
!
! -Other character constants 
      CHARACTER*40 for_argument
! 
! -Repository label       
      CHARACTER*4 RUN
      CHARACTER*6 RUNP
!
! -Date and time...
      CHARACTER*20 date, timc   
!
! -Some logical switches 
      LOGICAL LEX, LEXX(4) 
!      
! -Log file
      CHARACTER*30, parameter :: logfile="selen.log"
!
! -Water and ice density 
!  REAL*4,  PARAMETER :: RHO_EARTH=5511.57
      REAL*4,  PARAMETER :: RHO_WATER=1000.00 
      REAL*4,  PARAMETER :: RHO_ICE  = 916.70
!
! -Today 
      REAL*4,  PARAMETER :: time_today=2011. 
!
! -Scaling factors for ice thickness
      REAL*4 SCALFA_NUM
      CHARACTER*100 SCALFA_STRING
!
! ************************************************************************
! ************************************************************************
! ************************************************************************
! ************************************************************************
! ************************************************************************
!                           Execution starts here!
! ************************************************************************
! ************************************************************************
! ************************************************************************
! ************************************************************************
! ************************************************************************
!
!
      Write(*,*) " = = = = = = = = = = = = = = =" 
      Write(*,*) "---- ICE DENSITY is (kg/m^3): ", RHO_ICE
      Write(*,*) " = = = = = = = = = = = = = = =" 
!
      Open(88,file=logfile,status='unknown') 
      call DATE_AND_TIME (date,timc)      
      Write(88,*) " "
      Write(88,*) '# This is file "selen.log", created by "config.f90" on ', & 
             date(1:4), '.', date(5:6), '.', date(7:8), ' ', & 
	     timc(1:2), '.', timc(3:4), '.', timc(5:6) 
      Write(88,*) "  "
!
!
! ====================================
!
! 	Part #1 : 
! 	1/1: General settings 
! 	2/1: Output settings
! 	3/1: Organizing the outputs
!
! ====================================
!
!
! 	# Part 1/1: Reading the general settings "config.dat"   
!
!
  open(1,file='config.dat',status='unknown')
!
!
  Write(88,*) ""
  Write(88,*) "+----------------------------------------------+"
  Write(88,*) "|              Settings of SELEN               |"
  Write(88,*) "+----------------------------------------------+"
!
!
  DO 100 I=1,LARGE_INTEGER
!
  read(1,'(a200)',end=444) line 
!
!
!
! ###### WORKING DIRECTORY ######  
!
  IF(line(1:3)=='999') THEN 
	call scan_string (line, 1, ss, nout)
	wdir=ss(1) 
	Write(88,*) "The working directory is ", wdir 
	Open (55,file='working-directory.txt',status='unknown') 
	Write(55,'(a100)') wdir 
	close(55) 
  ENDIF
!
!
! ###### SYSTEM SETTINGS ######
!
!
!
! >>>>>>>> MPI setup -- NOT YET IMPLEMENTED !!!!
!
!
  option_mpi='n'
!  IF(line(1:3)=='996') THEN 
!	call scan_string (line, 2, ss, nout)
!	option_mpi=ss(1) 
!	ntask     =ss(2)
!        Call CHAR3_2_INT(ntask, nntask)
!	if(option_mpi=='y') Write(88,*) "MPI is enabled on ",ntask," tasks"
!	If(option_mpi/='y'.and.option_mpi/='n') then 
!		Write(* ,*) "For the MPI switch, only y/n are valid options"
!		Write(88,*) "For the MPI switch, only y/n are valid options"
!		Call stop_config 
!		Stop	
!	Endif
!  ENDIF
!
!
!
!
! >>>>>>>> OpenMP setup -- NOT YET IMPLEMENTED !!!!
!
  IF(line(1:3)=='997') THEN
	call scan_string (line, 2, ss, nout)
	option_omp=ss(1) 
	nthread   =ss(2)
        Call CHAR3_2_INT(nthread, nnthread)
	if(option_omp=='y') Write(88,*) "OpenMP is enabled with ",nthread," threads/task"
	If(option_omp/='y'.and.option_omp/='n') then 
		Write(* ,*) "For the OpenMP switch, only y/n are valid options"
		Write(88,*) "For the OpenMP switch, only y/n are valid options"
		Call stop_config 
		Stop	
	Endif
  ENDIF  
!
!
!
! >>>>>>>> Compiler setup
!
if(line(1:3)=="998") THEN
    call scan_string (line, 1, ss, nout)
    option_sys=ss(1)
!
    If( option_sys == "1" ) then     ! ------------------------------------------- MAC OS X + gfortran
!
       Write(88,*) "Configuring SELEN for gfortran on Mac OS X"
!
       CompileSeq = "gfortran -w -O3 -DGNU -m64 "      
       if( ( option_mpi=='n' ) .and. ( option_omp=='y' ) ) then
          CompileSmp = "gfortran -m64 -w -O3 -fopenmp -Wl,-stack_size,0x10000000 "   
	  CompileMpi = CompileSmp
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='n' ) ) then
          CompileMpi = "mpif90 -m64 -w -O3 -fopenmp -DMPI "   
	  CompileSmp = CompileSeq
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='y' ) ) then
          CompileMpi = "mpif90 -m64 -w -O3 -fopenmp -DMPI -Wl,-stack_size,0x10000000 " 
          CompileSmp = "gfortran -m64 -w -O3 -fopenmp -Wl,-stack_size,0x10000000 "   
       endif       
       if( ( option_mpi=='n' ) .and. ( option_omp=='n' ) ) then
          CompileMpi = CompileSeq
	  CompileSmp = CompileSeq
       endif
!
    ElseIf( option_sys == "2" ) then ! ------------------------------------------- MAC OS X + ifort
!
       Write(88,*) "Configuring SELEN for Intel ifort on Mac OS X"
!
       CompileSeq = "ifort "
       if( ( option_mpi=='n' ) .and. ( option_omp=='y' ) ) then
          CompileSmp = "ifort -openmp -Wl,-stack_size,0x10000000 "
	  CompileMpi = CompileSmp
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='n' ) ) then
          CompileMpi = "mpif90 -DMPI "    
	  CompileSmp = CompileSeq
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='y' ) ) then
          CompileMpi = "mpif90 -DMPI -openmp -Wl,-stack_size,0x10000000 "   
          CompileSmp = "ifort -openmp -Wl,-stack_size,0x10000000 "
       endif       
       if( ( option_mpi=='n' ) .and. ( option_omp=='n' ) ) then
          CompileMpi = CompileSeq
	  CompileSmp = CompileSeq
       endif
!
    ElseIf( option_sys == "3" ) then ! ------------------------------------------- Linux + gfortran
!
       Write(88,*) "Configuring SELEN for gfortran on Linux"
!
       CompileSeq = "gfortran -w -O3 -DGNU "
       if( ( option_mpi=='n' ) .and. ( option_omp=='y' ) ) then
          CompileSmp = "gfortran -w -fopenmp -O3 "   
	  CompileMpi = CompileSmp
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='n' ) ) then
          CompileMpi = "mpif90 -w -O3 -DMPI "  
	  CompileSmp = CompileSeq
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='y' ) ) then
          CompileMpi = "mpif90 -w -O3 -DMPI -fopenmp "     
          CompileSmp = "gfortran -w -fopenmp -O3 "   	  
       endif   
       if( ( option_mpi=='n' ) .and. ( option_omp=='n' ) ) then
          CompileMpi = CompileSeq
	  CompileSmp = CompileSeq
       endif
!    
    ElseIf( option_sys == "4" ) then ! ------------------------------------------- Linux + ifort
!
       Write(88,*) "Configuring SELEN for Intel ifort on Linux"
!
       CompileSeq = "ifort "
       if( ( option_mpi=='n' ) .and. ( option_omp=='y' ) ) then
          CompileSmp = "ifort -openmp "    	   
          CompileMpi = CompileSmp
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='n' ) ) then
          CompileMpi = "mpif90 -DMPI "
	  CompileSmp = CompileSeq
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='y' ) ) then
          CompileMpi = "mpif90 -openmp -DMPI "      
          CompileSmp = "ifort -openmp "    	   
       endif   
       if( ( option_mpi=='n' ) .and. ( option_omp=='n' ) ) then
          CompileMpi = CompileSeq
	  CompileSmp = CompileSeq
       endif
!
    ElseIf( option_sys == "9" ) then ! ------------------------------------------- Linux + ifort
!
       Write(88,*) "Configuring SELEN for Intel ifort on Linux"
!
       CompileSeq = "ifort "
       if( ( option_mpi=='n' ) .and. ( option_omp=='y' ) ) then
          CompileSmp = "ifort -openmp "
	  CompileMpi = CompileSmp    	   
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='n' ) ) then
          CompileMpi = "mpiifort -DMPI "   
	  CompileSmp = CompileSeq   
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='y' ) ) then
          CompileMpi = "mpiifort -openmp -DMPI " 
          CompileSmp = "ifort -openmp "	       
       endif   
       if( ( option_mpi=='n' ) .and. ( option_omp=='n' ) ) then
          CompileMpi = CompileSeq
	  CompileSmp = CompileSeq
       endif
!
    ElseIf( option_sys == "5" ) then ! ------------------------------------------- MinGW 64-bit
!
       Write(88,*) "Configuring SELEN for 64-bit CYGWIN/MinGW"
!
       CompileSeq = "x86_64-w64-mingw32-gfortran -w -O3 -DGNU "
       if( ( option_mpi=='n' ) .and. ( option_omp=='y' ) ) then
          CompileSmp = "x86_64-w64-mingw32-gfortran -openmp -w -O3 "
	  CompileMpi = CompileSmp    	   
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='n' ) ) then
          CompileMpi = "mpiifort -DMPI "   
	  CompileSmp = CompileSeq   
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='y' ) ) then
          CompileMpi = "mpiifort -openmp -DMPI " 
          CompileSmp = "ifort -openmp "	       
       endif   
       if( ( option_mpi=='n' ) .and. ( option_omp=='n' ) ) then
          CompileMpi = CompileSeq
	  CompileSmp = CompileSeq
       endif
!
    ElseIf( option_sys == "6" ) then ! ------------------------------------------- MinGW 32-bit
!
       Write(88,*) "Configuring SELEN for 32-bit CYGWIN/MinGW"
!
       CompileSeq = "i686-w64-mingw32-gfortran -w -O3 -DGNU "
       if( ( option_mpi=='n' ) .and. ( option_omp=='y' ) ) then
          CompileSmp = "i686-w64-mingw32-gfortran -openmp -w -O3 "
	      CompileMpi = CompileSmp    	   
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='n' ) ) then
          CompileMpi = "mpiifort -DMPI "   
	      CompileSmp = CompileSeq   
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='y' ) ) then
          CompileMpi = "mpiifort -openmp -DMPI " 
          CompileSmp = "ifort -openmp "	       
       endif   
       if( ( option_mpi=='n' ) .and. ( option_omp=='n' ) ) then
          CompileMpi = CompileSeq
	      CompileSmp = CompileSeq
       endif
!
    ElseIf( option_sys == "7" ) then ! ------------------------------------------- MinGW 32-bit
!
       Write(88,*) "Configuring SELEN for g95"
!
       CompileSeq = "g95 -O2 "
       if( ( option_mpi=='n' ) .and. ( option_omp=='y' ) ) then
          	Write(* ,*) "Error: The g95 compiler does not support OpenMP!"
          	Write(88 ,*) "Error: The g95 compiler does not support OpenMP!"
            call stop_config
            stop
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='n' ) ) then
          CompileMpi = "mpif90 -O2 -DMPI "   
	      CompileSmp = CompileSeq   
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='y' ) ) then
          CompileMpi = "mpiifort -openmp -DMPI " 
          CompileSmp = "ifort -openmp "	       
       endif   
       if( ( option_mpi=='n' ) .and. ( option_omp=='n' ) ) then
          	Write(* ,*) "Error: The g95 compiler does not support OpenMP!"
          	Write(88 ,*) "Error: The g95 compiler does not support OpenMP!"
            call stop_config
            stop
       endif
!
!
    Else   
	Write(* ,*) "Error: unknown platform"
	Write(88,*) "Error: unknown platform"
	Call stop_config 
	Stop	       
    EndIf
!
!
!
    if( option_mpi == "y" ) then
       if( option_sys == "9" ) then
           if( option_omp == "n" ) then
               MpiRunCmd = "mpirun -np "//ntask
           else
               MpiRunCmd = "mpirun -ppn 1 -np "//ntask//" -env OMP_NUM_THREADS "//trim(nthread)//" "
           endif
       else
           if( option_omp == "n" ) then
               MpiRunCmd = "mpirun -np "//ntask
           else
               MpiRunCmd = "mpirun -bynode -np "//ntask//" -x OMP_NUM_THREADS="//trim(nthread)//" "
	   endif
       endif
    else
       MpiRunCmd=""
    endif
!
endif
!
!
!
!
! >>>>>>>> Node-local storage 
!
IF(line(1:3)=="995") THEN 
	call scan_string (line, 2, ss, nout)
        option_nls = ss(1) 
        tmp_prefix = ss(2) 
        if((option_mpi=='n') .and. (option_nls=='y')) then
		     Write(88,*) "Node-local storage is available only for MPI jobs"
		     Write(88,*) "Option will be ignored"
		     option_nls='n'
		Endif   
	    if(option_nls=='y') then 
		     Write(88,*) "SELEN will store PX and SH files under: ", trim(adjustl(tmp_prefix))
      	 	 tmp_prefix=trim(adjustl(tmp_prefix))//'/'
        	 len_tmp=len(trim(adjustl(tmp_prefix))) 
      	 	 tmp_prefix="'"//trim(adjustl(tmp_prefix))//"'"
        Endif
ENDIF
!
!
!
!
! ###### Solution of the SLE ######  
!
! >>>>>>  Iterations 
!
IF(line(1:3)=="000") THEN 
	call scan_string (line, 2, ss, nout)
	iter   = ss(1)                      
	mode   = ss(2)	
	
	if(iter(1:1)=='0') iter_rev=iter(2:2)
	if(iter(1:1)/='0') iter_rev=iter(1:2)
	
	
	
!
	if    (mode=='1') then 
		Write(88,*) "The SLE will be solved in mode 1 (fully self-gravitating)"
	elseif(mode=='2') then 
		Write(88,*) "The SLE will be solved in mode 2 (elastic approximation)"
	elseif(mode=='3') then 
		Write(88,*) "The SLE will be solved in mode 3 (eustatic approximation)"	
	elseif(mode=='4') then 
		Write(88,*) "The SLE will be solved in mode 4 (Woodward approximation)"		
	elseif(mode=='5') then 
		Write(88,*) "The SLE will be solved in mode 5 (ice load neglected)"		
	endif					
!	
	If(iter/='0') then 
	Write(88,*) "The SLE will be solved with ", iter, " iteration(s)"
	else
	Write(88,*) "The SLE will *not* be solved"
        endif   
ENDIF
!
!
!
! ###### COASTLINES ######  
!
! >>>>>> Topo file & iterations 
!
IF(line(1:3)=="005") THEN 
	call scan_string (line, 3, ss, nout)
	OPTION_TOPO 	= ss(1) 
!
	FILE_TOPO='./DATA/'//ss(2)
!
	ITER_C          = ss(3)
	call CHAR10_2_INT(iter_c, iter_ci)             	
!	
 	if    (option_topo=='y') then 
		INQUIRE(FILE=file_topo,EXIST=lex)
                If(lex) then 
			len_file_topo=len(trim(adjustl(FILE_TOPO)))
			write(88,*) 'Solving the SLE with VARYING coastlines...'
 			write(88,*) 'The present time TOPO file is: ', trim(adjustl(file_topo))
			write(88,*) 'Number of external iterations: ', trim(adjustl(ITER_C))
		else
			write(88,*) 'File ', trim(adjustl(file_topo)), ' apparently does NOT exist'
			write(88,*) 'Check the input data & try again!'
	        	call stop_config 
	  		Stop		
		endif
!		
	elseif(option_topo=='n') then 		
 		write(88,*) 'Solving the SLE with FIXED coastlines...'
	endif	
ENDIF 
!
!
! >>>>>> Pixelized topography  
!
IF(line(1:3)=="006") THEN 
	call scan_string (line, 2, ss, nout)
        OPTION_PXTOPO = ss(1) 
	FILE_PXTOPO   = ss(2) 
!
	len_file_pxtopo=len(trim(adjustl(FILE_PXTOPO)))
!

	If      (option_topo=='y'.and.option_pxtopo=='y') then 
!	
        	Write(88,*) "A *NEW* pixelization of topography will be performed"
!	
	elseif  (option_topo=='y'.and.option_pxtopo=='n') then 
!
		INQUIRE(FILE=file_pxtopo,EXIST=lex)
		If(lex) then 
			Write(88,*) "Found a pre-built pixelized topography in file ", & 
		            trim(adjustl(FILE_PXTOPO)) 		
		else 
			write(88,*) 'File ', trim(adjustl(FILE_PXTOPO)), ' apparently does NOT exist'
			write(88,*) 'Check the input data & try again!'
	        	call stop_config 
	  		Stop				
		endif 	   
!	
	Endif 
!
ENDIF  
!
!
! ###### Maximum harmonic degree ###### 
!
IF(line(1:3)=="001") THEN 
	call scan_string (line, 1, ss, nout)
	degree=ss(1) 
 	call CHAR10_2_INT(degree, ndegree)
 	write(88,*) 'Maximum degree is: ', ndegree
 	write(88,*) 'Total number of harmonics:', (ndegree+1)*(ndegree+2)/2
ENDIF 
!
!
! ###### Dealing with harmonic degree "1" ###### 
!
IF(line(1:3)=="002") THEN 
!
	call scan_string (line, 2, ss, nout)
	OPTION_DEG1   = ss(1) 	
	OPTION_RFRAME = ss(2) 
!
 	If    (option_deg1=='y') then 
	write(88,*) 'Harmonic degree 1 IS included Yahouh!'	
!
 	If    (option_rframe=='CM') & 
	write(88,*) 'Reference frame of the MASS CENTER of the WHOLE EARTH'	
!
 	If    (option_rframe=='CE') & 
	write(88,*) 'Reference frame of the MASS CENTER of the SOLID EARTH'	
!
 	If    (option_rframe/='CE'.and.option_rframe/='CM') THEN 
		Write(88,*) "Please select a valid option (CE/CM) for the reference frame" 
		Write(*, *) "Please select a valid option (CE/CM) for the reference frame" 
        	CALL STOP_CONFIG 
        Endif
!	
	elseif(option_deg1=='n') then 
	write(88,*) 'Harmonic degree 1 is NOT included'
!
	elseif(option_deg1/='n'.and.option_deg1/='y') then 
	Write(88,*) "Please select a valid option (y/n) for harmonic degree 1" 
	Write(*, *) "Please select a valid option (y/n) for harmonic degree 1" 
        CALL STOP_CONFIG 
	Endif 
!
ENDIF 

!
!
! ###### Ice model ######
!
! >>>>>> Ice file name 
!
IF(line(1:3)=="020") THEN 
	call scan_string (line, 1, ss, nout)
	ice_file=ss(1)
! 
	Write(88,*) 'The ice file is: ', ice_file	
!
ice_flag=0
!
! -------------------------------------------------------------------------- ICE3G 
If    (ice_file=='ice3g.dat'    .or.ice_file=='ice3g_and.dat'.or.& 
       ice_file=='ice3g_ant.dat'.or.ice_file=='ice3g_bal.dat'.or.&
       ice_file=='ice3g_bar.dat'.or.ice_file=='ice3g_bri.dat'.or.& 
       ice_file=='ice3g_gro.dat'.or.ice_file=='ice3g_ice.dat'.or.&  
       ice_file=='ice3g_nam.dat'.or.ice_file=='ice3g_sib.dat')  then 
!
       ice_flag=1 ; ninc='18' ; header_lines=20 ; IMUL=2 ;  titlice='ICE3G'
!
! -------------------------------------------------------------------------- IMED1 
ELSEIF(ice_file=='imed1.dat')                                   then
!
       ice_flag=1 ; ninc='18' ; header_lines=20 ; IMUL=2 ;  titlice='IMED1'
!
! -------------------------------------------------------------------------- DISK 
ELSEIF(ice_file=='disk.dat'.or.ice_file=='disk_on.dat'.or.&
                               ice_file=='disk_off.dat')        then
!
       ice_flag=1 ; ninc='18' ; header_lines=20 ; IMUL=2 ;  titlice='DISK'
!
! -------------------------------------------------------------------------- ICAP 
ELSEIF(ice_file=='icap.dat'    .or.ice_file=='icap_on.dat'.or.&
       ice_file=='icap_off.dat'.or.ice_file=='icap_pm.dat')        then
!
       ice_flag=1 ; ninc='18' ; header_lines=20 ; IMUL=2 ;  titlice='ICAP'
!
! -------------------------------------------------------------------------- GREEN
ELSEIF(ice_file(1:4)=='gree')                                  then 
!
       ice_flag=1 ; ninc='1'  ; header_lines=10  ; IMUL=1 ; titlice='GREEN' 
!
! -------------------------------------------------------------------------- ANTA
ELSEIF(ice_file(1:4)=='anta')                                  then 
!
       ice_flag=1 ; ninc='1'  ; header_lines=10  ; IMUL=1 ; titlice='ANTA' 
!
! -------------------------------------------------------------------------- GLAC
ELSEIF(ice_file(1:4)=='glac')                                  then 
!
       ice_flag=1 ; ninc='1'  ; header_lines=10  ; IMUL=1 ; titlice='GLAC' 
!
! -------------------------------------------------------------------------- ICE1       
ELSEIF(ice_file=='ice1.dat'    .or.ice_file=='ice1_eup.dat'.or. & 
       ice_file=='ice1_gro.dat'.or.ice_file=='ice1_nam.dat')    then 
!
       ice_flag=1 ; ninc='18' ; header_lines=20 ; IMUL=1 ;  titlice='ICE1'
!
! -------------------------------------------------------------------------- ICE5G
ELSEIF(ice_file=='ice5g.dat'    .or.ice_file=='ice5g_and.dat'.or.&
       ice_file=='ice5g_ant.dat'.or.ice_file=='ice5g_fen.dat'.or.& 
       ice_file=='ice5g_gre.dat'.or.ice_file=='ice5g_icl.dat'.or.&
       ice_file=='ice5g_lau.dat'.or.ice_file=='ice5g_nwz.dat'.or.& 
       ice_file=='ice5gA1.dat')  then 
!
       ice_flag=1 ; ninc='21' ; header_lines=24 ; IMUL=1 ;  titlice='ICE5G'		
!
! -------------------------------------------------------------------------- ICE5G26
ELSEIF(ice_file=='ice5g26.dat')  then 
!
       ice_flag=1 ; ninc='26' ; header_lines=29 ; IMUL=1 ;  titlice='ICE5G26'
!
! -------------------------------------------------------------------------- ALPS
ELSEIF (ice_file=='alpst.dat'   .or.  ice_file=='alpsf.dat'  .or. & 
	ice_file=='alpsh.dat'   .or.  ice_file=='alpsc.dat'  .or. & 
	ice_file=='alpsi.dat'   .or.  ice_file=='alpsx.dat')    then 	 
! 
       ice_flag=1 ; ninc='18' ; header_lines=4  ; IMUL=1 ;  titlice='ALPS'
!      ice_flag=1 ; ninc='24' ; header_lines=4   ; IMUL=1 ;  titlice='ALPS'
!      ice_flag=1 ; ninc='21' ; header_lines=4   ; IMUL=1 ;  titlice='ALPS'
!
! -------------------------------------------------------------------------- FLORENCE ICE
ELSEIF (ice_file(1:4)=='flor')    then 	 
! 
       ice_flag=1 ; ninc='41' ;  header_lines=4  ; IMUL=1 ;  titlice='FLOR'

!
ELSEIF (ice_file=='alpsAK01.dat'.or.ice_file=='alpsg.dat') then 	
!
       ice_flag=1 ; ninc='20' ; header_lines=4  ; IMUL=1 ;  titlice='ALPS'
!
! -------------------------------------------------------------------------- IJMOD       
ELSEIF (ice_file=='ij05mod.dat')  then 
!
       ice_flag=1 ; ninc='18' ; header_lines=18 ; IMUL=2 ;  titlice='IJ05MOD'	
!
! -------------------------------------------------------------------------- ANU05 
ELSEIF (ice_file(1:5)=='anu05')   then 
!
       ice_flag=1 ; ninc='30' ; header_lines=26 ; IMUL=1 ;  titlice='ANU05'
!
! -------------------------------------------------------------------------- TRANS
ELSEIF (ice_file(1:4)=='tran')    then 
!
       ice_flag=1 ; ninc='50' ; header_lines=20 ; IMUL=6 ;  titlice='TRANS'
!
! -------------------------------------------------------------------------- RETT
ELSEIF (ice_file(1:4)=='rett') then 
!
       ice_flag=1 ; ninc='50' ; header_lines=20 ; IMUL=6 ;  titlice='RETT'
!
! -------------------------------------------------------------------------- "GRD" files  
ELSEIF (ice_file(1:3)=='grd') then 
!
       ice_flag=1 ; ninc='11' ; header_lines=10 ; IMUL=1 ;  titlice='GRD'
!
! -------------------------------------------------------------------------- "BEN0" files  
ELSEIF (ice_file(1:4)=='beno') then 
!
       ice_flag=1 ; ninc='200' ; header_lines=4 ; IMUL=1 ;  titlice='BENO'

ENDIF   
!
	  if(ice_flag==0) then
 	  Write(88,*) "The ice file apparently does not exist" 
 	  Write(*, *) "The ice file apparently does not exist" 
          call stop_config 
	  Stop
	  Endif
!
! ---- Counts the ice elements listed in the ice file 
	  xfile='./ICE-MODELS/'//trim(adjustl(ice_file))
          call ice_count(xfile, imul, header_lines, nice)
   	  Write(88,*) 'There are ', nice, 'elements in the ice file ', ice_file 
!
	  len_ice=len(trim(adjustl(ice_file))) 
!
	  ice_filename="'"//trim(adjustl(ice_file))//"'"
!
! ----  Completes the configuration of "alma.inc" 
	If(option_pw=='y') then 
		OPEN(120,file=trim(adjustl(wdir))//"/ALMA/alma.inc",status='old')		
		Call CHAR3_2_INT(NINC, NNINC)
		Write(120,*)"integer, parameter :: p= ",    NNINC, " ! Number of time points minus one" 
!Write(120,*)"real, parameter :: m1=0, m2=", NNINC*DELTA, " ! tmin and tmax"
!Write(120,*)"! "		 
!Write(120,*)"! ********* File closed by config.f90 *********"
!Write(120,*)"! "		 		 
!CLOSE(120)	
	Endif	  
!
ENDIF
!
!
!
IF(line(1:3)=="025") THEN 
	call scan_string (line, 1, ss, nout)
!
	ice_type=ss(1)
! 
	If(ice_type/='ho'.and.ice_type/='po'.and.ice_type/='pm') then 
	Write(88,*) "Kind of ice models available:"
	Write(88,*) "<<ho>>: Holocene"
	Write(88,*) "<<po>>: Present, one step"
	Write(88,*) "<<pm>>: Present, multi step"
	Write(88,*) "Invalid ice type" 
	Write(*, *) "Invalid ice type" 
        CALL STOP_CONFIG 
	Endif 
!
	if(ice_type=='ho') Write(88,*) 'The ice type is Holocene'	
	if(ice_type=='po') Write(88,*) 'The ice type is present (one step)'	
	if(ice_type=='pm') Write(88,*) 'The ice type is present (multi-step)'	
ENDIF
!
!
!
IF(line(1:3)=="026") THEN                       ! In progress ... In progress ... In progress ... 
	call scan_string (line, 1, ss, nout)
!
	ice_delta=ss(1)   
!
	Call CHAR10_2_REAL(ice_delta, delta)	
! 
! ....  Any warning here ???????
!
	Write(88,*) 'The ice time step (kyr) is :', delta	
!       Write(* ,*) 'The ice time step (kyr) is :', delta	
!
! ----  Completes the configuration of "alma.inc" 
	If(option_pw=='y') then 
		OPEN(120,file=trim(adjustl(wdir))//"/ALMA/alma.inc",status='old')		
		Call CHAR3_2_INT(NINC, NNINC)
!		Write(120,*)"integer, parameter :: p= ",    NNINC, " ! Number of time points minus one" 
		Write(120,*)"real, parameter :: m1=0, m2=", NNINC*DELTA, " ! tmin and tmax"
		Write(120,*)"! "		 
		Write(120,*)"! ********* File closed by config.f90 *********"
		Write(120,*)"! "		 		 
		CLOSE(120)	
	Endif	  
ENDIF                                           ! In progress ... In progress ... In progress ...
!
!
!
! >>>>>>  "Shape factors" file 
!
IF(line(1:3)=="030") THEN 
	call scan_string (line, 2, ss, nout)
	option_sf  = ss(1)
  	shape_file = ss(2) 
!
	if(option_sf=='y') then 
!		Write(88,*) "*New* shape factors will be computed"
!
        Write(88,*) "A *new* ice SH decomposition will be determined ..."
        Write(88,*) " - this accounts for the ice shape factors AND for the history of ice thickness "
        Write(88,*) " - This computation is time-consuming. Saving the data for future computations  "
        Write(88,*) "   based on the SAME ice and on the maximum harmonic DEGREE is higly recommended."
!
        else
!		Write(88,*) "Found a pre-built shape factors file ", & 
!		            trim(adjustl(shape_file)) 
!
        Write(88,*) "A pre-built file is possibly available for ice SHs: ", & 
                     trim(adjustl(shape_file)) 
!
	Endif     
ENDIF
!
!
!
!
!
! ###### EXTERNAl elastic Love numbers ######  
!
! >>>>>>  A priori Love numbers
!
IF(line(1:3)=="035") THEN 
!
	call scan_string (line, 2, ss, nout)
!
	option_lnap   = ss(1)
!
  	lnap_filename = ss(2) 
!	
        if(option_lnap=='y') THEN   !  External elastic Love numbers are chosen <<<<<<<<<<<
!
	lnap_file  = trim(adjustl(wdir))//'/LOVE_DEPOSIT/'//trim(adjustl(lnap_filename))
!
        INQUIRE(FILE=lnap_file,EXIST=lex)
!
	CODE='-2'
!
        If(lex) then 
           Write(88,*) "SELEN will use a priori LNs file: ", trim(adjustl(lnap_filename))
	   Write(2,*)  "cp ", trim(adjustl(lnap_file)), " ", "./"//trim(adjustl(lnap_filename))
        Else 
           Write(88,*) "File ", trim(adjustl(lnap_filename)), " has NOT been found"
           Write(88,*) "The program will STOP ------------------------"
           Write(* ,*) "File ", trim(adjustl(lnap_filename)), " has NOT been found"
           Write(*, *) "The program will STOP ------------------------"
           call Stop_Config
        Endif 
!
        IF (MODE/='2') THEN 
           Write(88,*) "From config.dat, SELEN will be solved in mode:", mode  
	   Write(88,*) "W A R N I N G --- you have chosen external elastic Love numbers --- " 
	   Write(88,*) "and I suspect you DO NOT want an elastic analysis (for which mode=2)" 
           Write(88,*) "The program will STOP ----------------------------------------------"
           Write(*, *) "From config.dat, SELEN will be solved in mode:", mode  
	   Write(*, *) "W A R N I N G --- you have chosen external elastic Love numbers --- " 
	   Write(*, *) "and I suspect you DO NOT want an elastic analysis (for which mode=2)" 
           Write(*, *) "The program will STOP ----------------------------------------------"
           call Stop_Config
        ENDIF 
!
!--- Reads the ELASTIC LNs from the user-supplied file & prepares an AD-HOC set of SLE input files ...  
!    These are the (elastic) Green functions for:
! 
!         - sea level, 
!         - vertical displacement, 
!         - geoid elevation, and 
!         - horizontal displacement 
!         
!    Their viscoelastic components are forced to ZERO. 
!
! Modified GS February 2012 for the implementation of the multi-step elastic rebound
	ZEROV(:)=ZERO
!
        open(13,file= 'ebs.dat', status='unknown') 
        open(14,file= 'ebu.dat', status='unknown') 
        open(15,file= 'ebn.dat', status='unknown') 
        open(16,file= 'ebv.dat', status='unknown') 
	open(17,file='ebfa.dat', status='unknown')    
	open(18,file='ebss.dat', status='unknown')    
!		
!    Reading the header of the external file 
	open(12,file=lnap_file,status='unknown') 
        do j=1, 25    
		read(12,'(a10)') AJUNK  
	enddo 
        do j=0, LN_MAX
                read (12,*,end=98813) II, LN_H, LN_L, LN_K 
		den=2*ii+1
!
! Modified GS February 2012 for the implementation of the multi-step elastic rebound 
!
                NNINC=50
! GS April 2014 
		NNINC=202
!		
                write(13,*) II, (1. + ln_k - ln_h)/float(DEN), (ZEROV(K),K=1,NNINC+1)
                write(14,*) II, (          + ln_h)/float(DEN), (ZEROV(K),K=1,NNINC+1)
                write(15,*) II, (1. + ln_k       )/float(DEN), (ZEROV(K),K=1,NNINC+1) 
                write(16,*) II, (          + ln_l)/float(DEN), (ZEROV(K),K=1,NNINC+1)		
                write(17,*) II, ZERO,                          (ZEROV(K),K=1,NNINC+1)
                write(18,*) II, ZERO,                          (ZEROV(K),K=1,NNINC+1)				
	enddo 
98813   close(12) 
        close(13) ;  close(14) ;  close(15) ; close(16) ; CLOSE(17) ; CLOSE(18)
!
ENDIF ! <<<<<<<<<<<<<<<<<<<<<<<<<<<< External elastic Love numbers are chosen <<<<<<<<<<<
!
ENDIF
!
!
!
! IN PROGRESS === IN PROGRESS === IN PROGRESS === IN PROGRESS === IN PROGRESS ===
! IN PROGRESS === IN PROGRESS === IN PROGRESS === IN PROGRESS === IN PROGRESS ===
! IN PROGRESS === IN PROGRESS === IN PROGRESS === IN PROGRESS === IN PROGRESS ===
!
!
IF(line(1:3)=="036") THEN 
!
	option_tidal_ext  = 'n'
!	
        option_tlove_ext  = 'n'
!
	call scan_string (line, 1, ss, nout)
!
	option_tidal_ext  = ss(1)
!
        if(option_tidal_ext=='y') then     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!
	option_tlove_ext = 'y'
!
        IF (option_lnap/='y') THEN 
	   Write(88,*) "You have requested tidal Love numbers, but...  " 
	   Write(88,*) "... the external switch appears to be OFF.  ==== The program will STOP ==== "
	   Write(*, *) "You have requested tidal Love numbers, but...  " 
	   Write(*, *) "... the external switch appears to be OFF.  ==== The program will STOP ==== "
           call Stop_Config
        ENDIF 
!
	endif
!
        IF (option_tlove_ext=='y') THEN 
!
	   Write(88,*) "  ========== External tidal Love numbers are requested =========   " 
	   Write(88,*) " Just one Note: They are directly provided in progrem S_ROT_PM.F90 " 
	   Write(88,*) "  and NOT read from external files (as it is for loading numbers)  " 
!
!	   Write(* ,*) "  *** External tidal Love numbers are requested ***   " 
!	   Write(* ,*) " They are directly provided in program S_ROT_PM.F90   " 
!	   Write(* ,*) " and NOT read from external files (as it is for LLNs) " 
!
        ENDIF 
!


ENDIF

!
! IN PROGRESS === IN PROGRESS === IN PROGRESS === IN PROGRESS === IN PROGRESS ===
! IN PROGRESS === IN PROGRESS === IN PROGRESS === IN PROGRESS === IN PROGRESS ===
! IN PROGRESS === IN PROGRESS === IN PROGRESS === IN PROGRESS === IN PROGRESS ===
!
!


!
!
!
!
! ###### Earth model parameters ######  
!
! >>>>>>  Number of mantle Maxwell VE layers and code 
!         Mantle viscosity profile & lithospheric thickness
!
        IF(line(1:3)=="040") THEN 
!
	option_nm  = 'n'
!
	call scan_string (line, 4, ss, nout)
	option_nm  = ss(1)
!
        if(option_nm=='y') then     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!
	nv         = ss(2)
	code       = ss(3) 
	visco_file = './VSC/'//ss(4)	
!
        IF (option_lnap=='y') THEN 
	   Write(88,*) "W A R N I N G --- you have chosen external elastic Love numbers --- " 
	   Write(88,*) "and at the same time an analysis by TABOO. The program will STOP -- "
	   Write(*, *) "W A R N I N G --- you have chosen external elastic Love numbers --- " 
	   Write(*, *) "and at the same time an analysis by TABOO. The program will STOP -- "
           call Stop_Config
        ENDIF 
!
        Write(88,*) "Love numbers & spectra will be computed by Normal Modes (TABOO)"		
	Write(88,*) "... number of mantle v/e layers (NV): ", NV 
        Write(88,*) "... Earth model CODE (see TABOO User guide): ", CODE
!
!	 
! ---   Number of lines in the viscosity file 
!
	IF(NV/='0') THEN 
!
        INQUIRE(FILE=visco_file,EXIST=lex)
!
        If(lex) then 
	Write(88,*) "The viscosity file is: ", trim(adjustl(visco_file))
	        Else 
	Write(88,*) "The file ", trim(adjustl(visco_file)), " has NOT been found in ./VSC"
        Write(* ,*) "The file ", trim(adjustl(visco_file)), " has NOT been found in ./VSC"
        call Stop_Config
        	Endif
!
	open(12,file=visco_file,status='unknown') 
!
	LEN_VISCO=len(trim(adjustl(visco_file))) 
	visco_filename="'"//trim(adjustl(visco_file))//"'"
	NLINES=0
		do j=1, large_integer
		read(12,*,end=88099) junk
!		write(*,'(a80)')JUNK
		NLINES=NLINES+1
		enddo
	88099   close(12) 	
	call CHAR3_2_INT(nv, nnv)		 
	if(nnv==NLINES-1) then 
		Write(88,*) "The number of lines in file ", trim(adjustl(visco_file)), & 
			   " is consistent with NV=", NNV
		           else
		Write(88,*) "The number of lines in file ", trim(adjustl(visco_file)), & 
		 	    "is NOT consistent with NV=", NNV
		Write(88,*) "Please check config.dat and ", trim(adjustl(visco_file))
		Write(*,*)  "The number of lines in file ", trim(adjustl(visco_file)), & 
		 	    " is NOT consistent with NV=", NNV
		Write(*,*)  "Please check config.dat and ", trim(adjustl(visco_file))			  
	       call Stop_config 			    
	endif
!	
! ---   Reading litho thickness and viscosity values from "visco_file"
!
	open(12,file=visco_file,status='unknown') 
	read(12,'(a200)') row 
	call scan_string (row, 1, ss, nout)
	litho=ss(1) 
!	
		do j=1, nnv 
			read(12,'(a200)') row
			call scan_string (row, 1, ss, nout)
			vsc(j)=ss(1)			
		enddo
	close(12)
!
! ---   Viscosity string for output on monitor & GMT plots
!
	open(44,file='visco.tmp',status='unknown') 
	if(nnv==1) then 
		   write(44,*) "/",trim(vsc(1)),"/" 
        	   else 
		   Write(44,*) "/",trim(vsc(1))," ",& 
		   		  (trim(vsc(j))," ",j=2,nnv-1),&
		   		   trim(vsc(nnv)),"/"		       	       
	endif
	close(44) 
!	
	open(44,file='visco.tmp',status='unknown') 
	read(44,'(a100)') vstring 
	close(44) 	     	
!
	ENDIF	
!
	ENDIF                !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!
ENDIF
!
!
!
!
!
!
IF(line(1:3)=="045") THEN 
!
	option_tidal  = 'n'
!	
        option_tlove  = 'n'
!
	call scan_string (line, 1, ss, nout)
	option_tidal  = ss(1)
!
        if(option_tidal=='y') then     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!
	option_tlove  = 'y'
!
        IF (option_nm/='y') THEN 
	   Write(88,*) "You have requested a tidal Love number analysis by TABOO but the " 
	   Write(88,*) "TABOO switch appears to be OFF.  ==== The program will STOP ==== "
	   Write(*,*) "You have requested a tidal Love number analysis by TABOO but the " 
	   Write(*,*) "TABOO switch appears to be OFF.  ==== The program will STOP ==== "
           call Stop_Config
        ENDIF 
!
	endif
!
ENDIF
!
!
!
!
!
IF(line(1:3)=="055") THEN 
!
	option_pw  = 'n'
!
	call scan_string (line, 2, ss, nout)
!
	option_pw        = ss(1)
	alma_file_mode_2 = trim(adjustl(wdir))//'/VSC/'//ss(2)
        short_visco_filename='./VSC/'//ss(2)
!
!
!if(option_nm=='n'.and.option_pw=='n') then 
!			write(88,*) 'Love numbers are required! Please choose one of the'
!			write(88,*) 'two computation methods (Normal modes or ALMA) from'			
!			write(88,*) 'file <<config.dat>>'			
!			write(* ,*) 'Love numbers are required! Please choose one of the'
!			write(* ,*) 'two computation methods (Normal modes or ALMA) from'			
!			write(* ,*) 'file <<config.dat>>'			
!	        	call stop_config 
!	Endif	
	if(option_nm=='y'.and.option_pw=='y') then 
			write(88,*) 'You are requiring to compute Love numbers by ALMA,'
			write(88,*) 'but it appears that the Normal Modes switch is ON '			
			write(88,*) '  ****** Please check file <<config.dat>> ******  '			
			write(*, *) 'You are requiring to compute Love numbers by ALMA,'
			write(*, *) 'but it appears that the Normal Modes switch is ON '			
			write(*, *) '  ****** Please check file <<config.dat>> ******  '			
	        	call stop_config 	
	Endif 
!
	If(option_pw=='y') then  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!
	Write(88,*) "The Earth model file is ", trim(adjustl(alma_file_mode_2))
!
        INQUIRE(FILE=alma_file_mode_2,EXIST=lex)
!	
        If(lex) then 
!
		Write(88,*) "Love numbers are computed by Post-Widder formula (ALMA)"	
!	
		write(88,*) 'File ', trim(adjustl(alma_file_mode_2)), ' effectively exists...'		
		LEN_VISCO=len(trim(adjustl(alma_file_mode_2))) 
!
		visco_filename="'"//trim(adjustl(alma_file_mode_2))//"'"		
!
		open(101,file=trim(adjustl(alma_file_mode_2)),status='unknown') 
!
			do k=1, 4 
				read(101,'(a10)')cjunk 
				if(cjunk(1:1).ne.'!') then 
					write(88,*) 'Possible format error in viscosity file:'			
					write(88,*) trim(adjustl(alma_file_mode_2))
					write(88,*) 'Please verify that the format is compatible with ALMA '
					write(*, *) 'Possible format error in viscosity file:'	 		
					write(*, *) trim(adjustl(alma_file_mode_2))
					write(*, *) 'Please verify that the format is compatible with ALMA '			
	        			call stop_config 		
				endif 
			enddo
			NLINES=0 
				do k=1, VERY_LARGE_INTEGER
					read(101,'(a10)',end=98210) cjunk 
!					write(*,'(a10)')cjunk 	
					NLINES=NLINES+1 
			enddo
			
98210 	close(101) 
!
		If(NLINES.LE.0)then 
			write(88,*) 'Possible format error in viscosity file:'			
			write(88,*) trim(adjustl(alma_file_mode_2))
			write(88,*) 'Please verify that the format is compatible with ALMA '
			write(*, *) 'Possible format error in viscosity file:'	 		
			write(*, *) trim(adjustl(alma_file_mode_2))
			write(*, *) 'Please verify that the format is compatible with ALMA '			
	        	call stop_config 
		endif

!
	write(88,*) 'Model ', trim(adjustl(alma_file_mode_2)), ' has ', NLINES-2, ' mantle layers' 			
	CALL INT_2_CHAR3(NLINES-2,NV)
	CODE='-1'
!
		else		
!
		write(88,*) 'File ', trim(adjustl(alma_file_mode_2)), ' does not exists'
		write(*, *) 'File ', trim(adjustl(alma_file_mode_2)), ' does not exists'
	        call stop_config 
!
	endif	
!
 	CALL ALMA_2_CONFIG(wdir,visco_filename, nlines-2, degree, option_deg1)	
!		
ENDIF   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!
ENDIF
!	
!
!
!
!
!
IF(line(1:3)=="056") THEN 
!
	option_pwa  = 'n'
!
	call scan_string (line, 1, ss, nout)
!
	option_pwa  = ss(1)
!
        If(option_pwa=='y') Write(88,*) "Pre-computed ALMA Love numbers are requested..."	
!
	if(option_nm=='y'.and.option_pwa=='y') then 
			write(88,*) 'You are requiring to use pre-computed ALMA Love numbers but it appears'
			write(88,*) 'that the TABOO switch is ON! *** Please check file <<config.dat>> ****' 		
			write(* ,*) 'You are requiring to use pre-computed ALMA Love numbers but it appears'
			write(* ,*) 'that the TABOO switch is ON! *** Please check file <<config.dat>> ****' 		
	        	call stop_config 	
	Endif 
!
	if(option_nm=='n'.and.option_pw=='n'.and.option_pwa=='n'.and.option_lnap=='n') then 
			write(88,*) 'Love numbers ARE required! Please choose one of the two computation'
			write(88,*) 'methods (Normal modes or ALMA) from file <<config.dat>> or select a'			
			write(88,*) 'preexisting set of ALMA-computed Love numbers... '			
			write(*, *) 'Love numbers ARE required! Please choose one of the two computation'
			write(*, *) 'methods (Normal modes or ALMA) from file <<config.dat>> or select a'			
			write(*, *) 'preexisting set of ALMA-computed Love numbers... '			
	        	call stop_config 
	Endif	
	if(option_pw=='y'.and.option_pwa=='y') then 
			write(88,*) 'You are requiring to use pre-computed ALMA Love numbers but it appears'
			write(88,*) 'that the ALMA switch is ON! **** Please check file <<config.dat>> ****' 		
			write(* ,*) 'You are requiring to use pre-computed ALMA Love numbers but it appears'
			write(* ,*) 'that the ALMA switch is ON! **** Please check file <<config.dat>> ****' 		
	        	call stop_config 	
	Endif 
!
	If(option_pwa=='y') then    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!	
	alma_file_pwa_h  = trim(adjustl(wdir))//'/LOVE_DEPOSIT/hpw_'//trim(adjustl(short_visco_filename(12:100)))
	alma_file_pwa_l  = trim(adjustl(wdir))//'/LOVE_DEPOSIT/lpw_'//trim(adjustl(short_visco_filename(12:100)))
	alma_file_pwa_k  = trim(adjustl(wdir))//'/LOVE_DEPOSIT/kpw_'//trim(adjustl(short_visco_filename(12:100)))
	alma_logfile_pwa = trim(adjustl(wdir))//'/LOVE_DEPOSIT/alma-logfile_'//trim(adjustl(short_visco_filename(12:100)))
!
        INQUIRE(FILE = alma_file_mode_2,EXIST = lex)
        INQUIRE(FILE = alma_file_pwa_h,  EXIST = lexx(1))
        INQUIRE(FILE = alma_file_pwa_l,  EXIST = lexx(2))
        INQUIRE(FILE = alma_file_pwa_k,  EXIST = lexx(3))
        INQUIRE(FILE = alma_logfile_pwa, EXIST = lexx(4))
!	
	If(lex.and.lexx(1).and.lexx(2).and.lexx(3).and.lexx(4)) then    ! ++++++++++++++++++++++++++++++++++++++++++++++
!
		write(88,*) 'File ', trim(adjustl(alma_file_mode_2)), ' exists'		
		write(88,*) 'Browsing the directory of ALMA pre-computed Love numbers...'
		write(88,*) 'Pre-computd ALMA Love numbers refer to model: ', trim(adjustl(alma_file_mode_2)) 		
		write(88,*) 'Love number file ', trim(adjustl(alma_file_pwa_h)), ' exists'		
		write(88,*) 'Love number file ', trim(adjustl(alma_file_pwa_l)), ' exists'		
		write(88,*) 'Love number file ', trim(adjustl(alma_file_pwa_k)), ' exists'		
		write(88,*) 'ALMA log file ',  trim(adjustl(alma_logfile_pwa )), ' exists'		
!
		LEN_VISCO=len(trim(adjustl(alma_file_mode_2))) 
!
		visco_filename="'"//trim(adjustl(alma_file_mode_2))//"'"		
!
		open(101,file=trim(adjustl(alma_file_mode_2)),status='unknown') 
!
			do k=1, 4 
				read(101,'(a10)')cjunk 
				if(cjunk(1:1).ne.'!') then 
					write(88,*) 'Possible format error in viscosity file:'			
					write(88,*) trim(adjustl(alma_file_mode_2))
					write(88,*) 'Please verify that the format is compatible with ALMA '
					write(*, *) 'Possible format error in viscosity file:'	 		
					write(*, *) trim(adjustl(alma_file_mode_2))
					write(*, *) 'Please verify that the format is compatible with ALMA '			
	        			call stop_config 		
				endif 
			enddo
			NLINES=0 
				do k=1, VERY_LARGE_INTEGER
					read(101,'(a10)',end=98217) cjunk 
					NLINES=NLINES+1 
			enddo
!			
98217 	close(101) 
!
		If(NLINES.LE.0)then 
			write(88,*) 'Possible format error in viscosity file:'			
			write(88,*) trim(adjustl(alma_file_mode_2))
			write(88,*) 'Please verify that the format is compatible with ALMA '
			write(*, *) 'Possible format error in viscosity file:'	 		
			write(*, *) trim(adjustl(alma_file_mode_2))
			write(*, *) 'Please verify that the format is compatible with ALMA '			
	        	call stop_config 
		endif

!
	write(88,*) 'Model ', trim(adjustl(alma_file_mode_2)), ' has ', NLINES-2, ' mantle layers' 			
	CALL INT_2_CHAR3(NLINES-2,NV)
	CODE='-1'	
!
	Else 
!
        write(88,*) 'At least one of the following file apparently does not exist:'
        write(88,*) trim(adjustl(alma_file_mode_2))
        write(88,*) trim(adjustl(alma_file_pwa_h)) 	
        write(88,*) trim(adjustl(alma_file_pwa_l)) 	
        write(88,*) trim(adjustl(alma_file_pwa_k)) 	
        write(* ,*) 'At least one of the following file apparently does not exist:'
        write(*, *) trim(adjustl(alma_file_mode_2))
        write(*, *) trim(adjustl(alma_file_pwa_h)) 	
        write(*, *) trim(adjustl(alma_file_pwa_l)) 	
        write(*, *) trim(adjustl(alma_file_pwa_k)) 	
!
        call stop_config 
!	
	Endif 	
!	
	Endif   ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  
!
!
!
ENDIF
!
!
!
!
!
!
 IF(line(1:3)=="007") THEN 
!
 option_pmtf  = 'n'
!
 call scan_string (line, 1, ss, nout)
!
 option_pmtf  = ss(1)
!
 if(option_pmtf=='y'.and.option_pw=='y') THEN  
 write(88,*) 'The Polar Motion Transfer Function analysis is only possible when  Love numbers are '
 write(88,*) 'computed by normal modes (TABOO) ********* Please check file <<config.dat>> ********'	       
 write(88,*) 'The program will STOP'	       
 write(*, *) 'The Polar Motion Transfer Function analysis is only possible when  Love numbers are '
 write(*, *) 'computed by normal modes (TABOO) ********* Please check file <<config.dat>> ********'	       
 write(*, *) 'The program will STOP'	       
 CALL STOP_CONFIG      
 Endif 
!
 if((option_pmtf=='y'.and.option_tlove=='n').and.(option_tlove_ext=='n')) THEN 
 write(88,*) 'The Polar Motion Transfer Function analysis is only possible if the tidal analysis'
 write(88,*) 'by TABOO or external Love numbers has been scheduled  ********* Please check file <<config.dat>> **********'	       
 write(88,*) 'The program will STOP'	       
 write(*, *) 'The Polar Motion Transfer Function analysis is only possible if the tidal analysis'
 write(*, *) 'by TABOO or external Love numbers has been scheduled  ********* Please check file <<config.dat>> **********' 
 write(*, *) 'The program will STOP'	       
 CALL STOP_CONFIG      
 Endif 
!
 If(option_pmtf=='y'.and.option_tlove_ext=='y') THEN 
 write(88,*) 'Note: For external elastic Love numbers, the Polar Motion Transfer Function is  ' 
 write(88,*) 'not computed by a special program. Rather, it is evaluated by S_ROT_PM.F90 using' 
 write(88,*) 'the numerical values of the tidal Love numbers for the PREM model (by P.Gegout) '  
 ENDIF
!
 ENDIF
!
!
!
 IF(line(1:3)=="008") THEN 
!
 option_pmd  = 'n'
!
 call scan_string (line, 1, ss, nout)
!
 option_pmd  = ss(1)
! 
 if(option_pmd=='y'.and.option_pmtf=='n') then  
 write(88,*) 'The direct effect on polar motion REQUIRES the polar motion transfer function'
 write(88,*) 'Please check file <<config.dat>> ========== The program will STOP ========== '	      
 write(*, *) 'The direct effect on polar motion REQUIRES the polar motion transfer function'
 write(*, *) 'Please check file <<config.dat>> ========== The program will STOP ========== '	      
 CALL STOP_CONFIG      
 Endif
!
 if(option_pmd=='y'.and.option_tlove_ext=='y') then
 write(88,*) 'Note: For external elastic Love numbers, output files containing predictions for'
 write(88,*) 'the direct polar motion are NOT yet available ======== Work in progress ========'
 endif
!
ENDIF
!
!
!
 IF(line(1:3)=="009") THEN 
!
 option_rfb  = 'n'
!
 call scan_string (line, 1, ss, nout)
!
 option_rfb  = ss(1)
! 
 if(option_rfb=='y'.and.option_pmd=='n') then  
 write(88,*) 'Modelling the polar motion feedback REQUIRES evaluation of the polar motion direct effect'
 write(88,*) 'Please check file <<config.dat>> +++++++========= The program will STOP ==========+++++++'	      
 write(*, *) 'Modelling the polar motion feedback REQUIRES evaluation of the polar motion direct effect'
 write(*, *) 'Please check file <<config.dat>> +++++++========= The program will STOP ==========+++++++'      
 CALL STOP_CONFIG      
 Endif
!
 if(option_rfb=='y'.and.option_tlove_ext=='y') then  
 write(88,*) 'Note: For external elastic Love numbers, the rotational feedback is implemented  ' 
 write(88,*) 'by program S_ROT_PM.F90, which also builds the Polar Motion Transfer Function ...' 
 Endif

!
ENDIF
!
!
!
!
! ###### Tegmark resolution ######
!
IF(line(1:3)=="015") THEN 
	call scan_string (line, 1, ss, nout)
        resolution=ss(1)
!
 	call CHAR10_2_INT(resolution, nresolution)
!
 	np=2*nresolution*(nresolution-1)*20+12
!
 	Write(88,*) "Tegmark's resolution, r: ", nresolution 
 	Write(88,*) "Number of pixels: ", np
!
!
! ---- A lower bound for Resolution 
 	If(nresolution.lt.res_min) then 
 	Write(*, *) "------ Resolution it TOO LOW -------"
 	Write(*, *) "  Please use a resolution RES>=12   "
 	Write(88,*) "------ Resolution it TOO LOW -------"
 	Write(88,*) "  Please use a resolution RES>=12   "
        Write(88,*) 
        call Stop_Config
	stop 
 	endif 

!
! ---- Testing Tegmark's condition 
 	If(np>=ndegree**2/3.and.np/=0) then 
 	Write(88,*) "------ 'Tegmark condition' is met ------"
 	Write(88,*) "(enough pixels for this harmonic degree)"
				       else 
 	Write(88,*) " ------- Tegmark condition is NOT met ------- "
 	Write(88,*) " (not enough pixels for this harmonic degree) "
 	endif 
ENDIF 
!
!
! ###### Style of computation of derivatives ######

IF(line(1:3)=="018") THEN   
        option_der='9'
	call scan_string (line, 1, ss, nout)
	option_der = ss(1) 
!	
 If(option_der=='1') Write(88,*) "Present time derivatives are computed in one step"	
 If(option_der=='2') Write(88,*) "Present time derivatives are computed in two steps"	
 If(option_der/='1'.and.option_der/='2') then 
        Write(88,*) 'Derivatives can be of type 1 or 2'	
        Write(*, *) 'Derivatives can be of type 1 or 2'	
        Call Stop_Config ; stop	 	
 Endif   
ENDIF
!
!
!
!
!
! ###### Pixel table ######
!
IF(line(1:3)=="016") THEN 
	call scan_string (line, 2, ss, nout)
        option_npx   = ss(1) 
        file_pxtable = ss(2) 
!
	if(option_npx=='y') then 
		 Write(88,*) "SELEN will prepare the new pixel table file: ", file_pxtable
        else
	         INQUIRE(FILE=file_pxtable,EXIST=lex)			   
	         If(lex) then 
		     Write(88,*) "SELEN will use the pre-built pixel table file: ", trim(adjustl(file_pxtable))
		 Else 
                 Write(88,*) "The pixel table file ", trim(adjustl(file_pxtable)), " has NOT been found"
	             Write(* ,*) "The pixel table file ", trim(adjustl(file_pxtable)), " has NOT been found"
	             call Stop_Config
	         Endif		   
        Endif 
ENDIF
!
!
!
! ###### Spherical harmonics (SH) file ######
!
IF(line(1:3)=="070") THEN 
	call scan_string (line, 2, ss, nout)
        option_sh = ss(1) 
        sh_file   = ss(2) 
	if(option_sh=='y') then 
		              Write(88,*) "SELEN will prepare the new SH file: ", sh_file
		           else
			      INQUIRE(FILE=sh_file,EXIST=lex)
			      If(lex) then 
			           Write(88,*) "SELEN will use the pre-built SH file: ", trim(adjustl(sh_file))
			      Else 
			           Write(88,*) "The SH file ", trim(adjustl(sh_file)), " has NOT been found"
			           Write(* ,*) "The SH file ", trim(adjustl(sh_file)), " has NOT been found"
			           call Stop_Config
			      Endif		   
       Endif 
ENDIF
!
!
!
! ###### Kind of OCEAN FUNCTION (realistic vs. zonal) ######
!
IF(line(1:3)=="075") THEN 
!
	call scan_string (line, 2, ss, nout)
!
        option_rof = ss(1) 
	radius_zof = ss(2)
!
        if(option_rof=='z'.and.(ice_type=='po'.or.ice_type=='pm'))then 
          		Write(88,*) "A zonal OF is still unavailable for present time ice sheets" 
		        Write(*, *) "A zonal OF is still unavailable for present time ice sheets"
	  		call Stop_Config	
        Endif
	If(option_topo=='n') then 
		if(option_rof=='r') Write(88,*) "Realistic Ocean Function"
		if(option_rof=='z') Write(88,*) "ZONAL Ocean Function with radius (degrees): ", & 
	                    trim(adjustl(radius_zof))
		if(option_rof/='r'.and.option_rof/='z')then 
          		Write(88,*) "Type 'r' for a realistic OF, 'z' for a zonal OF in 'config.dat'" 
          		Write(*, *) "Type 'r' for a realistic OF, 'z' for a zonal OF in 'config.dat'" 
	  		call Stop_Config	
		Endif
	Else 
			if(option_rof=='z') then 
			Write(88,*) "ZONAL Ocean Function conflicts with evolving coastlines"
          		Write(*, *) "ZONAL Ocean Function conflicts with evolving coastlines" 
	  		call Stop_Config
			endif 		
	Endif  
				    
ENDIF 
!
!
!
!
!
! ###### OF SH decomposition ######
!
IF(line(1:3)=="080") THEN 
	call scan_string (line, 2, ss, nout)
        option_oh = ss(1) 
        shof_file = ss(2) 
	if(option_oh=='y') then 
		 Write(88,*) "SELEN will prepare the new SH OF file: ", shof_file
		           else
			      INQUIRE(FILE=shof_file,EXIST=lex)			   
			      If(lex) then 
			           Write(88,*) "SELEN will use the pre-built SH OF file: ", trim(adjustl(shof_file))
			      Else 
			           Write(88,*) "The SH OF file ", trim(adjustl(shof_file)), " has NOT been found"
			           Write(* ,*) "The SH OF file ", trim(adjustl(shof_file)), " has NOT been found"
			           call Stop_Config
			      Endif		   
        Endif 
ENDIF
!
!
!
!
! ###### Reading the depot Label ######  		
!
IF(line(1:3)=="090") THEN
	call scan_string (line, 1, ss, nout)
		run=ss(1)
		write(88,*) 'The depot label is: ', run	
		depot="./depot-"//trim(adjustl(run))
		runp="'"//trim(adjustl(run))//"'"
ENDIF  
!
!
100 CONTINUE
!
444 close(1)
!
!
!
! # Part 2/1 Reading the output settings from "config.dat"   
!
!
  open(1,file='config.dat',status='unknown')
!
!
  Write(88,*) ""
  Write(88,*) "+----------------------------------------------+"
  Write(88,*) "|               Outputs of SELEN               |"
  Write(88,*) "+----------------------------------------------+"
!
!
  DO 200 I=1,LARGE_INTEGER
!
!
  read(1,'(a200)',end=555) line 
!
!
!
!
! ###### How to deal with GMT scripts ###### 
!
IF(line(1:3)=="666") THEN 
	call scan_string (line, 1, ss, nout)
	option_GMT = ss(1) 
	If(option_gmt=='n') & 
	Write(88,*) "Execution of GMT scripts is switched off"
	If(option_gmt/='y'.and.option_gmt/='n') then 
		Write(* ,*) "For the GMT switch, only y/n are valid options"
		Write(88,*) "For the GMT switch, only y/n are valid options"
		Call stop_config 
		Stop	
	Endif
ENDIF
!
!
!
!
! ###### Pixelization ###### 
!
! >>>>>> Producing pixelization maps
!
IF(line(1:3)=="100") THEN 
	call scan_string (line, 1, ss, nout)
	option_px = ss(1) 
	If(option_px=='y') & 
	Write(88,*) "Maps of the pixelizations will be drawn"
ENDIF
!
! >>>>>> Computing and drawing the function 
!
IF(line(1:3)=="110") THEN 
	call scan_string (line, 1, ss, nout)
	option_wi = ss(1) 
	If(option_wi=='y') &
	Write(88,*) "The window function will be computed & plotted"
ENDIF
!
!
!
!
! ###### Ocean function reconstruction and mapping ###### 
!
IF(line(1:3)=="120") THEN 
	call scan_string (line, 1, ss, nout)
	option_of = ss(1) 
	If(option_of=='y') & 
	Write(88,*) "The OF will be reconstructed & mapped"
ENDIF
!
! ###### Plot of Ocean Function "degree variance" ###### 
!
IF(line(1:3)=="125") THEN 
	call scan_string (line, 1, ss, nout)
	option_ofdv = ss(1) 
	If(option_ofdv=='y') & 
	Write(88,*) "The OF DV will be computed & plotted"
ENDIF
!
!
! ###### Paleo-topography maps ###### 
!
IF(line(1:3)=="126") THEN 
	call scan_string (line, 1, ss, nout)
	option_ptmap = ss(1) 
!
	If(option_ptmap=='y') then
	If(option_topo=='y')  then 
		Write(88,*) "Paleotopography maps will be prepared"
	else 
		Write(88,*) "WARNING/ Paleo-topography maps can only be obtained for varying coastlines!"
		option_ptmap='n'	
	endif	
!
	endif
!
ENDIF
!
!
!
! ###### Ice sheets graphical outputs ###### 
!
! >>>>>> Maps of original ice sheets  
!
IF(line(1:3)=="130") THEN 
	call scan_string (line, 1, ss, nout)
	option_or = ss(1)
	If(option_or=='y') & 
	Write(88,*) "Maps of original ice sheets will be drawn"
ENDIF
!
! >>>>>> Plot of Equivalent Sea Level (ESL) vs time for "ho" ice models  
!
IF(line(1:3)=="135") THEN 
	call scan_string (line, 1, ss, nout)
	option_esl = ss(1)
!	
	if(option_esl=='y'.and.ice_type.ne.'ho')THEN
		Write(88,*) "This ESL option is only available for 'ho' ice models"
		Write(*, *) "This ESL option is only available for 'ho' ice models"
        CALL STOP_CONFIG 
	ENDIF
!
	If(option_esl=='y') & 
	Write(88,*) "A plot of Equivalent Sea level ('ho' ice model) will be drawn"
ENDIF
!
IF(line(1:3)=="137") THEN 
	call scan_string (line, 1, ss, nout)
	option_esl_pm = ss(1)
!	
	if(option_esl_pm.eq.'y'.and.ice_type.ne.'pm')THEN
		Write(88,*) "This ESL option is only available for 'pm' ice models"
		Write(*, *) "This ESL option is only available for 'pm' ice models"
        CALL STOP_CONFIG 
	ENDIF
!
	If(option_esl_pm=='y') & 
	Write(88,*) "A plot of Equivalent Sea level ('pm' ice model) will be drawn"
ENDIF
!
! >>>>>> Ice sheets reconstruction from SH coefficients & mapping
!
IF(line(1:3)=="140") THEN 
	call scan_string (line, 1, ss, nout)
	option_ri = ss(1)
	If(option_ri=='y') & 
	Write(88,*) "The ice sheets will be reconstructed from SH coefficients and mapped"
ENDIF
!
!
!
!
! ###### Earth model Spectral properties ###### 
!
! >>>>>> Diagrams of LDCs, relaxation spectrum & residues
!	
IF(line(1:3)=="150") THEN 
	call scan_string (line, 1, ss, nout)
	option_ln = ss(1) 
	If(option_ln=='y') & 
	Write(88,*) "The Love numbers will be plotted"
ENDIF
!
!
!
!
! ###### Plots of polar motion and rate of polar motion ###### 
!	
IF(line(1:3)=="310") THEN 
	call scan_string (line, 1, ss, nout)
	option_pmplot = ss(1) 
	If(option_pmplot=='y') Write(88,*) "PM and rate of PM will be plotted"
ENDIF
!
!
!
! ###### Relative Sea level (RSL) analysis ###### 
! [Only available for LGM to Holocene ice sheets]
!
! >>>>>>  RSL database
!
IF(line(1:3)=="160") THEN 
!
	call scan_string (line, 3, ss, nout)
!
	OPTION_RSLA = SS(1) 
!	
        if(ice_type=='po'.or.ice_type=='pm') option_rsla='n'
!
        If(option_rsla=='y') then 
!
		RSL_FILE = './DATA/'//SS(2)
!		
		RSL_DATABASE_FORMAT = SS(3)		
!
        	INQUIRE(FILE=RSL_FILE,EXIST=lex)
!
        	if(lex)then 
!
		len_rsl=len(trim(adjustl(rsl_file))) 
!
		rsl_database="'"//trim(adjustl(rsl_file))//"'"
!
! ---   Counts the number of RSL sites in database  
! 
! NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW DEC 2010
! NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW DEC 2010
!
	IF(RSL_DATABASE_FORMAT=='3')    then 
!
! ...   Counting the header lines 
!
        open(10,file=rsl_file,status='old') 
        N_HEADER_LINES=0 
        do k=1, VERY_LARGE_INTEGER
			read(10,'(a80)',end=23132)  ANOTHERROW
			if(ANOTHERROW(1:1)=='!') N_HEADER_LINES = N_HEADER_LINES + 1
		enddo
23132   close(10) 
!
! ...   Counting the number of RSL sites 
!
        open(10,file=rsl_file,status='old') 
	do k=1, N_HEADER_LINES
		read(10,'(a80)')  ANOTHERROW
	enddo
		nrsl=0
		do 33473 j=1, very_large_integer
		read(10,'(a200)',end=41491) row 
			do k=1, 200
			do n=1, NALFA
				if(row(k:k)==alfa(n)) then 
				nrsl=nrsl+1
				goto 33473
				endif 
			enddo
			enddo
		33473 continue
41491   close(10)
        nrsl=nrsl/2   
!
!
! NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW DEC 2010
! NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW DEC 2010
!
	ELSEIF(RSL_DATABASE_FORMAT=='0') then 
!
		open(10,file=rsl_file,status='old') 
!
		nrsl=0
		do 33433 j=1, very_large_integer
		read(10,'(a200)',end=41492) row 
			do k=1, 200
			do n=1, NALFA
				if(row(k:k)==alfa(n)) then 
				nrsl=nrsl+1
				goto 33433
				endif 
			enddo
			enddo
		33433 continue
41492	close(10) 
!
	Elseif(RSL_DATABASE_FORMAT=='1'.or.RSL_DATABASE_FORMAT=='2') then 
!	
		open(10,file=rsl_file,status='old') 
!
		nrsl=0
		do 33533 j=1, very_large_integer
		read(10,'(a200)',end=484) row 
!		
!               write(*,'(a100)') row 
!		
			do k=1, 200
			do n=1, NALFA
				if(row(k:k)==alfa(n)) then 
				nrsl=nrsl+1
				goto 33533
				endif 
			enddo
			enddo
			
		33533 continue
484	close(10)		
	If(RSL_DATABASE_FORMAT=='1')nrsl=nrsl/2 
	If(RSL_DATABASE_FORMAT=='2')nrsl=nrsl/3 
!			
	Endif
!
        Write(88,*) 'The RSL database is: ', trim(rsl_file), ' ---', nrsl, ' RSL sites'
!
	If(nrsl==0) then 
		Write(88,*) "The RSL database has no elements"
		Write(*, *) "The RSL database has no elements"
		Call stop_config 
		Stop	
	Endif
!
!
! --- File for scattered RSl data 
 		open(82,file='scatter-data.dat',status='unknown') 
!
! --- Reading the coordinates of the RSL sites
!
	        If    (RSL_DATABASE_FORMAT=='0') then 	
!
  		OPEN(44,FILE=RSL_FILE,STATUS='old')
!
     		ii=1 
299  		READ (44, 699, END=199) j, lats, lons, nrsl_data     
   		IF(lons<=0.) lons=360.+lons		    
   		do k=1, nrsl_data 
           		READ (44 ,*) time_rsl, rjunk, rsl_datum, d_rsldatum
	   	Write(82,*) time_rsl/1000., rsl_datum, d_rsldatum 
   		enddo
   		ii=ii+1
   		IF(ii <= nrsl) GOTO 299
199  CLOSE(44)
699  FORMAT(1X,I3,1X,F5.1,1X,F6.1,1X,I2,1X,A22)
!
!
!
! NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW DEC 2010
! NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW DEC 2010
!
		Elseif(RSL_DATABASE_FORMAT=='3') then
!
! Header lines ...
  		OPEN(44,FILE=RSL_FILE,STATUS='old')
		do k=1, N_HEADER_LINES
			read(44,'(a80)')  ANOTHERROW	
!write(*,'(a80)')  ANOTHERROW	
		enddo
! 
     		ii=1 
2991  		READ (44, *, END=1991) jj 
	        READ (44,'(a80)')      ANOTHERROW	  		       
	        READ (44,'(a80)')      ANOTHERROW	  		       
		READ (44, *)           LONS, LATS, NRSL_DATA   
		
!write(*,*) jj
!write (*,'(a80)')      ANOTHERROW	  		       
!write (*,'(a80)')      ANOTHERROW	  		       
!write (*, *)           LONS, LATS, NRSL_DATA   
		
		 
		IF(LONS<=0.) LONS=360.+LONS		    
!		
! Data ... 
!
   		do k=1, nrsl_data 
!
           		READ (44 ,*) j, jj, time_rsl, rjunk, rsl_datum, d_rsldatum
!
			if(j==1) time_rsl = time_today - time_rsl 
			if(j==2) time_rsl =              time_rsl 		        
!
	   		Write(82,*) time_rsl/1000., rsl_datum, d_rsldatum 
!                       Write(* ,*) time_rsl/1000., rsl_datum, d_rsldatum 
   		enddo
!
   		ii=ii+1
!
   		IF(ii <= nrsl) GOTO 2991
!
1991  CLOSE(44)
!
! NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW DEC 2010
! NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW DEC 2010
!
!
		Elseif(RSL_DATABASE_FORMAT=='1') then
!
  		OPEN(44,FILE=RSL_FILE,STATUS='old')
!
		Do ii=1, nrsl 
			READ(44,'(a10)')CJUNK 
			READ(44,'(a10)')AJUNK 
			READ(44,'(a10)')CJUNK 
			READ(44,*) LATS, LONS 
!									
			IF     (AJUNK(1:1)=="-") THEN 			
				READ(44,*) TIME_RSL, RJUNK, RSL_DATUM, D_RSLDATUM 
							
			ELSEIF (AJUNK(1:1)=="*") THEN 	
					
                                READ(44,*) TIME_RSL1, TIME_RSL2, RSL_DATUM, D_RSLDATUM 
				
			        TIME_RSL =  MIN(TIME_RSL1,TIME_RSL2)+ &
				           (MAX(TIME_RSL1,TIME_RSL2)-MIN(TIME_RSL1,TIME_RSL2))/2.
					   			
			ELSEIF (AJUNK(1:1)/="*".AND.AJUNK(1:1)/="-") THEN 
!						
				Write(88,*) "The RSL database is badly configured"
				Write(*, *) "The RSL database is badly configured"
				Call stop_config 
				Stop					
			ENDIF			 		
			READ(44,'(a10)')CJUNK 	
 			WRITE(82,*) TIME_RSL/1000., RSL_DATUM, D_RSLDATUM		
!			
		Enddo
		Close(44) 
!
		Elseif(RSL_DATABASE_FORMAT=='2') then	
! 
  		OPEN(44,FILE=RSL_FILE,STATUS='old')
!
		Do ii=1, nrsl 
			READ(44,'(a10)')CJUNK 
			READ(44,'(a10)')AJUNK 
			READ(44,'(a10)')CJUNK 
!
			READ(44,'(a200)') linep 
!			
			call scan_string (linep, 2, ss, nout)	
!			
			LATSC10=trim(adjustl(ss(1))) 
			LONSC10=trim(adjustl(ss(2)))
!							
			call FROMDDMMSS_2_DEGREES (LATSC10, LATS) 
			call FROMDDMMSS_2_DEGREES (LONSC10, LONS) 
!									
			IF     (AJUNK(1:1)=="-") THEN 			
				READ(44,*) TIME_RSL, RJUNK, RSL_DATUM, D_RSLDATUM 
							
			ELSEIF (AJUNK(1:1)=="*") THEN 	
					
                                READ(44,*) TIME_RSL1, TIME_RSL2, RSL_DATUM, D_RSLDATUM 
				
			        TIME_RSL =  MIN(TIME_RSL1,TIME_RSL2)+ &
				           (MAX(TIME_RSL1,TIME_RSL2)-MIN(TIME_RSL1,TIME_RSL2))/2.
					   			
			ELSEIF (AJUNK(1:1)/="*".AND.AJUNK(1:1)/="-") THEN 
			 
			
				Write(88,*) "The RSL database is badly configured"
				Write(*, *) "The RSL database is badly configured"
				Call stop_config 
				Stop					
			ENDIF			 		
			READ(44,'(a10)')CJUNK 	
 			WRITE(82,*) TIME_RSL/1000., RSL_DATUM, D_RSLDATUM		
!			
		Enddo
		Close(44) 
!
	ENDIF 
! 
     CLOSE(82) 	     
!
    	 Else
		Write(88,*) "The RSL database does not exists"
		Write(*, *) "The RSL database does not exists"
		Call stop_config 
		Stop	
!     
     	 Endif 	
!
     Else 
		Write(88,*) "NO RSL analysis will be performed "
     Endif     
!
ENDIF
!
!
!
!
IF(option_rsla=='y') THEN 
!
!
!
! >>>>>>  Plotting the distribution of RSL sites in database 
!
IF(line(1:3)=="170") THEN 
	call scan_string (line, 1, ss, nout)
	option_rsldb   = ss(1)
	If(option_rsldb=='y') Write(88,*) "The distribution of RSL stations will be plotted"
ENDIF
!
!
! >>>>>> Computing RSL predictions and individual plots of predictions vs. observations 
!
IF(line(1:3)=="180") THEN 
	call scan_string (line, 2, ss, nout)
	option_rsl   = ss(1) 
	option_rslp  = ss(2) 
	If(option_rsl =='y') Write(88,*) "RSL predictions will be obtained"
	If(option_rslp=='y') Write(88,*) "RSL predictions will be plotted against data"
	If(option_rslp=='y'.and.option_rsl=='n') then 
			Write(88,*) "Plotting RSl predictions vs. data requires computing ..." 
			Write(88,*) "... RSL predictions -- please modify file config.dat' --"
                        Write(* ,*) "Plotting RSl predictions vs. data requires computing ..." 
			Write(* ,*) "... RSL predictions -- please modify file config.dat' --"
			call stop_config 			    
			stop 
			endif
ENDIF
!
!
! >>>>>> Scatterplot of RSL observations vs. predictions 
!
IF(line(1:3)=="190") THEN 
	call scan_string (line, 1, ss, nout)
	option_rslsca   = ss(1) 
	If(option_rslsca=='y') Write(88,*) "A RSL scatterplot (global data vs predictions) will be drawn"
	If(option_rslsca=='y'.and.option_rsl=='n') then 
			Write(88,*) "Plotting a RSL 'scatterplot' showing data & predictions requires ..."
			Write(88,*) "... computing RSL predictions first. Please modify file 'config.dat'"
			Write(*,*)  "Plotting a RSL 'scatterplot' showing data & predictions requires ..."
			Write(*,*)  "... computing RSL predictions first. Please modify file 'config.dat'"
			call stop_config 			    			
			stop
			endif				
ENDIF
!
!
! >>>>>> Misfit analysis for RSL   
!
IF(line(1:3)=="200") THEN 
	call scan_string (line, 1, ss, nout)
	option_rslmf   = ss(1) 
	If(option_rslmf=='y') Write(88,*) "A misfit analysis for RSL will be performed"
	If(option_rslmf=='y'.and.option_rsl=='n') then 
			Write(88,*) "Plotting the RSL <<misfit histogram>> requires computing "
			Write(88,*) "... RSL predictions first. Please modify file config.dat'"
                        Write(*, *) "Plotting the RSL <<misfit histogram>> requires computing "
			Write(*, *) "... RSL predictions first. Please modify file config.dat'"
			Call Stop_Config 
			stop
			endif		
ENDIF
!
!
!
! >>>>>> RSL table with all observations & predictions   
!
IF(line(1:3)=="205") THEN 
	call scan_string (line, 1, ss, nout)
	option_rsltab   = ss(1) 
	If(option_rsltab=='y') Write(88,*) "A table with all data and predictions will be printed" 
	If(option_rsltab=='y'.and.option_rsl=='n') then 
			Write(88,*) "Printing a RSL table with all predictions requires computing ..."
			Write(88,*) "... RSL predictions first. ** Please modify file config.dat ** '"
			Write(* ,*) "Printing a RSL table with all predictions requires computing ..."
			Write(* ,*) "... RSL predictions first. ** Please modify file config.dat ** '"
			Call Stop_Config 
			stop
			endif	
ENDIF
!
ENDIF    ! on "option_rsla"
!
!
!
! >>>>>> RSL "regions"    
!
! [Only available for LGM to Holocene ice sheets]
!
! >>> Global RSL zones
!
IF(line(1:3)=="210") THEN 
	call scan_string (line, 1, ss, nout)
	option_rslz   = ss(1) 
!
	if(ice_type=='po'.or.ice_type=='pm') option_rslz='n'
!
	If(option_rslz=='y') Write(88,*) "<<RSL zones>> will be determined"
ENDIF
!
! >>> Regional RSL contour lines 
!
IF(line(1:3)=="215") THEN 
	call scan_string (line, 2, ss, nout)
	option_rslc     = ss(1) 
!
	file_region     = './DATA/'//ss(2) 
!
	if(option_rslc=='y') then 	
!
        INQUIRE(FILE=file_region,EXIST=lex)
!
	If(lex) then
!
                Write(88,*) "<<RSL contour lines>> will be drawn" 
		Write(88,*) "The RSL region is described in file ", file_region
!
! --- Counting the number of points ("virtual RSL sites") in rectangle
!
	        file_region_lonlat='lonlat_rslc.dat'
!
! --- Sbr. Scan_region scans the user-supplied input file "file_region" in order to:
!		- Check lon-lat bounds 
!		- Check the contour interval 
!		- Count the number of "Virtual RSL" sites
! 	        - Report lon-lat coordinates on file "file_region_lonlat"
!		---> NOTICE that the rectangle is enlarged by 2 
!		     degrees NSEW to avoid border effects with GMT 
!
  	        Call scan_region (file_region, file_region_lonlat, nrslc, time_bpc, time_bpcc, & 
				  LONMINC, LONMAXC, LATMINC, LATMAXC, & 
				  min_rslc, max_rslc, RSL_INT, & 
				  name_of_region)				 	
!
	        LEN_RSLC = len(trim(adjustl(file_region_lonlat)))
	        RSLC_LONLAT_FILE     = "'"//trim(adjustl(file_region_lonlat))//"'"	
!
                Write(88,*) "There are ", nrslc, " points in file ", trim(adjustl(file_region)) 
                Write(88,*) "Time of analysis is ", time_bpc, " ka" 
		Write(88,*) "Coordinates of RSL <<sites>> are reported on file ", & 
		             trim(adjustl(file_region_lonlat))
!		
   	        else 
                	Write(*, *) "File ", trim(adjustl(file_region)), " apparently"
			Write(*, *) "does not exist -  Please check file <config.dat>"
                	Write(88,*) "File ", trim(adjustl(file_region)), " apparently"
			Write(88,*) "does not exist -  Please check file <config.dat>"
			Call Stop_Config 
			stop		 
		endif 
	endif 
ENDIF
!
!
!
!
! ###### Sea level change at tide-gauge stations ###### 
!
!
! >>>>>> Tide-gauges database 
!
IF(line(1:3)=="220") THEN 
!
	call scan_string (line, 3, ss, nout)
!	
 	OPTION_TGA    = SS(1)
!
	TGAUGES_FILE  = './DATA/'//SS(2) 
!
        if(option_tga=='y') then 
!
	TGAUGES_DATABASE_FORMAT = SS(3)
	
! ---   Counts the number of tide-gauge sites in database  
	
        INQUIRE(FILE=TGAUGES_FILE,EXIST=lex)
!
	If(lex)then 
!
		Write(88,*) 'The tide gauges database is: ', TGAUGES_FILE
!
		LEN_TGAUGES=len(trim(adjustl(TGAUGES_FILE))) 
		TGAUGES_DATABASE="'"//trim(adjustl(TGAUGES_FILE))//"'"
!
		open(10,file=TGAUGES_FILE,status='old') 
		NTIDEGAUGES=0
		do 33833 j=1, very_large_integer
			read(10,'(a200)',end=464) row 
			if(row(1:1)/='#') NTIDEGAUGES=NTIDEGAUGES+1
		33833 continue
		464 continue 
        	Write(88,*) 'There are ', NTIDEGAUGES, 'tide-gauges records in file ', TGAUGES_FILE
 		close(10)
!
		else
!	
		Write(88,*) 'No tide-gauge database has been found...'	
		Write(*, *) 'No tide-gauge database has been found...'	
                Call Stop_Config 
                stop		 		
!			
		endif	
!
	endif
!
ENDIF
!
IF(option_tga=='y') THEN 
!
! >>>>>> Plotting tide-gauge stations distribution 
!
IF(line(1:3)=="230") THEN 
	call scan_string (line, 1, ss, nout)
	option_tgplot = ss(1) 
	If(option_tgplot=='y')&
	Write(88,*) 'The tide gauge sites distribution will be mapped'
ENDIF
!
! >>>>>> Tide gauge data scatterplot & averages
!
IF(line(1:3)=="240") THEN 
	call scan_string (line, 1, ss, nout)
	option_tgsca = ss(1) 
!
  	if(TGAUGES_DATABASE_FORMAT=='1'.and.option_tgsca=='y') then 
        Write(88,*) 'The scatterplot option is not available for tide gauge data files of format 1'
        Write(*, *) 'The scatterplot option is not available for tide gauge data files of format 1'
        Call Stop_Config ; stop 		
	Endif 
!
	If(option_tgsca=='y') Write(88,*) 'Tide gauge data scatterplots will be drawn'
ENDIF
!
!
! >>>>>> Predictions at tide-gauges 
!
IF(line(1:3)=="250") THEN 
	call scan_string (line, 1, ss, nout)
	option_tg = ss(1) 
	If(option_tg=='y')& 
	Write(88,*) 'Predictions of SLC at tide-gauges will be obtained'
ENDIF
!
ENDIF
!
!
!
!
! ###### Elastic rebound ###### (one-step melting ("po" ice type) 
!
! [Only available for present ice sheets]
!
! >>>>>> Global and Regional maps of S, U, and N  
!
IF(line(1:3)=="700") THEN 
!
	call scan_string (line, 1, ss, nout)
	option_reb = ss(1) 	
!
	if(ice_type=='ho')then 
		option_reb='n'
		goto 90009 
	endif
	If(option_reb=='y') Write(88,*) "Some elastic rebound maps will be drawn."	
	If(option_reb=='n') Write(88,*) "NO elastic rebound maps will be drawn."	
!
	If(mode/='2'.and.option_reb=='y') then  
		Write(88,*) "For elastic rebound, the SLE mode must be set to 2 (Elastic, GSC)"	
		Write(*, *) "For elastic rebound, the SLE mode must be set to 2 (Elastic, GSC)"	
		Call Stop_Config 
	Endif 
!
ENDIF
!
! IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - 
! IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - 
!
IF(line(1:3)=="702") THEN 
!
	option_xygrid = '1' 
!
	call scan_string (line, 1, ss, nout)
	option_xygrid = ss(1) 
!
	if(option_xygrid=='1') Write(88,*) "Elastic rebound is computed on the pixelized grid"
	if(option_xygrid=='2') Write(88,*) "Elastic rebound is computed on a 1 deg x 1 deg grid"
!
	If(option_xygrid/='1'.and.option_xygrid/='2') then  
		Write(88,*) "For elastic rebound, the grid option should be <<1>> or <<2>>"	
		Write(*, *) "For elastic rebound, the grid option should be <<1>> or <<2>>"	
		Call Stop_Config 
	Endif 

ENDIF
!
!
IF(line(1:3)=="703") THEN 
!
!	call scan_string (line, 2, ss, nout)
!	option_scalfa = ss(1) 
!	scalfa_string = trim(adjustl(ss(2))) 
!
!	If(option_scalfa/='y'.and.option_scalfa/='n') then  
!		Write(88,*) "The scaling factor option should be <<y>> or <<n>>"	
!		Write(*, *) "The scaling factor option should be <<y>> or <<n>>"	
!		Call Stop_Config 
!	Endif 
!
!	if(option_scalfa=='y') then 
!	Write(88,*) "A scaling factor is applied to ice thickness"
!	Write(88,*) "The scaling factor is: ", trim(adjustl(scalfa_string))
!        endif
!
!	if(option_scalfa=='n') scalfa_string='1.0' 
!
!        CALL CHAR100_2_REAL(scalfa_string, scalfa_num)
!
ENDIF
!
! IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - 
! IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - 
!
!
IF(line(1:3)=="705") THEN 
	call scan_string (line, 2, ss, nout)
	option_reb_gg = ss(1) 	
	option_reb_gr = ss(2) 
	If((option_reb_gg=='y'.or.option_reb_gr=='y').and.option_reb=='n') then
		       Write(88,*) 'The general switch for elastic rebound is turned OFF - try again...'  
		       Write(*, *) 'The general switch for elastic rebound is turned OFF - try again...'  		     	       
!		       
!		       write(*,*) option_reb
!		       write(*,*) option_reb_gg
!		       write(*,*) option_reb_gr
!		       
		       Call Stop_Config ; stop 	
	Endif		
	If(option_reb_gg=='y') Write(88,*) "Global elastic rebound maps for Greenland"	
	If(option_reb_gr=='y') Write(88,*) "Regional elastic rebound maps for Greenland"
	if((option_reb_gg=='y'.or.option_reb_gr=='y').and.ice_file(1:4).ne.'gree') then 
        	Write(88,*) 'You have selected an elastic rebound analysis for Greenland but ...'  
        	Write(88,*) '... the ice filename does not start by *gree* ==== Please check ===='  
        	Write(* ,*) 'You have selected an elastic rebound analysis for Greenland but ...'  
        	Write(* ,*) '... the ice filename does not start by *gree* ==== Please check ===='  
                Call Stop_Config ; stop 
	endif

ENDIF
!
IF(line(1:3)=="706") THEN 
    if( option_reb_gr == 'y' ) then
        call scan_string (line, 1, ss, nout)
	    res_reg_grc = ss(1) 	
        Write(88,*) "Resolution of regional elastic rebound maps for Greenland is "// &
                       trim(res_reg_grc),' km'
    else
        res_reg_grc='-1.0'
    endif
ENDIF
!
IF(line(1:3)=="715") THEN 
	call scan_string (line, 2, ss, nout)
	option_reb_ag = ss(1) 		
	option_reb_ar = ss(2) 	
	If((option_reb_ag=='y'.or.option_reb_ar=='y').and.option_reb=='n') then
		       Write(88,*) 'The general switch for elastic rebound is turned OFF - try again...'  
		       Write(*, *) 'The general switch for elastic rebound is turned OFF - try again...'  
		       Call Stop_Config ; stop 	
	Endif		
        If(option_reb_ag=='y') then 
		       Write(88,*) "Global elastic rebound maps for Antarctica" 		       
		       If(option_reb_gg=='y'.or.option_reb_gr=='y') then 
		       Write(88,*) 'Only ONE elastic rebound analysis is permitted!!!'  
		       Write(*, *) 'Only ONE elastic rebound analysis is permitted!!!'  
		       Call Stop_Config ; stop 
		       Endif		       
        Endif		       
	If(option_reb_ar=='y') then 
		       Write(88,*) "Regional elastic rebound maps for Antarctica"	
        	       If(option_reb_gg=='y'.or.option_reb_gr=='y') then 
		       Write(88,*) 'Only ONE elastic rebound analysis is permitted!!!'  
		       Write(*, *) 'Only ONE elastic rebound analysis is permitted!!!'  
		       Call Stop_Config ; stop 
		       Endif		       
	Endif
	if((option_reb_ag=='y'.or.option_reb_ar=='y').and.ice_file(1:4).ne.'anta') then 
        	Write(88,*) 'You have selected an elastic rebound analysis for Antarctica but ...'  
        	Write(88,*) '... the ice filename does not start by *anta* ==== Please check ===='  
        	Write(* ,*) 'You have selected an elastic rebound analysis for Antarctica but ...'  
        	Write(* ,*) '... the ice filename does not start by *anta* ==== Please check ===='  
                Call Stop_Config ; stop 
	endif
ENDIF
!
IF(line(1:3)=="716") THEN 
    if( option_reb_ar == 'y' ) then
        call scan_string (line, 1, ss, nout)
	    res_reg_grc = ss(1) 	
        Write(88,*) "Resolution of regional elastic rebound maps for Antarctica is "// &
                       trim(res_reg_anc),' km'
    else
        res_reg_anc='-1.0'
    endif
ENDIF
!
IF(line(1:3)=="730") THEN 
	call scan_string (line, 1, ss, nout)
	option_reb_sm = ss(1) 	
	If(option_reb_sm=='y'.and.option_reb=='n') then
		       Write(88,*) 'The general switch for elastic rebound is turned OFF - try again...'  
		       Write(*, *) 'The general switch for elastic rebound is turned OFF - try again...'  
		       Call Stop_Config ; stop 	
	Endif		
        If(option_reb_sm=='y') then 
        Write(88,*) "Global elastic rebound maps for small glaciers"    
        	      If(option_reb_ag=='y'.or.option_reb_ar=='y'.or. &       
			 option_reb_gg=='y'.or.option_reb_gr=='y')       then 
		      Write(88,*) 'Only ONE elastic rebound analysis is permitted!!!'  
		      Write(*, *) 'Only ONE elastic rebound analysis is permitted!!!'  
		      Call Stop_Config ; stop 
        	      Endif
        Endif
	if((option_reb_sm=='y').and.ice_file(1:4).ne.'glac') then 
        	Write(88,*) 'You have selected an elastic rebound analysis for glaciers but ...'  
        	Write(88,*) '... the ice filename does not start by *glac* ==== Please check ===='  
        	Write(* ,*) 'You have selected an elastic rebound analysis for glaciers but ...'  
        	Write(* ,*) '... the ice filename does not start by *glac* ==== Please check ===='  
                Call Stop_Config ; stop 
	endif

ENDIF
!
IF(line(1:3)=="740") THEN   		    
!
	call scan_string (line, 3, ss, nout)
!
	option_3d_reb = 'n'    
!
	option_3d_reb = ss(1) 
!
!
	If(option_3d_reb=='y') THEN 

        file_3d = './DATA/'//ss(2)	
!
	If(mode/='2') then  
		Write(88,*) "For elastic rebound, the SLE mode must be set to 2 (Elastic, GSC)"	
		Write(*, *) "For elastic rebound, the SLE mode must be set to 2 (Elastic, GSC)"	
		Call Stop_Config 
	Endif 
!	
	If(option_reb/='y') then 
		Write(88,*) "For elastic rebound at geodetic sites, the general swith must be ON"	
		Write(*, *) "For elastic rebound at geodetic sites, the general swith must be ON"	
		Call Stop_Config 
	Endif 
!	
	ELAREB_DATABASE_FORMAT=ss(3)
	
		Write(88,*) "Point estimates for elastic rebound"	
		Write(88,*) "Elastic rebound 3D velocities and \dot S and N at sites on file ", & 
		     	     trim(adjustl(file_3d))       
!
! ---   Counts the number of sites in database  
!
        	INQUIRE(FILE=FILE_3D,EXIST=lex)
!
        If(lex)then 
!
		LEN_GEOD=len(trim(adjustl(FILE_3D))) 
		GEO_DATABASE_ELAREB="'"//trim(adjustl(FILE_3D))//"'"
!
        	open(10,file=file_3d,status='old') 
!
!	
	If     (ELAREB_DATABASE_FORMAT=='0') THEN     ! ---------- Database of type '0'
!
        	ngeod=0
		do 33783 j=1, very_large_integer
		read(10,'(a200)',end=47401) row 
			do k=1, 200
			do n=1, NALFA
				if(row(k:k)==alfa(n)) then 
				ngeod=ngeod+1
				goto 33783
				endif 
			enddo
			enddo
		33783 continue
47401		close(10) 
!
	ELSEIf (ELAREB_DATABASE_FORMAT=='1') THEN    ! ---------- Database of type '1'
!
		ngeod=0
		do 33133 j=1, very_large_integer
			read(10,'(a200)',end=4741) row 
			if(row(1:1)/='#') ngeod=ngeod+1
		33133 continue
		4741 continue 
		CLOSE(10)
!
        ENDIF                                       ! ----- End of Database switch ...
!	
!
		Write(88,*) 'There are ', NGEOD, 'geodetic points in file ', FILE_3D
!
	else 
!
        Write(88,*) 'No geodetic database has been found...'	
        Write(*, *) 'No geodetic database has been found...'	
        Call Stop_Config ; stop	
!
	ENDIF   ! On the existence of the database ... 
!
	ENDIF   ! On the option ... 
!
			    ! WORK IN PROGRESS WORK IN PROGRESS WORK IN PROGRESS WORK IN PROGRESS 
ENDIF                       ! WORK IN PROGRESS WORK IN PROGRESS WORK IN PROGRESS WORK IN PROGRESS 
!
!
90009 CONTINUE
!






! ###### Elastic rebound ###### (multi-step melting ("pm" ice type) 
!
! [Only available for present ice sheets]
!
! >>>>>> Maps of (\dot S-U-N) during melting
!
!
IF(line(1:3)=="800") THEN 
!
	option_ela_ms = 'n'    
!
	call scan_string (line, 1, ss, nout)
!
	option_ela_ms = ss(1) 
!
	If(option_ela_ms=='y') then 
!	
        If(option_reb=='y') then 
                  Write(88,*) "The elastic (one-step) switch is ON- and it should be off"	
                  Write(* ,*) "The elastic (one-step) switch is ON- and it should be off"	
                  Call Stop_Config ; stop		
	Endif
!
	     If     (ice_type=='pm') then 
	          Write(88,*) "Multi-step analysis selected"	
	     elseif (ice_type/='pm') then 
                  Write(88,*) "For a multi-step analysis, the ice type should be <<pm>>"	
                  Write(* ,*) "For a multi-step analysis, the ice type should be <<pm>>"	
                  Call Stop_Config ; stop		
	     endif
!
        ENDIF
ENDIF
!
!
IF(line(1:3)=="830") THEN 
!
	option_gia_corr='n'
!
	call scan_string (line, 1, ss, nout)
!
	option_gia_corr= ss(1)

        If(option_ela_ms=='y')  then 
!		
		If    (option_gia_corr=='y') then
                  Write(88,*) "The multi-step analyses will be GIA-corrected"	
!                 Write(* ,*) "The multi-step analyses will be GIA-corrected"	             
	        elseif(option_gia_corr=='n') then 
                  Write(88,*) "The multi-step analyses will be NOT GIA-corrected"	
!                 Write(* ,*) "The multi-step analyses will be NOT GIA-corrected"	
		Endif		 			     
!
        Endif
!
ENDIF
!
!
IF(line(1:3)=="811") THEN    
!
	option_gm_ela_ms_s_var = 'n'     
	option_gm_ela_ms_s_dot = 'n'   
!
	call scan_string (line, 2, ss, nout)
!
	option_gm_ela_ms_s_var = ss(1)   
	option_gm_ela_ms_s_dot = ss(2)   
!
        If(option_gm_ela_ms_s_var=='y') then  
		Write(88,*) "Maps of RSL will be drawn"
		Endif
        If(option_gm_ela_ms_s_dot=='y') then   
		Write(88,*) "Maps of rate of RSL change will be drawn "
                Endif
!
ENDIF
IF(line(1:3)=="812") THEN    
!
	option_gm_ela_ms_n_var = 'n'     
	option_gm_ela_ms_n_dot = 'n'   
!
	call scan_string (line, 2, ss, nout)
!
	option_gm_ela_ms_n_var = ss(1)   
	option_gm_ela_ms_n_dot = ss(2)   
!
        If(option_gm_ela_ms_n_var=='y') then  
		Write(88,*) "Maps of ASL will be drawn"
		Endif
        If(option_gm_ela_ms_n_dot=='y') then   
		Write(88,*) "Maps of rate of ASL change will be drawn "
                Endif
ENDIF
IF(line(1:3)=="813") THEN       
!
	option_gm_ela_ms_u_var = 'n'     
	option_gm_ela_ms_u_dot = 'n'   
!
	call scan_string (line, 2, ss, nout)
!
	option_gm_ela_ms_u_var = ss(1)   
	option_gm_ela_ms_u_dot = ss(2)   
!
        If(option_gm_ela_ms_u_var=='y') then  
		Write(88,*) "Maps of U will be drawn"
		Endif
        If(option_gm_ela_ms_u_dot=='y') then   
		Write(88,*) "Maps of rates of U will be drawn "
                Endif
ENDIF
IF(line(1:3)=="814") THEN    
!
	option_gm_ela_ms_se_var = 'n'     
	option_gm_ela_ms_se_dot = 'n'   
!
	call scan_string (line, 2, ss, nout)
!
	option_gm_ela_ms_se_var = ss(1)   
	option_gm_ela_ms_se_dot = ss(2)   
!
        If(option_gm_ela_ms_se_var=='y') then  
		Write(88,*) "Maps of Eustatic Sea Level change will be drawn"
		Endif
        If(option_gm_ela_ms_se_dot=='y') then   
		Write(88,*) "Maps of rates of Eustatic Sea Level change will be drawn"
                Endif
!
ENDIF
!
!
IF(line(1:3)=="820") THEN 
!
! ... In progress - in progress - in progress - in progress - in progress  
!
	call scan_string (line, 3, ss, nout)
!
	OPTION_TG_ELA_MS = ss(1) 
!
        FILE_TG_ELA_MS = './DATA/'//ss(2)	
!
	If(OPTION_TG_ELA_MS=='y') THEN 
!
	If(mode/='2') then  
		Write(88,*) "For elastic rebound (in steps), the SLE mode must be set to 2 (Elastic, GSC)"	
		Write(*, *) "For elastic rebound (in steps), the SLE mode must be set to 2 (Elastic, GSC)"	
		Call Stop_Config 
	Endif 
	If(option_ela_ms/='y') then 
		Write(88,*) "For elastic rebound (in steps) at TGs, the general swith must be ON"	
		Write(*, *) "For elastic rebound (in steps) at TGs, the general swith must be ON"	
		Call Stop_Config 
	Endif 
!
	ELAREB_MS_DATABASE_FORMAT=ss(3)
!	
        Write(88,*) "Point estimates for elastic rebound (multi-steps)"	
        Write(88,*) "Name of database: ", trim(adjustl(FILE_TG_ELA_MS))	
	Write(88,*) "Format of database: ", ELAREB_MS_DATABASE_FORMAT
!
	if(ELAREB_MS_DATABASE_FORMAT/='1') THEN 
		Write(88,*) "Only files with format <<1>> are accepted for the moment"		
		Write(* ,*) "Only files with format <<1>> are accepted for the moment"		
        Call Stop_Config ; stop	
	ENDIF
!
!
! ---   Counts the number of sites in database  
!
        INQUIRE(FILE=FILE_TG_ELA_MS,EXIST=lex)
!
        If(lex)then 
!
		LEN_FILE_TG_ELA_MS = len(trim(adjustl(FILE_TG_ELA_MS))) 
		DATABASE_TG_ELA_MS = "'"//trim(adjustl(FILE_TG_ELA_MS))//"'"
!
        	open(10,file=FILE_TG_ELA_MS,status='old') 
!
	IF (ELAREB_MS_DATABASE_FORMAT=='1') THEN       !  Database of type '1'
!
		ngeod=0
		do 73733 j=1, very_large_integer
			read(10,'(a200)',end=4781) row 
			if(row(1:1)/='#') ngeod=ngeod+1
		73733 continue
		4781 continue 
		CLOSE(10)
!
        ENDIF                                       ! End of Database switch ...
!
        Write(88,*) 'There are ', NGEOD, 'geodetic points in file ', FILE_TG_ELA_MS

! ... In progress - in progress - in progress - in progress - in progress  
!
	else 
!
        Write(88,*) 'No database for the elastic multi-step analysis has been found...'	
        Write(*, *) 'No database for the elastic multi-step analysis has been found...'	
        Call Stop_Config ; stop	

        Endif
!
ENDIF
!
ENDIF
!
!
!
!
! ###### Present-day rates ###### 
!
! [Only available for LGM to Holocene ice sheets]
!
! >>>>>> *Global* maps of dot S, U, and N  
!
IF(line(1:3)=="290") THEN 
	call scan_string (line, 1, ss, nout)
	option_gm = ss(1) 
!	
	if(ice_type=='po'.or.ice_type=='pm')then 
		option_gm='n'
		goto 93039 
	endif
!
	If(option_gm=='y') Write(88,*) "Some GLOBAL maps will be drawn."	
ENDIF
!
!
! IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - 
! IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - 
!
IF(line(1:3)=="300") THEN 
!
	option_gmaps_xygrid = '1' 
!
	call scan_string (line, 1, ss, nout)
	option_gmaps_xygrid = ss(1) 
!
	if(option_gmaps_xygrid=='1') Write(88,*) "Global maps are computed on the pixelized grid"
	if(option_gmaps_xygrid=='2') Write(88,*) "Global maps are computed on a 1 deg x 1 deg grid"
!
!       write(*,*) option_gmaps_xygrid
!
	If(option_gmaps_xygrid/='1'.and.option_gmaps_xygrid/='2') then  
		Write(88,*) "For global maps, the grid option should be <<1>> or <<2>>"	
		Write(*, *) "For global maps, the grid option should be <<1>> or <<2>>"	
		Call Stop_Config 
	Endif 

ENDIF
!
!
! IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - 
! IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - 
!
IF(line(1:3)=="291") THEN 
	call scan_string (line, 1, ss, nout)
	option_dots = ss(1) 
	If(option_dots=='y'.and.option_gm=='y') Write(88,*) "A global map for dot-S will be drawn."
ENDIF
IF(line(1:3)=="292") THEN 
	call scan_string (line, 1, ss, nout)
	option_dotu = ss(1) 
	If(option_dotu=='y'.and.option_gm=='y') Write(88,*) "A global map for dot-U will be drawn."
ENDIF
IF(line(1:3)=="293") THEN 
	call scan_string (line, 1, ss, nout)
	option_dotn = ss(1) 
	If(option_dotn=='y'.and.option_gm=='y') Write(88,*) "A global map for dot-N will be drawn."
ENDIF
IF(line(1:3)=="299") THEN 
	call scan_string (line, 1, ss, nout)
	option_dotg = ss(1) 
	If(option_dotg=='y'.and.option_gm=='y') Write(88,*) "A global map for dot-Phi/gamma will be drawn."
ENDIF
IF(line(1:3)=="301") THEN 
	call scan_string (line, 1, ss, nout)
	option_dotw = ss(1) 
	If(option_dotw=='y'.and.option_gm=='y') & 
	Write(88,*) "A global map for the <<equivalent water height>> is drawn "   ! ------- In progress ...
ENDIF
IF(line(1:3)=="294") THEN 
	call scan_string (line, 1, ss, nout)
	option_dotfa = ss(1) 
	If(option_dotfa=='y'.and.option_gm=='y') Write(88,*) "A global map for dot-FA will be drawn."
ENDIF
IF(line(1:3)=="295") THEN 
	call scan_string (line, 1, ss, nout)
	option_dotss = ss(1) 
	If(option_dotss=='y'.and.option_gm=='y') Write(88,*) "A global map for dot-SS will be drawn."
ENDIF
IF(line(1:3)=="296") THEN 
	call scan_string (line, 1, ss, nout)
	option_dotloi = ss(1) 
	If(option_dotloi=='y'.and.option_gm=='y')Write(88,*) "A global map for dot-LOAD-Ice will be drawn."
ENDIF
IF(line(1:3)=="297") THEN 
	call scan_string (line, 1, ss, nout)
	option_dotloo = ss(1) 
	If(option_dotloo=='y'.and.option_gm=='y')Write(88,*) "A global map for dot-LOAD-Ocean will be drawn."
ENDIF
IF(line(1:3)=="298") THEN 
	call scan_string (line, 1, ss, nout)
	option_dotlot = ss(1) 
	If(option_dotlot=='y'.and.option_gm=='y')Write(88,*) "A global map for dot-LOAD will be drawn."
ENDIF
93039   Continue
!
!
! >>>>>> *Regional* maps of dot S, U, and N  
!
! [Only available for LGM to Holocene ice sheets]
!
IF(line(1:3)=="500") THEN 
	call scan_string (line, 1, ss, nout)
	option_rm(0) = ss(1) 
!
	if(ice_type=='po'.or.ice_type=='pm')then 
		option_rm(0)='n'
		goto 91019 
	endif
!
	If(option_rm(0)=='y'.and.option_rof=='r') then 
	Write(88,*) "\dot S, U, and N at present-time will be computed..."
	Write(88,*) "... and some REGIONAL maps will be drawn."
	Endif
!
	If(option_rm(0)=='y'.and.option_rof=='z') then 
	Write(88,*) "NO REGIONAL maps for a zonal OF"
	option_rm(0)='n'	
	Endif
ENDIF
!
!
IF(line(1:3)=="501") THEN 
	call scan_string (line, 1, ss, nout)
	option_rm(1) = ss(1) 
	If(option_rm(1)=='y'.and.option_rm(0)=='y') & 
	Write(88,*) "Regional Map of \dot S, U, and N for ITALY"
ENDIF
IF(line(1:3)=="502") THEN 
	call scan_string (line, 1, ss, nout)
	option_rm(2) = ss(1) 
	If(option_rm(2)=='y'.and.option_rm(0)=='y') & 
        Write(88,*) "Regional Map of \dot S, U, and N for the MEDITERRANEAN"
ENDIF
IF(line(1:3)=="503") THEN 
	call scan_string (line, 1, ss, nout)
	option_rm(3) = ss(1) 
	If(option_rm(3)=='y'.and.option_rm(0)=='y') & 
	Write(88,*) "Regional Map of \dot S, U, and N for EUROPE"
ENDIF
IF(line(1:3)=="504") THEN 
	call scan_string (line, 1, ss, nout)
	option_rm(4) = ss(1) 
	If(option_rm(4)=='y'.and.option_rm(0)=='y') & 	
	Write(88,*) "Regional Map of \dot S, U, and N for FENNOSCANDIA"
ENDIF
IF(line(1:3)=="505") THEN 
	call scan_string (line, 1, ss, nout)
	option_rm(5) = ss(1) 
        If(option_rm(5)=='y'.and.option_rm(0)=='y') & 
	Write(88,*) "Regional Map of \dot S, U, and N for GREENLAND"
ENDIF
IF(line(1:3)=="506") THEN 
	call scan_string (line, 1, ss, nout)
	option_rm(6) = ss(1) 
	If(option_rm(6)=='y'.and.option_rm(0)=='y') & 
        Write(88,*) "Regional Map of \dot S, U, and N for NORTH AMERICA"
ENDIF
IF(line(1:3)=="507") THEN 
	call scan_string (line, 1, ss, nout)
	option_rm(7) = ss(1) 
	If(option_rm(7)=='y'.and.option_rm(0)=='y') & 
	Write(88,*) "Regional Map of \dot S, U, and N for ANTARCTICA"
ENDIF
91019   CONTINUE
!
!
! ###### 3D VELOCITY ###### 
!
! >>>>>> 3D velocity and dot S & N at user-supplied geodetic sites...  
!
!
IF(line(1:3)=="270") THEN 
!
	call scan_string (line, 2, ss, nout)
!
	option_3d = 'n'
!
	option_3d = ss(1) 
!	
        If(option_3d=='y') THEN 
!
	If(ice_type/='ho') then   
		Write(88,*) "3D velocity at sites is only available for Holocene ice"	
		Write(*, *) "3D velocity at sites is only available for Holocene ice"	
		Call Stop_Config 
	Endif 
!	
	file_3dd = './DATA/'//trim(adjustl(ss(2))) 
	
!write(*,*) trim(adjustl(file_3dd)) 
!
        Write(88,*) "Point estimates for Glacial Isostatic Adjustment (GIA)"	   
	Write(88,*) "3D velocities and \dot S and N at sites on file ", & 
		     trim(adjustl(file_3dd))
!
! ---   Counts the number of sites in database  
!
        INQUIRE(FILE=FILE_3DD,EXIST=lex)
!
!write(*,*) LEX
!
        If(lex)then 
!
		LEN_GEOD=len(trim(adjustl(FILE_3DD))) 
		GEO_DATABASE="'"//trim(adjustl(FILE_3DD))//"'"
!
        	open(10,file=file_3dd,status='old') 
!
        	ngeod=0
		do 33733 j=1, very_large_integer
		read(10,'(a200)',end=474) row 
			do k=1, 200
			do n=1, NALFA
				if(row(k:k)==alfa(n)) then 
				ngeod=ngeod+1
				goto 33733
				endif 
			enddo
			enddo
		33733 continue
474	close(10) 
!
        Write(88,*) 'There are ', NGEOD, 'geodetic points in file ', FILE_3DD
!
        else
!	
        Write(88,*) 'No geodetic database has been found...'	
        Write(*, *) 'No geodetic database has been found...'	
	
	write(*,*) file_3dd
	
        Call Stop_Config ; stop	
!
	ENDIF   ! On the existence of the database ... 
!
	ENDIF   ! On the option ... 
!
 ENDIF 
!
! >>>>>> 3D velocity on maps of interesting regions
!
	IF(line(1:3)=="275") THEN 
!
		call scan_string (line, 2, ss, nout)
!
		option_3d_regions = ss(1) 
!		
		file_3d_regions = './DATA/'//ss(2)
!		
		If(option_3d_regions=='y') then  
!		
                INQUIRE(FILE=FILE_3D_REGIONS,EXIST=lex)     
!
                      If(lex)then 
!		      
	              LEN_3D_REGIONS=len(trim(adjustl(FILE_3D_REGIONS))) 
!		   
	              TRED_REGIONS_DATABASE="'"//trim(adjustl(FILE_3D_REGIONS))//"'"		      
		  		
			Write(88,*) "3D velocity fields maps across interesting regions will be produced"
			Write(88,*) "The regions are listed in file ", trim(adjustl(file_3d_regions))
!
		        Call scan_3d_regions (file_3d_regions,	&    
					     n_3d_regions,	& 
					     name_3d_regions,	& 
					     np_3d_regions,	& 
					     pix_3d_filenames)
!
        		If(n_3d_regions==0) then 
				Write(88,*) 'No regions selected in file ', trim(adjustl(FILE_3D_REGIONS)) 
				Write(*, *) 'No regions selected in file ', trim(adjustl(FILE_3D_REGIONS)) 
        			Call Stop_Config 
			Endif
!
			do j=1, n_3d_regions
				TRED_REGIONS_NAME(j)     = "'"//trim(adjustl(name_3d_regions(j)))//"'"
				LEN_TRED_REGIONS_NAME(j) = len(trim(adjustl(name_3d_regions(j))))
			enddo
!
			endif
		endif
!
	ENDIF
!
!
!
! ###### Stokes coefficients (SC) ###### 
!
! >>>>>> Rate of change of SC, and range of degrees  
!
IF(line(1:3)=="280") THEN 
	call scan_string (line, 3, ss, nout)
	option_st     = ss(1)
	degree_st_min = ss(2)
	degree_st_max = ss(3)
	If(option_st=='y') Write(88,*) 'Range of degrees for the Stokes coefficients: ', & 
	    	  	   trim(degree_st_min), '-', trim(degree_st_max)
ENDIF
!
!
200 CONTINUE
!
555 close(1)
!
!
!
!
!
! ====================================
!
! 	Part #2 : building "selen.sh"
!
! 	1/2: Compilation commands 
! 	2/2:  Execution commands
!
! ====================================
!
!
!
  Write(88,*) ""
  Write(88,*) "Compilation and execution commands are written on selen.sh"
!
  open(2,file='selen.sh',status='unknown')
!
! 
! >>>>>> A time stamp on "selen.sh"
!
  call DATE_AND_TIME (date,timc)      
  Write(2,*) " "
  Write(2,*) '# This is file "selen.sh", created by "config.f90" on ', & 
             date(1:4), '.', date(5:6), '.', date(7:8), ' ', & 
	     timc(1:2), '.', timc(3:4), '.', timc(5:6) 
  Write(2,*) "  "
!
!
! >>>>>>  # Part 1/2: compilation 
!
Write(2,*) "echo ''"
Write(2,*) "echo '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -'"
Write(2,*) "echo '     SELEN, a Sea levEL EquatioN solver, Version 3.2     '"
Write(2,*) "echo '                                                         '" 
Write(2,*) "echo '      Web page: http://fcolleoni.free.fr/SELEN.html      '"	    
Write(2,*) "echo '   Send comments, requests of help and suggestions to:   '" 
Write(2,*) "echo '                <giorgio.spada@gmail.com>                '"
Write(2,*) "echo '                            -                            '"
Write(2,*) "echo '                    Copyright(C) 2008                    '"    
Write(2,*) "echo '     Giorgio Spada, Florence Colleoni & Paolo Stocchi    '"
Write(2,*) "echo '                          * * *                          '"
Write(2,*) "echo '     This programs comes with  ABSOLUTELY NO WARRANTY    '"
Write(2,*) "echo ' This is free software and you are welcome to distribute '"
Write(2,*) "echo '              it under certain conditions.               '"
Write(2,*) "echo '    For details, visit  <http://www.gnu.org/licenses/>   '"
Write(2,*) "echo '                  or edit file COPYING                   '"
Write(2,*) "echo '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -'"
Write(2,*) "echo ''"
!
Write(2,*) " "
Write(2,*) " "
Write(2,*) "echo '----------------------------------------'"
Write(2,*) "echo ' >>> 0. Compiling the SELEN programs ...'"
Write(2,*) "echo '----------------------------------------'"
!
!
! --- Compilation: SHTOOLS module
 Write(2,*) "echo '--- SHTOOLS.F90'"
 Write(2,*) CompileSeq, " -c shtools.f90" 
!
!
! --- Compilation: HARMONICS module
 Write(2,*) "echo '--- HARMONICS.F90'"
 Write(2,*) CompileSeq, " -c harmonics.f90 -o harmonics.o" 
 Write(2,*) CompileSmp, " -c harmonics.f90 -o harmonics-smp.o" 
!
!
!
!
! --- Communicating with ALMA 
Write(2,*) "cp working-directory.txt ./ALMA"
!


!If(option_rfb=='y'.and.option_tlove_ext=='n') then 
!	Write(2,*) "echo '--- S_ROT.F90'"
!	Write(2,*) CompileSeq, " s_rot.f90 harmonics.o -o s_rot.exe" 
!    Write(2,*) "echo '--- U_ROT.F90'"
!    Write(2,*) CompileSeq, " u_rot.f90 harmonics.o -o u_rot.exe" 
!    Write(2,*) "echo '--- N_ROT.F90'"
!    Write(2,*) CompileSeq, " n_rot.f90 harmonics.o -o n_rot.exe" 
!Endif	

!If(option_rfb=='y'.and.option_tlove_ext=='y') then 
!	Write(2,*) "echo '--- S_ROT_PM.F90'"
!	Write(2,*) CompileSeq, " s_rot_pm.f90 harmonics.o -o s_rot_pm.exe" 
!    Write(2,*) "echo '--- U_ROT_PM.F90'"
!    Write(2,*) CompileSeq, " u_rot_pm.f90 harmonics.o -o u_rot_pm.exe" 
!    Write(2,*) "echo '--- N_ROT_PM.F90'"
!    Write(2,*) CompileSeq, " n_rot_pm.f90 harmonics.o -o n_rot_pm.exe" 
!    Endif	

!
! --- Compilation: SLE
!
If(option_topo=='n') then
!
If(option_rfb=='n') then 
!
	if((ice_type=='ho').or.(ice_type=='po')) then 
	     Write(2,*) "echo '--- SLE.F90'"
	     Write(2,*) CompileSmp, " sle.f90 harmonics-smp.o -o sle.exe -O3" 
        endif
	if(ice_type=='pm') then 
	     Write(2,*) "echo '--- SLE_PM.F90'"
	     Write(2,*) CompileSmp, " sle_pm.f90 harmonics-smp.o -o sle.exe -O3" 
        endif
!
else
!	
	if(ice_type=='ho') then 
	     Write(2,*) "echo '--- SLE_ROTAZ_NEW.F90'"
	     Write(2,*) CompileSmp, " sle_rotaz_new.f90 harmonics-smp.o -o sle.exe -O3" 
        endif
	if(ice_type=='pm') then 
	     Write(2,*) "echo '--- SLE_ROTAZ_NEW_PM.F90'"
	     Write(2,*) CompileSmp, " sle_rotaz_new_pm.f90 harmonics-smp.o -o sle.exe -O3" 
        endif

Endif	
!
Else ! If(option_topo=='y') then  
!
! --- Compilation: Evolving coastlines 
!
!
If(option_pxtopo=='y') then 
    Write(2,*) "echo '--- PX_TOPO.F90'" 
    Write(2,*) CompileSmp, " px_topo.f90 harmonics-smp.o -o pxtopo.exe -O3" 
endif 
!
Write(2,*) "echo '--- ICE_MASK.F90'" 
Write(2,*) CompileSmp, " ice_mask.f90 harmonics-smp.o -o icemask.exe" 
!
Write(2,*) "echo '--- RSL_CPX.F90'" 
Write(2,*) CompileSmp, " rsl_cpx.f90 harmonics-smp.o -o rslcpx.exe" 
!
Write(2,*) "echo '--- SH_OF_VAROC.F90'" 
Write(2,*) CompileSmp, " sh_of_varoc.f90 harmonics-smp.o -o shofvaroc.exe" 
!
!
If(option_rfb=='n') then
!
Write(2,*) "echo '--- SLE.F90'"
Write(2,*) CompileSmp, " sle.f90 harmonics-smp.o -o sle.exe -O3"
!
Write(2,*) "echo '--- SLE_VAROC.F90'"
Write(2,*) CompileSmp, " sle_varoc.f90 harmonics-smp.o -o slevaroc.exe" 
!
Else
!
Write(2,*) "echo '--- SLE_ROTAZ_NEW.F90'"
Write(2,*) CompileSmp, " sle_rotaz_new.f90 harmonics-smp.o -o sle.exe -O3" 
!
Write(2,*) "echo '--- SLE_VAROC_ROTAZ.F90'" 
Write(2,*) CompileSmp, " sle_varoc_rotaz.f90 harmonics-smp.o -o slevaroc.exe" 
!
Endif
!
endif
!
! --- Compilation: Ice sheets ----
!
! --- Computation of ice Shape factors and SH decomposition 
If(option_sf=='y') then    
!
Write(2,*) "echo '--- SHAPE_FACTORS_MOD.F90'" 
Write(2,*) CompileSmp, " shape_factors_mod.f90 harmonics-smp.o -o shapefactors.exe" 
!
Write(2,*) "echo '--- SHICE_MOD.F90'"
Write(2,*) CompileSmp, " shice_mod.f90 harmonics-smp.o -o shice.exe"
!
! Write(2,*) "echo '--- SHAPE_FACTORS.F90'" 
! Write(2,*) CompileSmp, " shape_factors.f90 -o shapefactors.exe" 
!
Endif
!
! --- Equivalent sea level for "ho" models 
If     (option_esl=='y'.and.option_topo=='n') then 
	Write(2,*) "echo '--- ESL.F90'"
	Write(2,*) CompileSeq, " esl.f90 harmonics.o -o esl.exe"
elseif (option_esl=='y'.and.option_topo=='y') then
	Write(2,*) "echo '--- ESL.F90'"
	Write(2,*) CompileSeq, " esl.f90 harmonics.o -o esl.exe"
endif 
!
!
!
! --- Equivalent sea level for "pm" models 
If     (option_esl_pm=='y') THEN 
                      Write(2,*) "echo '--- ESL_PM.F90'"
                      Write(2,*) CompileSeq, " esl_pm.f90 harmonics.o -o esl_pm.exe"
endif 
!
!
! --- SH decomposition 
!  Write(2,*) "echo '--- SHICE.F90'"
!  Write(2,*) Compile, " SHICE.F90 -o shice.exe"
!
!
! --- Ice sheets contours 
If(option_or=='y') then 
Write(2,*) "echo '--- MS.F90'"
Write(2,*) CompileSeq, " ms.f90 harmonics.o -o ms.exe"
endif
!
!=======================================
! --- Compilation of Love numbers tools 
!=======================================
If    (option_nm=='y'.and.option_tlove=='n') then 
!
! --- TABOO, Loading analysis 
!
Write(2,*) "echo '--- TB_MOD.F90' - standard TABOO environment"
Write(2,*) CompileSeq, " tb_mod.F90 harmonics.o -o tb.exe"  
!
Endif
!
If    (option_nm=='y'.and.option_tlove=='y') then 
!
! --- TABOO, Loading analysis & Tidal analysis, using a HP TABOO version...
!
Write(2,*) "echo '--- TBHP.F90' - MULTI-precision TABOO environment"
Write(2,*) "export FORTCMD=""",trim(CompileSeq),""""
Write(2,*) "sh FMLIB_compile.sh"
Write(2,*) CompileSeq, " tbhp.f90 -c -o tbhp.o"
Write(2,*) CompileSeq, " tbhp.o harmonics.o FMSAVE.o FMZM90.o FM.o -o tb.exe"
!
Endif
!
!
If(option_pw=='y') then 
!
! --- ALMA 
Write(2,*) "echo '--- alma.f90'"
!
Write(2,*) "if [ ! -f ", trim(adjustl(wdir))//"/ALMA/alma.inc", " ]" 
Write(2,*) "then"
Write(2,*) "/bin/rm -v", trim(adjustl(wdir))//"/ALMA/alma.inc"
Write(2,*) "fi"
! --------------- DM 06.11.2012 ---------------
Write(2,*) "export FORTCMD=""",trim(CompileSeq),""""
Write(2,*) "cd ALMA ; sh FMLIB_compile.sh ; cd .."
! ---------------------------------------------

Write(2,*) CompileSeq, trim(adjustl(wdir))//"/ALMA/alma.f90 -c -o ", trim(adjustl(wdir))//"/ALMA/alma.o"
Write(2,*) CompileSeq, trim(adjustl(wdir))//"/ALMA/alma.o ",   & 
                       trim(adjustl(wdir))//"/ALMA/FMSAVE.o ", & 
		       trim(adjustl(wdir))//"/ALMA/FMZM90.o ", & 
		       trim(adjustl(wdir))//"/ALMA/FM.o -o ", & 
		       trim(adjustl(wdir))//"/ALMA/alma.exe"
!
Write(2,*) "echo '--- AM.F90'"
Write(2,*) CompileSeq, " am.f90 -o am.exe"  
!
elseif(option_pwa=='y') then 
!
Write(2,*) "echo '--- AM.F90'"
Write(2,*) CompileSeq, " am.f90 -o am.exe"  
!
endif 
!
!
! --- Compilation: Polar Motion Transfer Function 



If(option_pmtf=='y'.and.option_tlove_ext=='n') then  
Write(2,*) "echo '--- PM_v6.F90'"      
Write(2,*) CompileSeq,  " PM_v6.f90 -o pmtf.exe" 
Endif
!
If(option_pmtf=='y'.and.option_tlove_ext=='y') then  
!
! Nothing to do, in this case... 
!
Endif
!
If(option_pmd=='y'.and.option_tlove_ext=='n') then  
Write(2,*) "echo '--- PMD_MOD.F90'"      
Write(2,*) CompileSeq,  " pmd_mod.f90 harmonics.o -o pmd.exe" 
!
Write(2,*) "echo '--- LOAD_RE.F90'"      
Write(2,*) CompileSeq,  " load_re.f90 harmonics.o -o loadre.exe" 
Endif
!
!
If(option_pmd=='y'.and.option_tlove_ext=='y') then  
!
Write(2,*) "echo '--- LOAD_RE.F90'"      
Write(2,*) CompileSeq,  " load_re.f90 harmonics.o -o loadre.exe" 
!
Endif
!
if (option_npx=='y') then   ! -------------- A NEW PIXEL TABLE
!
! --- Compilation: Pixelization (i) 
Write(2,*) "echo '--- PX.F90'"      
Write(2,*) CompileSeq,  " px.f90 harmonics.o -o px.exe" 
!
! --- Compilation: Pixelization (ii) 
Write(2,*) "echo '--- PX_REC.F90'"      
Write(2,*) CompileSeq,  " px_rec.f90 -o pxrec.exe" 
!
else                        ! -------------- AN EXISTING PIXEL TABLE
!
! --- Compilation: Retrieve pixelization info from existing table
Write(2,*) "echo '--- PX_REBUILD.f90'"      
Write(2,*) trim(CompileSeq), " px_rebuild.f90 harmonics.o -o pxrebuild.exe" 
!
endif
!
! --- Compilation: Spherical harmonics 
If(option_sh=='y') then 
Write(2,*) "echo '--- SH.F90'"      
Write(2,*) CompileSmp, " sh.f90 harmonics-smp.o -o sh.exe" 
Endif
!
! --- Compilation: Window function
If(option_wi=='y') then 
Write(2,*) "echo '--- WNW.F90'"      
Write(2,*) CompileSeq, " wnw.f90 harmonics.o -O3 -o wnw.exe" 
Endif
!     
! --- Compilation: Ocean function (OF) harmonics 
If(option_oh=='y') then 
Write(2,*) "echo '--- SH_OF.F90'"   
Write(2,*) CompileSmp, " sh_of.f90 harmonics-smp.o -o shof.exe"   
Endif 
!
! --- Compilation: OF DV computation  
If(option_ofdv=='y') then 
Write(2,*) "echo '--- OF_DV.F90'"  
Write(2,*) CompileSeq, " of_dv.f90 harmonics.o -o ofdv.exe"
Endif
!
! --- Compilation: OF Reconstruction 
If(option_of=='y') then 
Write(2,*) "echo '--- REC_OF.F90'"  
Write(2,*) CompileSmp, " rec_of.f90 harmonics-smp.o -o recof.exe"
Endif
!
! --- Compilation: Ice sheets Reconstruction 
If(option_ri=='y') then 
Write(2,*) "echo '--- REC_ICE.F90'" 
Write(2,*) CompileSmp, " rec_ice.f90 harmonics-smp.o -o recice.exe"
Endif 
!
! --- Compilation: Elastic rebound  
If(option_reb=='y') then 
Write(2,*) "echo '--- ELA_REB.F90'"
Write(2,*) CompileSmp, " ela_reb.f90 harmonics-smp.o -o elareb.exe"
Endif
!
! --- Compilation: Elastic rebound with geodetic quantities 
If(option_3d_reb=='y') then 
Write(2,*) "echo '--- GEO_ELA_REB.F90'"
Write(2,*) CompileSeq, " geo_ela_reb.f90 harmonics.o -o geoelareb.exe"
Endif
!
! WORK IN PROGRESS WORK IN PROGRESS WORK IN PROGRESS WORK IN PROGRESS WORK IN PROGRESS - - - - - - - - - 
!
! --- Compilation: Global maps of elastic rebound with multi-steps (=variations=)
!
 If(option_ela_ms=='y') then 
!
 If(option_gm_ela_ms_s_var=='y'.or. & 
    option_gm_ela_ms_n_var=='y'.or. &
    option_gm_ela_ms_u_var=='y') then   
    Write(2,*) "echo '--- GMAPS_ELA_MS_VAR.F90'"
    Write(2,*) CompileSeq, " gmaps_ela_ms_var.f90 harmonics.o -o gmaps_ela_ms_var.exe"  
 Endif
!
! --- Compilation: Global maps of elastic rebound with multi-steps (=rates=)
!
 If(option_gm_ela_ms_s_dot=='y'.or. & 
    option_gm_ela_ms_n_dot=='y'.or. &
    option_gm_ela_ms_u_dot=='y') then   
    Write(2,*) "echo '--- GMAPS_ELA_MS_DOT.F90'"
    Write(2,*) CompileSeq, " gmaps_ela_ms_dot.f90 harmonics.o -o gmaps_ela_ms_dot.exe"  
 Endif
!
 Endif
!
! --- Compilation: SH at tide gauges           (for elastic rebound with multi-steps)
If(OPTION_TG_ELA_MS=='y') then
   Write(2,*) "echo '--- SH_TGAUGES.F90'" 
   Write(2,*) CompileSeq, " sh_tgauges.f90 harmonics.o -o shtgauges.exe" 
!
! --- Compilation: Predictions at tide gauges  (for elastic rebound with multi-steps)
   Write(2,*) "echo '--- TGAUGES_ELA_MS.F90'" 
   Write(2,*) CompileSeq, " tgauges_ela_ms.f90 harmonics.o -o tgaugeselams.exe"
endif
!
! WORK IN PROGRESS WORK IN PROGRESS WORK IN PROGRESS WORK IN PROGRESS WORK IN PROGRESS - - - - - - - - - -
!
! --- Compilation: Global maps 
If(option_gm=='y') then 
   Write(2,*) "echo '--- GMAPS.F90'"
   Write(2,*) CompileSmp, " gmaps.f90 harmonics-smp.o -o gmaps.exe"
Endif
!
! --- Compilation: Regional maps 
If(option_rm(0)=='y') then 
   Write(2,*) "echo '--- RMAPS.F90'"
   Write(2,*) CompileSmp, " rmaps.f90 harmonics-smp.o -o rmaps.exe"
Endif
!
! --- Compilation: 3D veolcity and \dot S, N, U at sites
If(option_3d=='y') then 
   Write(2,*) "echo '--- GEO.F90'"
   Write(2,*) CompileSeq, " geo.f90 harmonics.o -o geo.exe"
Endif
!
! --- Compilation: 3D veolcity and \dot U at pixels on maps
If(option_3d_regions=='y') then 
   Write(2,*) "echo '--- GEO_MAPS.F90'"
   Write(2,*) CompileSeq, " geo_maps.f90 harmonics.o -o geomaps.exe"
Endif
!
! --- Compilation: SH at Relative Sea Level (RSL) sites 
If(option_rsl=='y') then
   Write(2,*) "echo '--- SH_RSL.F90'"
   Write(2,*) CompileSeq, " sh_rsl.f90 harmonics.o -o shrsl.exe" 
!
! --- Compilation: RSL at RSL sites 
   Write(2,*) "echo '--- RSL.F90'" 
   Write(2,*) CompileSeq, " rsl.f90 harmonics.o -o rsl.exe"
Endif
!
! --- Compilation: RSL Zones 
If(option_rslz=='y') then
   Write(2,*) "echo '--- RSL_ZONES.F90'"
   Write(2,*) CompileSeq, " rsl_zones.f90 harmonics.o -o rslzones.exe" 
Endif
!
! --- Compilation: SH at virtual sites for RSL contours
If(option_rslc=='y') then
   Write(2,*) "echo '--- SH_RSLC.F90'"
   Write(2,*) CompileSeq, " sh_rslc.f90 harmonics.o -o shrslc.exe" 
!
! --- Compilation: RSL at virtual sites for RSL contours 
   Write(2,*) "echo '--- RSLC.F90'" 
   Write(2,*) CompileSeq, " rslc.f90 harmonics.o -o rslc.exe"
Endif
!
! --- Compilation: SH at tide gauges
If(option_tg=='y') then
   Write(2,*) "echo '--- SH_TGAUGES.F90'" 
   Write(2,*) CompileSeq, " sh_tgauges.f90 harmonics.o -o shtgauges.exe" 
!
! --- Compilation: SL change at tide gauges
   Write(2,*) "echo '--- TGAUGES.F90'" 
   Write(2,*) CompileSeq, " tgauges.f90 harmonics.o -o tgauges.exe"
endif
!
! --- Stokes coefficients
If(option_st=='y') then
   Write(2,*) "echo '--- STOKES.F90'"
   Write(2,*) CompileSeq, " stokes.f90 harmonics.o -o stokes.exe" 
Endif
!	      
! --- End of compilation
!
!
!
Write(2,*) " "
Write(2,*) " "
Write(2,*) "echo ''"
Write(2,*) "echo '----------------------------'"
Write(2,*) "echo ' >>> 1. Executing SELEN  ...'"
Write(2,*) "echo '----------------------------'"
!
!
!
! ===================================================
! --- Creating folders into the ./depot directory ---
! ===================================================
!
! ================================
! EXE 00 --- Working directory --- 
! ================================
Write(2,*) " "
Write(2,*) "#echo --------------------------------------------------------"
Write(2,*) "echo                                                          "
Write(2,*) " echo '---> Working directory: '", trim(adjustl(wdir)) 
Write(2,*) "#echo --------------------------------------------------------"


Write(2,*) " "
Write(2,*) "#echo ------------------------------------------------------"
Write(2,*) "echo"
Write(2,*) " echo '---> Output data will be stored into directory '", trim(adjustl(depot))
Write(2,*) "#echo ------------------------------------------------------"
Write(2,*) "if [ ! -d ", depot, " ]"
Write(2,*) "then"
Write(2,*) "mkdir ",depot//"/" 
Write(2,*) "mkdir ",depot//"/MODULES/" 
Write(2,*) "mkdir ",depot//"/bin/" 
Write(2,*) "mkdir ",depot//"/px/" 
Write(2,*) "mkdir ",depot//"/wnw/"
Write(2,*) "mkdir ",depot//"/log/"	
Write(2,*) "mkdir ",depot//"/of/"
Write(2,*) "mkdir ",depot//"/of/degree_variance/" 
Write(2,*) "mkdir ",depot//"/paleo-topography/" 
Write(2,*) "mkdir ",depot//"/paleo-topography/data/" 
Write(2,*) "mkdir ",depot//"/paleo-topography/more-data/" 
Write(2,*) "mkdir ",depot//"/paleo-topography/ps/" 		 		
Write(2,*) "mkdir ",depot//"/paleo-topography/pdf/" 		 			
Write(2,*) "mkdir ",depot//"/"//trim(adjustl(titlice))//"/"
Write(2,*) "mkdir ",depot//"/"//trim(adjustl(titlice))//"/original/"
Write(2,*) "mkdir ",depot//"/"//trim(adjustl(titlice))//"/esl/"
Write(2,*) "mkdir ",depot//"/"//trim(adjustl(titlice))//"/reconstructed/"	
Write(2,*) "mkdir ",depot//"/"//trim(adjustl(titlice))//"/sh/"
Write(2,*) "mkdir ",depot//"/"//trim(adjustl(titlice))//"/esl-pm/"
Write(2,*) "mkdir ",depot//"/TABOO/"
Write(2,*) "mkdir ",depot//"/Love-Numbers-by-TABOO/"
Write(2,*) "mkdir ",depot//"/ALMA/"
Write(2,*) "mkdir ",depot//"/Love-Numbers-by-ALMA/"
Write(2,*) "mkdir ",depot//"/PMTF/"
Write(2,*) "mkdir ",depot//"/PM-direct-effect/"
Write(2,*) "mkdir ",depot//"/rsl/"
Write(2,*) "mkdir ",depot//"/rsl/rsl-sites/"
Write(2,*) "mkdir ",depot//"/rsl/rsl-curves/"
Write(2,*) "mkdir ",depot//"/rsl/rsl-curves/ps/"
Write(2,*) "mkdir ",depot//"/rsl/rsl-curves/pdf/"
Write(2,*) "mkdir ",depot//"/rsl/rsl-scplot/"
Write(2,*) "mkdir ",depot//"/rsl/rsl-misfit/"
Write(2,*) "mkdir ",depot//"/rsl/rsl-table-obs-vs-predictions/"
Write(2,*) "mkdir ",depot//"/rsl/rsl-zones/"
Write(2,*) "mkdir ",depot//"/rsl/rsl-contours/"
Write(2,*) "mkdir ",depot//"/tgauges/"
Write(2,*) "mkdir ",depot//"/tgauges/tgauges-sites/"
Write(2,*) "mkdir ",depot//"/tgauges/tgauges-scplots/"
Write(2,*) "mkdir ",depot//"/tgauges/tgauges-predictions/"
Write(2,*) "mkdir ",depot//"/elastic-rebound/"
Write(2,*) "mkdir ",depot//"/tgauges-ela-ms/"
Write(2,*) "mkdir ",depot//"/elastic-rebound/Geodetic/"
Write(2,*) "mkdir ",depot//"/elastic-rebound/Greenland/"
Write(2,*) "mkdir ",depot//"/elastic-rebound/Antarctica/"
Write(2,*) "mkdir ",depot//"/elastic-rebound/Small_glaciers/"
Write(2,*) "mkdir ",depot//"/gmaps/"
Write(2,*) "mkdir ",depot//"/gmaps-ela-ms/"
Write(2,*) "mkdir ",depot//"/gmaps-ela-ms/var"
Write(2,*) "mkdir ",depot//"/gmaps-ela-ms/dot"
Write(2,*) "mkdir ",depot//"/geod/"
Write(2,*) "mkdir ",depot//"/geod/sites/"
Write(2,*) "mkdir ",depot//"/geod/3dmaps/"
Write(2,*) "mkdir ",depot//"/geod/3dmaps/northamerica/"
Write(2,*) "mkdir ",depot//"/geod/3dmaps/fennoscandia/"
Write(2,*) "mkdir ",depot//"/geod/3dmaps/mediterranean/"
Write(2,*) "mkdir ",depot//"/rmaps/"
Write(2,*) "mkdir ",depot//"/rmaps/sun_data"
Write(2,*) "mkdir ",depot//"/rmaps/Italy"
Write(2,*) "mkdir ",depot//"/rmaps/Mediterranean"
Write(2,*) "mkdir ",depot//"/rmaps/Europe"
Write(2,*) "mkdir ",depot//"/rmaps/Fennoscandia"
Write(2,*) "mkdir ",depot//"/rmaps/Greenland"
Write(2,*) "mkdir ",depot//"/rmaps/North_America"
Write(2,*) "mkdir ",depot//"/rmaps/Antarctica"
Write(2,*) "mkdir ",depot//"/stokes/"
Write(2,*) "else"
Write(2,*) "echo"
Write(2,*) "echo '+++> WARNING/ a repository already exists with name: '", depot 
Write(2,*) "    "
Write(2,*) "echo '+++> [existing data will be overwritten]'"
Write(2,*) "    "
Write(2,*) "fi"
!
!
!
! =======================
! --- System settings --- 
! ========================
!
! >>>>>>>> Setting the number of threads for OpenMP execution 
!
if( (option_omp=='y') ) then
  Write(2,*) " "
  Write(2,*) "#echo --------------------------------------------------------"
  Write(2,*) "echo                                                          "
  Write(2,*) " echo '---> Number of threads for OpenMP execution: '", nthread 
  Write(2,*) "#echo --------------------------------------------------------"
  Write(2,*) " echo"
  Write(2,*) "export OMP_NUM_THREADS=",trim(nthread)
Endif
!
! >>>>>>>> Adjusting record size for the Intel Fortran Compiler
!
if( (option_sys=='2') .or. (option_sys=='4') .or. (option_sys=='9') ) then
  Write(2,*) " "
  Write(2,*) "#echo --------------------------------------------------------"
  Write(2,*) "echo                                                          "
  Write(2,*) " echo '---> Setting the FORT_FMT_RECL environment variable'   "
  Write(2,*) "#echo --------------------------------------------------------"
  Write(2,*) " echo"
  Write(2,*) "export FORT_FMT_RECL=1024"
Endif
!
!
!
if (option_npx=='y') then          ! ++++++++++++++++++++ A new pixelization
!
!
!
! ====================================
! EXE 01 --- Pixelizing the sphere --- 
! ====================================
!
!
Write(2,*) " "
Write(2,*) "#echo --------------------------------------------------------"
Write(2,*) "echo                                                          "
Write(2,*) " echo '---> PX.F90: Hicosahedral pixelization of the sphere  '"
Write(2,*) "#echo --------------------------------------------------------"
Write(2,*) "./px.exe"
Write(2,*) "cp px.dat ", depot//"/px"
Write(2,*) "cp pxa.dat ", depot//"/px"
Write(2,*) "cp px-lat.dat ", depot//"/px"
!
!
! ========================================================
! EXE 02 --- Pixelizing the PRESENT-DAY ocean function ---
! ========================================================
!
 Write(2,*) " "
 Write(2,*) "#echo ---------------------------------------"
 Write(2,*) "echo" 
 Write(2,*) " echo '---> px.gmt: Separating wet from dry pixels'"
 Write(2,*) "#echo ---------------------------------------"
!
!--- A realistic Ocean Function
 if(option_rof=='r') then 
	     Write(2,*) "#echo -------------------------------------"
             Write(2,*) " echo '     - Realistic ocean function'"
	     Write(2,*) "#echo -------------------------------------"
             Call OF_pixelization_real
 Endif
! 
!--- A zonal Ocean Function
 if(option_rof=='z') then 
	     Write(2,*) "#echo ---------------------------------"
             Write(2,*) " echo '     - ZONAL ocean function'"
	     Write(2,*) "#echo ---------------------------------"
 	     FILE_CAP='shore.dat'	  
	     Call OF_pixelization_zonal(file_cap, radius_zof) 
 Endif
!
 Write(2,*) "sh px.gmt"
 Write(2,*) "cp weta.dat ", depot//"/px" 
 Write(2,*) "cp drya.dat ", depot//"/px" 
! 
 Write(2,*) " "
 Write(2,*) "#echo ---------------------------------------"
 Write(2,*) "echo" 
 Write(2,*) " echo '---> PX_REC.F90: Merging the wet & dry pixels tables'"
 Write(2,*) "#echo ---------------------------------------" 
 Write(2,*) "./pxrec.exe"
 Write(2,*) "cp px-table.dat ", depot//"/px" 
 Write(2,*) "mv px.gmt ", depot//"/px"
 Write(2,*) "mv -f  px-table.dat ",trim(file_pxtable)
 Write(2,*) "ln -sf ",trim(file_pxtable)," px-table.dat"!
!
!
elseif (option_npx=='n') then        ! ++++++++++++++++++++ An existing pixelization
!
!
 Write(2,*) " "
 Write(2,*) "#echo --------------------------------------------------------"
 Write(2,*) "echo                                                          "
 Write(2,*) " echo '---> PXREBUILD.F90: Retrieving information from pixel table file  '"
 Write(2,*) "#echo --------------------------------------------------------"
 Write(2,*) "ln -sf ",trim(file_pxtable)," px-table.dat"
 Write(2,*) "./pxrebuild.exe"
 Write(2,*) "cp px.dat ", depot//"/px"
 Write(2,*) "cp pxa.dat ", depot//"/px"
 Write(2,*) "cp px-lat.dat ", depot//"/px"
 Write(2,*) "cp weta.dat ", depot//"/px" 
 Write(2,*) "cp drya.dat ", depot//"/px" 
!Write(2,*) "mv px.gmt ", depot//"/px"
!
!
endif
!
!
!
	If(option_px=='y') then 
	Write(2,*) " "
	Write(2,*) "#echo -----------------------------------"
 	Write(2,*) "echo" 	
        Write(2,*) " echo '---> pxmap.gmt: Producing pixelization maps'"
	Write(2,*) "#echo -----------------------------------"
		file_gmt="pxmap.gmt"
		Call make_pxmap (nresolution, file_gmt)
	    if(option_gmt=='y') Write(2,*) "sh pxmap.gmt"
		if(option_gmt=='y') Write(2,*) "ps2pdf px.ps"
		if(option_gmt=='y') Write(2,*) "ps2pdf px-sphere.ps"			
		if(option_gmt=='y') Write(2,*) "mv px.pdf ", depot//"/px"
		if(option_gmt=='y') Write(2,*) "mv px.ps ", depot//"/px"
		if(option_gmt=='y') Write(2,*) "ps2pdf px-sphere.ps"			
		if(option_gmt=='y') Write(2,*) "mv px-sphere.pdf ", depot//"/px"
		if(option_gmt=='y') Write(2,*) "mv px-sphere.ps ", depot//"/px"
		Write(2,*) "mv ", trim(adjustl(file_gmt)), " ", depot//"/px"
	Endif	
!
!
! ==================================================
! EXE 03 --- Computing the spherical harmonics  --- 
! ==================================================
! 
If(option_sh=='y') then 
	Write(2,*) " "
	Write(2,*) "#echo ----------------------------------------"
	Write(2,*) "echo" 
	Write(2,*) " echo '---> SH.F90: Building the spherical harmonics'"
	Write(2,*) "#echo ----------------------------------------"
	Write(2,*) "if [ -e sh.bin ]; then /bin/rm sh.bin ; fi"
	Write(2,*) "if [ -h sh.bin ]; then /bin/rm sh.bin ; fi"
	Write(2,*) "./sh.exe"
	Write(2,*) "mv sh.bin ", trim(adjustl(sh_file)) 
	Write(2,*) "ln -sf ", trim(adjustl(sh_file))," sh.bin" 	
!	Write(2,*) "cp sh.bin ", trim(adjustl(sh_file)) 
!
	           else
	Write(2,*) " "		   
	Write(2,*) "#echo ------------------------------------------------------"
	Write(2,*) "echo" 	
	Write(2,*) " echo '+++> A SH file already exists with name:' ", trim(adjustl(sh_file)) 
	Write(2,*) "#echo ------------------------------------------------------"	
!	Write(2,*) "cp ", trim(adjustl(sh_file)), " sh.bin"
	Write(2,*) "ln -sf ", trim(adjustl(sh_file)), " sh.bin"
!
endif	    
!
!	 
! =========================================================
! EXE 04 ---  Evaluating & plotting the WINDOW function ---
! =========================================================
!
If(option_wi=='y') then
	Write(2,*) " " 
	Write(2,*) "#echo -----------------------------------------------------"
	Write(2,*) "echo" 
	Write(2,*) " echo '--->  WNW.F90: SH orthonormality evaluating the window function'"
	Write(2,*) "#echo -----------------------------------------------------"
	Write(2,*) "./wnw.exe"
	file_gmt="wnw.gmt"
	call make_wnw (resolution, degree, file_gmt)
	if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file_gmt))		
		Write(2,*) "mv ", trim(adjustl(file_gmt)), " "//depot//"/wnw"
		Write(2,*) "mv wnw.dat ", depot//"/wnw"
		if(option_gmt=='y') Write(2,*) "ps2pdf wnw.ps "
		if(option_gmt=='y') Write(2,*) "mv wnw.pdf ", depot//"/wnw"
		if(option_gmt=='y') Write(2,*) "mv wnw.ps ", depot//"/wnw"
		Write(2,*) "mv tmpp4.dat ", depot//"/wnw"
		Write(2,*) "mv tmpp5.dat ", depot//"/wnw"
		Write(2,*) "mv tmpp6.dat ", depot//"/wnw"
		Write(2,*) "mv tmpp7.dat ", depot//"/wnw"
		Write(2,*) "/bin/rm  tmpw0.dat"
		Write(2,*) "/bin/rm  tmpw1.dat"
		Write(2,*) "/bin/rm  tmpw2.dat"
		Write(2,*) "/bin/rm  tmpw3.dat"
Endif
!
!
! ===================================================================
! EXE 05 --- Evaluating the SH coefficients of the PRESENT DAY OF ---
! ===================================================================
!
If(option_oh=='y') then 
!
! A new SH decomposition of the present-day OF is performed 
	Write(2,*) " "
	Write(2,*) "#echo ------------------------------------------"
	Write(2,*) "echo                                         " 
	Write(2,*) " echo '---> SH_OF.F90: SH expansion of the present-day ocean function'"
	Write(2,*) "#echo ------------------------------------------"
	Write(2,*) "./shof.exe"
	Write(2,*) "cp shof.dat ", trim(adjustl(shof_file)) 
	Write(2,*) "cp ", trim(adjustl(shof_file)), " ", depot//"/of"
		   else
!		   
! An existing SH decomposition is employed  
!
	Write(2,*) " "
	Write(2,*) "#echo -----------------------------------------------------------"
	Write(2,*) "echo" 	
	Write(2,*) " echo '+++> An OF SH file already exists with name:' ", trim(adjustl(shof_file))
	Write(2,*) "#echo -----------------------------------------------------------" 
	Write(2,*) "cp ", trim(adjustl(shof_file)), " shof.dat"
	Write(2,*) "cp ", trim(adjustl(shof_file)), " ", depot//"/of"
Endif
!
!
! ==================================================
! EXE 06 --- Present-day OF Degree Variance (DV) ---
! ==================================================
!
If(option_ofdv=='y') then 
	Write(2,*) " "
	Write(2,*) "#echo ---------------------------"
	Write(2,*) "echo"
	Write(2,*) " echo '---> Computing the present-day OF DV'"
	Write(2,*) "#echo ---------------------------"
!
! Computing the DV 
	Write(2,*) "./ofdv.exe"  
!
! Producing a script for plotting the degree variance 
	file_gmt="ofdv.gmt"
 	Call make_ofdvmap (degree, file_gmt) 
	if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file_gmt))
	if(option_gmt=='y') Write(2,*) "ps2pdf ofdv.ps"
	Write(2,*) "mv ofdv.dat ",      depot//"/of/degree_variance"
	if(option_gmt=='y') Write(2,*) "mv ofdv.ps ",      depot//"/of/degree_variance"
	if(option_gmt=='y') Write(2,*) "mv ofdv.pdf ",     depot//"/of/degree_variance"						
	Write(2,*) "mv dv-power.tmp ", depot//"/of/degree_variance"
	Write(2,*) "mv ofdv.gmt ",     depot//"/of/degree_variance"
Endif
!
!
! ============================================================
! EXE 07 --- Reconstruction and mapping the present-day OF ---
! ============================================================
!
If(option_of=='y') then 
	Write(2,*) " "
	Write(2,*) "#echo -----------------------------------------------------"
	Write(2,*) "echo" 
	Write(2,*) " echo '---> REC_OF.F90: Reconstructing and mapping the ocean function'"
	Write(2,*) "#echo -----------------------------------------------------"
!
! Reconstruction 
	Write(2,*) "./recof.exe"  
!		
! Mapping 	
	file_gmt="of.gmt"
 	Call make_ofmap (degree, option_rof, radius_zof, file_cap, file_gmt)    		
	if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file_gmt)) 
		Write(2,*) "mv recof.dat ", depot//"/of"
		if(option_gmt=='y') Write(2,*) "ps2pdf of.ps"
		If(option_rof=='z')& 
		  Write(2,*) "cp ", trim(adjustl(file_cap)), " "//depot//"/of" 
		if(option_gmt=='y') Write(2,*) "mv of.ps ", depot//"/of"
		if(option_gmt=='y') Write(2,*) "mv of.pdf ", depot//"/of"
		Write(2,*) "mv ", trim(adjustl(file_gmt)), " "//depot//"/of"
		If(option_gmt=='y'.and.option_rof=='r')&
		  Write(2,*) "mv pale.cpt ", " "//depot//"/of"		  
Endif 
!
!!
! ===============================================
! --- Copying the ice file from the ice store --- 
! ===============================================
	Write(2,*) " "
        Write(2,*) "#echo --------------------------------------------------------------------"
	Write(2,*) "echo"   
        Write(2,*) " echo '---> Importing '", trim(adjustl(ice_file))," from ICE-MODELS/"   
        Write(2,*) "#echo --------------------------------------------------------------------"  	
	Write(2,*) "cp ./ICE-MODELS/"//trim(adjustl(ice_file)), " ", "./"//trim(adjustl(ice_file))
!
!
! ======================================
! EXE 09 --- Computing shape factors --- 
! ======================================
!
IF(option_sf=='y') THEN 
!
  Write(2,*) " "
  Write(2,*) "#echo ------------------------------------------------------------"  
  Write(2,*) "echo                                                                   " 
  Write(2,*) " echo '---> SHAPE_FACTORS_MOD.F90: Computing the shape factors for model: '", trim(adjustl(ice_file))   
  Write(2,*) "#echo ------------------------------------------------------------"  
  Write(2,*) "./shapefactors.exe"
!
  Write(2,*) " "
  Write(2,*) "#echo -------------------------------------------"
  Write(2,*) "echo                                                                "   
  Write(2,*) " echo '---> SHICE_MOD.F90: Computing SH coefficients for the ice model'"  
  Write(2,*) "#echo -------------------------------------------"  
  Write(2,*) "./shice.exe" 
!
  Write(2,*) "cp shice.dat ", depot//"/"//trim(adjustl(titlice))//"/sh/"
!
  Write(2,*) "cp shice.dat ", trim(adjustl(shape_file))
!  
  If(option_ri=='n') Write(2,*) "mv sht*.dat ", depot//"/"//trim(adjustl(titlice))//"/sh/"
!
  		   ELSE
!
   Write(2,*) " "
   Write(2,*) "#echo ------------------------------------------------------------" 
   Write(2,*) "echo                                                                    "  
   Write(2,*) " echo '---> Ice harmonics are pre-computed in file: '", shape_file
   Write(2,*) "#echo ------------------------------------------------------------"  
!
   Write(2,*) "cp ", trim(adjustl(shape_file)), " shice.dat"
!
  If(option_ri=='y') then 
   Write(2,*) "#echo ------------------------------------------------------------" 
   Write(2,*) "echo                                                              "  
   Write(2,*) " echo '---> Copying from '", depot, " the ice harmonics for reconstruction"
   Write(2,*) "#echo ------------------------------------------------------------"  
!  
  Write(2,*) "cp  ", depot//"/"//trim(adjustl(titlice))//"/sh/"//"sht*.dat  ."
!
  Endif	
!
  ENDIF 
!
!
! Old piece of script: 
!
! IF(option_sf=='y') THEN 
!
!  Write(2,*) " "
!  Write(2,*) "#echo ------------------------------------------------------------"  
!  Write(2,*) "echo                                                                   " 
!  Write(2,*) " echo '---> SHAPE_FACTORS.F90: Computing the shape factors for model: '", trim(adjustl(ice_file))   
!  Write(2,*) "#echo ------------------------------------------------------------"  
!
!  Write(2,*) "./shapefactors.exe"
!
!  Write(2,*) "cp shicec.bin ", trim(adjustl(shape_file))   
!
!  ELSE
!
!   Write(2,*) " "
!   Write(2,*) "#echo ------------------------------------------------------------" 
!   Write(2,*) "echo                                                                    "  
!   Write(2,*) " echo '---> SHAPE_FACTORS.F90: Shape factors are pre-computed in file: '", shape_file
!   Write(2,*) "#echo ------------------------------------------------------------"  
!   Write(2,*) "cp ", trim(adjustl(shape_file)), " shicec.bin"
!
!  ENDIF 
!
!
!
!
!
! ============================================================
! EXE 10 --- Computing the SH coefficients from the shape factors --- 
! ============================================================
!
!       Write(2,*) " "
!       Write(2,*) "#echo -------------------------------------------"
! 	Write(2,*) "echo                                                                "   
!   	Write(2,*) " echo '---> SHICE.F90: Computing SH coefficients for the ice model'"  
!       Write(2,*) "#echo -------------------------------------------"  
!   	Write(2,*) "./shice.exe" 
!
!        Write(2,*) "cp shice.dat ", depot//"/"//trim(adjustl(titlice))//"/sh/"
!
!	If(option_ri=='n') Write(2,*) "mv sht*.dat ", depot//"/"//trim(adjustl(titlice))//"/sh/"
!
!
! ================================================== 
! EXE 12 --- Mapping of the original ice sheets ---
! ==================================================
!
If(option_or=='y') then 
!
Write(2,*) " "
Write(2,*) "#echo ---------------------------------------------" 
Write(2,*) "echo"    
Write(2,*) "echo '---> MS.F90: Creating multi-segment files for ice sheets maps'"
Write(2,*) "#echo ---------------------------------------------"  
Write(2,*) "./ms.exe" 
!
!
Write(2,*) " "
Write(2,*) "#echo ------------------------------------------------"  
Write(2,*) "echo"   
Write(2,*) " echo '---> mapice.gmt: Creating ps images of original ice sheets'"
Write(2,*) "#echo ------------------------------------------------"  	   	
!	   	
!--- Creates a "rainbow palette" suitable for plotting the 
!    reconstructed ice distribution, named "ice-pal.cpt": the 
!    same palette is also used for the reconstructed ice 
!    distribution...(see below)
!
!
   file_gmt="mapice.gmt"
   Call make_icemap (ninc, titlice, option_rof, FILE_CAP, file_gmt)
!		
   if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file_gmt))
   Write(2,*) "mv msg* ", depot//"/"//trim(adjustl(titlice))//"/original/"
   Write(2,*) "mv rtmp*.dat ", depot//"/"//trim(adjustl(titlice))//"/original/"   
   if(option_gmt=='y') Write(2,*) "mv mapice*.ps ", depot//"/"//trim(adjustl(titlice))//"/original/"
   if(option_gmt=='y') Write(2,*) "mv mapice*.pdf ", depot//"/"//trim(adjustl(titlice))//"/original/"
   if(option_gmt=='y') Write(2,*) "mv pice*.cpt ", depot//"/"//trim(adjustl(titlice))//"/original/"
   Write(2,*) "mv ",  trim(adjustl(file_gmt)), & 
   	       " "//depot//"/"//trim(adjustl(titlice))//"/original/"
   If(option_rof=='z')& 
      Write(2,*)"cp ", trim(adjustl(file_cap)), " ", & 
      depot//"/"//trim(adjustl(titlice))//"/original/"
ENDIF
!
!
! ============================================================
! EXE 13---  Reconstruction and mapping of the ice sheets ---
! ============================================================
!
If(option_ri=='y') then 
!
Write(2,*) " "
Write(2,*) "#echo -------------------------------------------------"  
Write(2,*) "echo"   
Write(2,*) " echo '---> REC_ICE.F90: Reconstruction of ice sheets distribution'"
Write(2,*) "#echo -------------------------------------------------"  
!
! --- Reconstruction ...
!
   Write(2,*) "./recice.exe"
!	 
! --- ... and mapping of reconstruction 
! 
Write(2,*) " "
Write(2,*) "#echo ---------------------------------------------------------"  
Write(2,*) "echo"   
Write(2,*) " echo '---> recice.gmt: Creating ps images of reconstructed ice sheets'"
Write(2,*) "#echo ---------------------------------------------------------"  
!

!	   	
!--- Creates a "rainbow palette" suitable for plotting the 
!    reconstructed ice distribution, named "ice-pal.cpt": the 
!    same palette is also used for the original ice sheets
!    distribution...(see above). 
!
!
   File_gmt="recice.gmt"	
   Call make_recicemap (ninc, degree, titlice, option_rof, file_cap, file_gmt)
!
   if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file_gmt))
   Write(2,*) "mv rect*.dat ", depot//"/"//trim(adjustl(titlice))//"/reconstructed/"
   Write(2,*) "mv tmp0*.dat ", depot//"/"//trim(adjustl(titlice))//"/reconstructed/"
   if(option_gmt=='y') Write(2,*) "mv recice*.ps ", depot//"/"//trim(adjustl(titlice))//"/reconstructed/"
   if(option_gmt=='y') Write(2,*) "mv recice*.pdf ", depot//"/"//trim(adjustl(titlice))//"/reconstructed/"
   if(option_gmt=='y') Write(2,*) "mv pice*.cpt ", depot//"/"//trim(adjustl(titlice))//"/reconstructed/"
   Write(2,*) "mv ", trim(adjustl(file_gmt)), & 
              " "//depot//"/"//trim(adjustl(titlice))//"/reconstructed/" 
!
   Write(2,*) "mv sht*.dat ", depot//"/"//trim(adjustl(titlice))//"/sh/"            
!
ENDIF
!
!
!
!
! =========================================
! EXE 14 --- Computing the Love numbers ---
! =========================================
! 
!************************
IF(OPTION_NM=='y') THEN 
!************************
!
Write(2,*) "#echo  ---------------------------------------------------------------------" 
Write(2,*) "echo " 
Write(2,*) " echo  '---> Load-deformation coefficients by TABOO (Normal Modes)'"
Write(2,*) "#echo  ---------------------------------------------------------------------"  
Write(2,*) "./tb.exe"         
!
if(option_tidal=='y') then 
!
	Write(2,*) "mv h_hp.dat h.dat"
	Write(2,*) "mv k_hp.dat k.dat"
	Write(2,*) "mv l_hp.dat l.dat" 
!
	Write(2,*) "mv hh_hp.dat hh.dat" 
	Write(2,*) "mv kk_hp.dat kk.dat"
	Write(2,*) "mv ll_hp.dat ll.dat"
!
	Write(2,*) "mv ih_hp.dat ih.dat"
	Write(2,*) "mv ik_hp.dat ik.dat"
	Write(2,*) "mv il_hp.dat il.dat"
!
	Write(2,*) "mv ihh_hp.dat ihh.dat"
	Write(2,*) "mv ikk_hp.dat ikk.dat"
	Write(2,*) "mv ill_hp.dat ill.dat"
!
	Write(2,*) "mv h_tidal_hp.dat h_tidal.dat"
	Write(2,*) "mv k_tidal_hp.dat k_tidal.dat"
	Write(2,*) "mv l_tidal_hp.dat l_tidal.dat"
!
	Write(2,*) "mv ih_tidal_hp.dat ih_tidal.dat"
	Write(2,*) "mv ik_tidal_hp.dat ik_tidal.dat"
	Write(2,*) "mv il_tidal_hp.dat il_tidal.dat"
!
	Write(2,*) "mv spectrum_hp.dat spectrum.dat"
	Write(2,*) "mv ss_hp.dat ss.dat"
!
endif
!
!
Write(2,*) "echo " 
Write(2,*) "echo '+++> WARNING/ in SELEN 3, the average density of the Earth is computed by the'"
Write(2,*) "echo '            / density structure of the input model, NOT using an a-priori'"
Write(2,*) "echo '            / value as done in previous versions - GS & FC July 27 2009 - '"
!
!
Write(2,*) "/bin/rm visco.tmp"
!
Write(2,*) "cp ", trim(adjustl(visco_file)), " ", depot//"/TABOO/"
Write(2,*) "cp ", trim(adjustl(visco_file)), " ", depot//"/Love-Numbers-by-TABOO/"  
!
Write(2,*) "cp task_1.dat ", depot//"/log/"
Write(2,*) "mv task_1.dat ", depot//"/TABOO/"
!
Write(2,*) "mv spectrum.dat ", depot//"/Love-Numbers-by-TABOO/"
Write(2,*) "mv ih.dat ",       depot//"/Love-Numbers-by-TABOO/"
Write(2,*) "mv ik.dat ",       depot//"/Love-Numbers-by-TABOO/"
Write(2,*) "mv il.dat ",       depot//"/Love-Numbers-by-TABOO/"
Write(2,*) "mv h.dat ",        depot//"/Love-Numbers-by-TABOO"
Write(2,*) "mv l.dat ",        depot//"/Love-Numbers-by-TABOO"
Write(2,*) "cp k.dat ",        depot//"/Love-Numbers-by-TABOO"     ! Changed (it was "mv")  
!
if(option_tidal=='y') then 
	Write(2,*) "cp task_1_tidal.dat ",    depot//"/log/"
	Write(2,*) "mv task_1_tidal.dat ",    depot//"/TABOO/"
	Write(2,*) "mv ih_tidal.dat ", depot//"/Love-Numbers-by-TABOO/"
	Write(2,*) "mv ik_tidal.dat ", depot//"/Love-Numbers-by-TABOO/"
	Write(2,*) "mv il_tidal.dat ", depot//"/Love-Numbers-by-TABOO/"
	Write(2,*) "mv h_tidal.dat ",  depot//"/Love-Numbers-by-TABOO"
	Write(2,*) "mv l_tidal.dat ",  depot//"/Love-Numbers-by-TABOO"
	Write(2,*) "mv k_tidal.dat ",  depot//"/Love-Numbers-by-TABOO"
endif
!
If(option_ln=='y') then 
Write(2,*) " "
Write(2,*) "#echo -----------------------------------------------------------"  
Write(2,*) "echo "
Write(2,*) " echo '---> ldcs.gmt: Plotting Love numbers and other spectral quantities'"
Write(2,*) "#echo -----------------------------------------------------------"  
!
	File_gmt="ldcs.gmt"
	CALL make_plot_ldc (nv, code, degree, vstring, file_gmt)
!
    if(option_gmt=='y') Write(2,*) "sh ", file_gmt
	Write(2,*) "mv ", trim(adjustl(file_gmt)), " "//depot//"/Love-Numbers-by-TABOO/"
	Write(2,*) "mv ss.dat ", depot//"/Love-Numbers-by-TABOO/"
	Write(2,*) "mv ihh.dat ", depot//"/Love-Numbers-by-TABOO/"
	Write(2,*) "mv ikk.dat ", depot//"/Love-Numbers-by-TABOO/"
	Write(2,*) "mv ill.dat ", depot//"/Love-Numbers-by-TABOO/"
	Write(2,*) "mv hh.dat ", depot//"/Love-Numbers-by-TABOO/"
	Write(2,*) "mv kk.dat ", depot//"/Love-Numbers-by-TABOO/"
	Write(2,*) "mv ll.dat ", depot//"/Love-Numbers-by-TABOO/"
	if(option_gmt=='y') Write(2,*) "mv h.tmp ", depot//"/Love-Numbers-by-TABOO/"
	if(option_gmt=='y') Write(2,*) "mv k.tmp ", depot//"/Love-Numbers-by-TABOO/"
	if(option_gmt=='y') Write(2,*) "mv l.tmp ", depot//"/Love-Numbers-by-TABOO/"
	if(option_gmt=='y') Write(2,*) "mv spe.tmp ", depot//"/Love-Numbers-by-TABOO/"
	Write(2,*) "mv tmpg*.dat ", depot//"/Love-Numbers-by-TABOO"
	if(option_gmt=='y') Write(2,*) "ps2pdf ela-flu.ps"
	if(option_gmt=='y') Write(2,*) "ps2pdf spectrum.ps"
	if(option_gmt=='y') Write(2,*) "ps2pdf n-residues.ps"
	if(option_gmt=='y') Write(2,*) "mv ela-flu.ps ", depot//"/Love-Numbers-by-TABOO"
	if(option_gmt=='y') Write(2,*) "mv spectrum.ps ", depot//"/Love-Numbers-by-TABOO"
	if(option_gmt=='y') Write(2,*) "mv n-residues.ps ", depot//"/Love-Numbers-by-TABOO"
	if(option_gmt=='y') Write(2,*) "mv ela-flu.pdf ", depot//"/Love-Numbers-by-TABOO"
	if(option_gmt=='y') Write(2,*) "mv spectrum.pdf ", depot//"/Love-Numbers-by-TABOO"
	if(option_gmt=='y') Write(2,*) "mv n-residues.pdf ", depot//"/Love-Numbers-by-TABOO"
Endif
!
!****************************
ELSEIF(OPTION_PW=='y') THEN 
!****************************
Write(2,*) "#echo  --------------------------------------------------------------------" 
Write(2,*) "echo " 
Write(2,*) " echo  '---> alma.f90: Computing Love numbers by the PW formula'"
Write(2,*) "#echo  --------------------------------------------------------------------" 
!
Write(2,*) "echo " 
Write(2,*) "echo '+++> WARNING/ in SELEN 3, the average density of the Earth is computed by the'"
Write(2,*) "echo '            / density structure of the input model, NOT using an a-priori'"
Write(2,*) "echo '            / value as done in previous versions - GS & FC July 27 2009 - '"
Write(2,*) "echo " 
!
Write(2,*) "echo '     +++++++++++++++++++++++++++++++++++++++++++++++'" 
Write(2,*) "echo '                        A L M A                     '" 
Write(2,*) "echo '         Version 1.1 - last modified July 2009      '"
Write(2,*) "echo '                - ported to SELEN 3.2 -             '"
Write(2,*) "echo '                  Copyright  (C) 2007               '" 
Write(2,*) "echo '      Giorgio Spada -Email: giorgio.spada@gmail.com '"
Write(2,*) "echo '     +++++++++++++++++++++++++++++++++++++++++++++++'"
Write(2,*) trim(adjustl(wdir))//"/ALMA/alma.exe"
!
Write(2,*) "echo '     +++++++++++++++++++++++++++++++++++++++++++++++'" 
Write(2,*) "echo '                        A L M A                     '" 
Write(2,*) "echo '         Version 1.1 - last modified July 2009      '"
Write(2,*) "echo '                - ported to SELEN 3.2 -             '"
Write(2,*) "echo '                  Copyright  (C) 2007               '" 
Write(2,*) "echo '      Giorgio Spada -Email: giorgio.spada@gmail.com '"
Write(2,*) "echo '     +++++++++++++++++++++++++++++++++++++++++++++++'"
Write(2,*) "cp ./ALMA/alma-logfile.dat ." 
Write(2,*) "cp ./ALMA/h.dat hpw.dat"
Write(2,*) "cp ./ALMA/l.dat lpw.dat"
Write(2,*) "cp ./ALMA/k.dat kpw.dat"
Write(2,*) "echo " 
Write(2,*) " echo  '---> AM.F90: Building the Green functions by the Love numbers '"
Write(2,*) "./am.exe"
Write(2,*) "mv hpw.dat ", depot//"/Love-Numbers-by-ALMA/"
Write(2,*) "mv lpw.dat ", depot//"/Love-Numbers-by-ALMA/"
Write(2,*) "mv kpw.dat ", depot//"/Love-Numbers-by-ALMA/"
Write(2,*) "cp ", trim(adjustl(wdir))//"/ALMA/alma.inc ",  depot//"/Love-Numbers-by-ALMA/alma.inc" 
Write(2,*) "cp ", trim(adjustl(wdir))//"/ALMA/alma.inc ",  depot//"/ALMA/alma.inc" 
!
!****************************
ELSEIF(OPTION_PWA=='y') THEN 
!****************************
!
Write(2,*) "#echo  --------------------------------------------------------------------" 
Write(2,*) "echo " 
Write(2,*) " echo  '---> alma.f90: Using pre-computed Love numbers by ALMA '"
Write(2,*) "#echo  --------------------------------------------------------------------" 
!
Write(2,*) "echo '     +++++++++++++++++++++++++++++++++++++++++++++++'" 
Write(2,*) "echo '                        A L M A                     '" 
Write(2,*) "echo '         Version 1.1 - last modified July 2009      '"
Write(2,*) "echo '                - ported to SELEN 3.2 -             '"
Write(2,*) "echo '                  Copyright  (C) 2007               '" 
Write(2,*) "echo '      Giorgio Spada -Email: giorgio.spada@gmail.com '"
Write(2,*) "echo '     +++++++++++++++++++++++++++++++++++++++++++++++'"
!
Write(2,*) "cp ./LOVE_DEPOSIT/alma-logfile_"//trim(adjustl(short_visco_filename(12:100))), " ./alma-logfile.dat" 
Write(2,*) "cp ", trim(adjustl(alma_file_pwa_h)),  " hpw.dat"
Write(2,*) "cp ", trim(adjustl(alma_file_pwa_l)),  " lpw.dat"
Write(2,*) "cp ", trim(adjustl(alma_file_pwa_k)),  " kpw.dat"
Write(2,*) "./am.exe"
Write(2,*) "mv hpw.dat ", depot//"/Love-Numbers-by-ALMA/"
Write(2,*) "mv lpw.dat ", depot//"/Love-Numbers-by-ALMA/"
Write(2,*) "mv kpw.dat ", depot//"/Love-Numbers-by-ALMA/"
!
!******
ENDIF 
!******
!
!
!
! ================================================================================
! --- Determining the spectral form of the Polar Motion Transfer function PMTF ---
! ================================================================================
!
If(option_pmtf=='y'.and.option_tlove_ext=='n') then 
Write(2,*) " "
Write(2,*) "#echo  --------------------------------------------------------------------- " 
Write(2,*) "echo " 
Write(2,*) " echo  '---> PM_v6.F90: Determining the Polar Motion Transfer Function, PMTF'"
Write(2,*) "#echo  --------------------------------------------------------------------- "  
Write(2,*) "./pmtf.exe" 
Write(2,*) "cp PM_spectrum.dat        ", depot//"/PMTF/"	
Write(2,*) "cp PM_time_domain.dat     ", depot//"/PMTF/"	
Write(2,*) "cp Deg_2_Love_numbers.dat ", depot//"/PMTF/"	
Endif
!
!
!
! ============================================================
! --- Computing the direct effect of load on polar motion  ---
! ============================================================
!
If(option_pmd=='y') then 
!
 Write(2,*) " "
 Write(2,*) "#echo  ---------------------------------------------------------- " 
 Write(2,*) "echo " 
 Write(2,*) " echo  '---> LOAD_RE.F90: Total load (I + S) for a rigid Earth'"
 Write(2,*) "#echo  ---------------------------------------------------------- "  
!
 Write(2,*) "./loadre.exe" 
 Write(2,*) "cp shload_re.bin   ", depot//"/PM-direct-effect/"	
 Write(2,*) "cp shload_rice.bin ", depot//"/PM-direct-effect/"	
 Write(2,*) "cp shload_roce.bin ", depot//"/PM-direct-effect/"	
!

 if(option_tlove_ext=='n') then 

!
 Write(2,*) " "
 Write(2,*) "#echo  ---------------------------------------------------------- " 
 Write(2,*) "echo " 
 Write(2,*) " echo  '---> PMD_MOD.F90: Polar motion driven by the surface load'"
 Write(2,*) "#echo  ---------------------------------------------------------- "  
! 
! 
 Write(2,*) "./pmd.exe" 
 Write(2,*) "cp m.dat     ", depot//"/PM-direct-effect/"	
 Write(2,*) "cp m.dot     ", depot//"/PM-direct-effect/"	
 Write(2,*) "cp Rotational_data.bck ", depot//"/PM-direct-effect/"
!
!
 If(option_pmplot=='y')then 
!
 file_gmt='pm.gmt'
 Call MAKE_PM (RUN, & 
               NV,  & 
	       CODE, & 
	       TITLICE, & 
	       RESOLUTION, & 
	       ITER, & 
               MODE, & 
	       DEGREE, & 
	       VSTRING, & 
	       SHORT_VISCO_FILENAME, & 
               NINC,     & 
	       FILE_GMT)
!
 if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file_gmt))	
 if(option_gmt=='y') Write(2,*) "ps2pdf m.ps     " 
 if(option_gmt=='y') Write(2,*) "ps2pdf m-dot.ps " 	
!
 Write(2,*) "cp ", trim(adjustl(file_gmt)), " ", depot//"/PM-direct-effect/"	
 if(option_gmt=='y') Write(2,*) "cp m.ps      ", depot//"/PM-direct-effect/"	
 if(option_gmt=='y') Write(2,*) "cp m-dot.ps  ", depot//"/PM-direct-effect/"	
 if(option_gmt=='y') Write(2,*) "cp m.pdf     ", depot//"/PM-direct-effect/"	
 if(option_gmt=='y') Write(2,*) "cp m-dot.pdf ", depot//"/PM-direct-effect/"	
 Write(2,*) "cp mx.tmp     my.tmp     mod.tmp         ", depot//"/PM-direct-effect/"	
 Write(2,*) "cp mdot-x.tmp mdot-y.tmp mdot-mod.tmp    ", depot//"/PM-direct-effect/"	
 Write(2,*) "cp pm-1.tmp   pm-2.tmp                   ", depot//"/PM-direct-effect/"	
!
Endif
!
Endif 
!
ENDIF
!     
!
!
!
! ================================================================================
! --- Solving the Sea Level Equation (SLE) in the USUAL way (FIXED coastlines) ---
! ================================================================================
!
If(OPTION_TOPO=='n'.and.OPTION_RFB=='n') then 
	Write(2,*) " "
	Write(2,*) "#echo  ------------------------------------------------------------------- " 
	Write(2,*) "echo " 
	Write(2,*) " echo  '---> SLE.F90: Solving the SLE - for FIXED coastlines with NO rotational feedback'"
	Write(2,*) "#echo  ------------------------------------------------------------------- "  
	Write(2,*) "./sle.exe" 
Endif 
!
If(OPTION_TOPO=='n'.and.OPTION_RFB=='y') then 
	If(option_tlove_ext=='n') then 
! Write(2,*) " "
! Write(2,*) "#echo  ------------------------------------------------------------------- " 
! Write(2,*) "echo " 
! Write(2,*) " echo  '---> S_ROT.F90: Amplitude of the degree 2 order 1 polar motion sea level tide'"
! Write(2,*) "#echo  ------------------------------------------------------------------- "  	
! Write(2,*) "./s_rot.exe"
! Write(2,*) " "
! Write(2,*) "#echo  ------------------------------------------------------------------- " 
! Write(2,*) "echo " 
! Write(2,*) " echo  '---> U_ROT.F90: Amplitude of the degree 2 order 1 rotational displacement'"
! Write(2,*) "#echo  ------------------------------------------------------------------- "  	
! Write(2,*) "./u_rot.exe"
! Write(2,*) " "
! Write(2,*) "#echo  ------------------------------------------------------------------- " 
! Write(2,*) "echo " 
! Write(2,*) " echo  '---> N_ROT.F90: Amplitude of the degree 2 order 1 sea surface polar motion tide'"
! Write(2,*) "#echo  ------------------------------------------------------------------- "  	
! Write(2,*) "./n_rot.exe"

	Endif 
	If(option_tlove_ext=='y') then 
!
!	Write(2,*) " "
!	Write(2,*) "#echo  ------------------------------------------------------------------- " 
!	Write(2,*) "echo " 
!	Write(2,*) " echo  '---> S_ROT_PM.F90: Amplitude of the degree 2 order 1 polar motion sea level tide'"
!	Write(2,*) "#echo  ------------------------------------------------------------------- "  	
!	Write(2,*) "./s_rot_pm.exe"
!	Write(2,*) " "
! Write(2,*) "#echo  ------------------------------------------------------------------- " 
! Write(2,*) "echo " 
! Write(2,*) " echo  '---> U_ROT_PM.F90: Amplitude of the degree 2 order 1 rotational displacement'"
! Write(2,*) "#echo  ------------------------------------------------------------------- "  	
! Write(2,*) "./u_rot_pm.exe"
! Write(2,*) " "
! Write(2,*) "#echo  ------------------------------------------------------------------- " 
! Write(2,*) "echo " 
! Write(2,*) " echo  '---> N_ROT_PM.F90: Amplitude of the degree 2 order 1 sea surface polar motion tide'"
! Write(2,*) "#echo  ------------------------------------------------------------------- "  	
! Write(2,*) "./n_rot_pm.exe"
!
	Endif 

!
	Write(2,*) "#echo  ------------------------------------------------------------------- " 
	Write(2,*) "echo " 
	Write(2,*) " echo  '---> SLE.F90: Solving the SLE - for FIXED coastlines WITH rotational feedback'"
	Write(2,*) "#echo  ------------------------------------------------------------------- "  
	Write(2,*) "./sle.exe" 
Endif 
!


!

!
! ===================================================================
! --- Solving the Sea Level Equation (SLE) for VARYING coastlines ---  
! ===================================================================
!
If(OPTION_TOPO=='y') then 
!
! --------------------------------------------------
!  0) Preliminaries (topo pixelization & ice masking) 
! --------------------------------------------------
!
If(option_pxtopo=='y') then 
Write(2,*) " " 
Write(2,*) "#echo ------------------------------------------"
Write(2,*) "echo "
Write(2,*) " echo '---> PX_TOPO.F90: Pixelizing the present-day topography'"
Write(2,*) "#echo ------------------------------------------"
Write(2,*) "./pxtopo.exe" 
Else 
Write(2,*) " " 
Write(2,*) "#echo ------------------------------------------"
Write(2,*) "echo "
Write(2,*) " echo '---> Present-day pixelized topography is read from file '", & 
trim(adjustl(FILE_PXTOPO))
Write(2,*) "#echo ------------------------------------------"
Endif 
!
Write(2,*) " " 
Write(2,*) "#echo ------------------------"
Write(2,*) "echo "
Write(2,*) " echo '---> ICE_MASK.F90: Producing Ice masks'"
Write(2,*) "#echo ------------------------"
Write(2,*) "./icemask.exe" 
!
!
!
! -----------------------------------------
!  1) Initialization (present day OF & SLE) 
! -----------------------------------------
!
!   The present-day OF is already decomposed on the SH basis ... 
!
Write(2,*) " "
Write(2,*) "#echo  --------------------------------------------------------- "  
Write(2,*) "echo "
Write(2,*) " echo  '---> SLE.F90: Solving the SLE with present-day ocean function'" 
Write(2,*) "#echo  --------------------------------------------------------- "  
Write(2,*) "./sle.exe" 

!
! -----------------------------------------
!  2) Iterating the SLE  
! -----------------------------------------
!
Write(2,*) " "
Write(2,*) "#echo  ----------------------------------------------------- " 
Write(2,*) "echo " 
Write(2,*) " echo '---> Iterating consistently with paleo-topography'" 
Write(2,*) "#echo  ----------------------------------------------------- " 
! 
for_argument = "((i=1; i<="//ITER_C//";i+=1))"
Write(2,*     ) "for "//trim(adjustl(for_argument))
Write(2,'(a2)') "do"
Write(2,*) " echo"
Write(2,*) " echo '===> Performing iteration ' $i ' of' ", trim(adjustl(iter_c))
       Write(2,*) " echo"  
       Write(2,*) " echo  '---> RSL_CPX.F90: Building the paleo-topography by the history of RSL'" 
Write(2,*) "./rslcpx.exe" 
       Write(2,*) " echo"
       Write(2,*) " echo  '---> SH_OF_VAROC.F90: SH analysis of the paleo OFs'" 
Write(2,*) "./shofvaroc.exe"
       Write(2,*) " echo"
       if( option_rfb=='n' ) then
         Write(2,*) " echo  '---> SLE_VAROC.F90: Solving the SLE with varying coastlines'" 
       else
         Write(2,*) " echo  '---> SLE_VAROC_ROTAZ.F90: Solving the SLE with varying coastlines AND rotational feedback'" 
       endif
Write(2,*) "./slevaroc.exe"
Write(2,'(a4)') "done"
Write(2,*) " echo"
Write(2,*) " echo '===> End of iterations!'"
!
ENDIF  
!
!
!
!

!
! =======================================================
! EXE 08 --- Reconstruction and mapping the PALEO OF ---
! =======================================================
!
If(option_ptmap=='y') then 
	Write(2,*) " "
	Write(2,*) "#echo -----------------------------------------------------"
	Write(2,*) "echo" 
	Write(2,*) " echo '---> paleo_topo.gmt: Mapping the paleo-topography'"
	Write(2,*) "#echo -----------------------------------------------------"
!
! Mapping ...
!	
	file_gmt="paleo_topo.gmt"
 	Call make_paleo_topo_map (ninc,file_gmt)    		
	if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file_gmt)) 
!
! Moving files in the repository  ...
!	
        Write(2,*) "mv ",  trim(adjustl(file_gmt)), " ", " "//depot//"/paleo-topography/data"

   	Write(2,*) "mv ptmp-* ",         " "//depot//"/paleo-topography/data"
        if(option_gmt=='y') Write(2,*) "mv map-topo-*.ps ",  " "//depot//"/paleo-topography/ps"
        if(option_gmt=='y') Write(2,*) "mv map-topo-*.pdf ", " "//depot//"/paleo-topography/pdf"  
        if(option_gmt=='y') Write(2,*) "mv pale_ice.cpt ",    " "//depot//"/paleo-topography/data"
        if(option_gmt=='y') Write(2,*) "mv pale_topo.cpt ",    " "//depot//"/paleo-topography/data"    	
Endif 


!
! ========================================================= ==================
! EXE 11 --- Plotting Equivalent Sea Level (ESL) curves for "ho" ice models---
! ========================================================= ==================
!
If(option_esl=='y') then 
!
    if        (option_topo=='n') then 
    Write(2,*) " "
    Write(2,*) "#echo -------------------------------"  
    Write(2,*) "echo "   
    Write(2,*) " echo '---> ESL.F90: Plotting the ESL curve'"
    Write(2,*) "#echo -------------------------------"  
    Write(2,*) "./esl.exe" 
    elseif    (option_topo=='y') then 
    Write(2,*) " "
    Write(2,*) "#echo -------------------------------"  
    Write(2,*) "echo "   
    Write(2,*) " echo '---> ESL.F90: Plotting the ESL curve (assuming present-day shorelines)'"
    Write(2,*) "#echo -------------------------------"      
    Write(2,*) "./esl.exe" 
    endif  
!
   file_gmt='eslplot.gmt'
   Call make_eslplot (ninc, titlice, shof_file, file_gmt)
!
   if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file_gmt))
   if(option_gmt=='y') Write(2,*) "ps2pdf esl.ps"
   Write(2,*) "mv ", trim(adjustl(file_gmt)), " "//depot//"/"//trim(adjustl(titlice))//"/esl/"
   if(option_gmt=='y') Write(2,*) "mv esl.ps ", " "//depot//"/"//trim(adjustl(titlice))//"/esl/"
   if(option_gmt=='y') Write(2,*) "mv esl.pdf ", " "//depot//"/"//trim(adjustl(titlice))//"/esl/"
   Write(2,*) "mv esl.dat ", " "//depot//"/"//trim(adjustl(titlice))//"/esl/"
   if(option_gmt=='y') Write(2,*) "mv esl.tmp ", " "//depot//"/"//trim(adjustl(titlice))//"/esl/"
   Write(2,*) "mv esl-thin.dat ", " "//depot//"/"//trim(adjustl(titlice))//"/esl/"
   Write(2,*) "mv esl-tot.dat ", " "//depot//"/"//trim(adjustl(titlice))//"/esl/"
!
Endif 
!
!
!
! ========================================================= ==================
! EXE 11b --- Plotting Equivalent Sea Level (ESL) curves for "pm" ice models--
! ========================================================= ==================
!
!
If(option_esl_pm=='y') then 
!
    Write(2,*) " "
    Write(2,*) "#echo -------------------------------"  
    Write(2,*) "echo "   
    Write(2,*) " echo '---> ESL.F90: Plotting the ESL curve for <<pm>> ice'"
    Write(2,*) "#echo -------------------------------"  
    Write(2,*) "./esl_pm.exe" 
!
   file_gmt='eslplot_pm.gmt'
   Call make_eslplot_pm (titlice,file_gmt)
!
   if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file_gmt))
   if(option_gmt=='y') Write(2,*) "ps2pdf esl-pm.ps"
   Write(2,*) "mv ", trim(adjustl(file_gmt)), " "//depot//"/"//trim(adjustl(titlice))//"/esl-pm/"
   if(option_gmt=='y') Write(2,*) "mv esl-pm.ps ", " "//depot//"/"//trim(adjustl(titlice))//"/esl-pm/"
   Write(2,*) "mv esl-pm.dat ", " "//depot//"/"//trim(adjustl(titlice))//"/esl-pm/"
   Write(2,*) "mv esl-pm-thin.dat ", " "//depot//"/"//trim(adjustl(titlice))//"/esl-pm/"
!
Endif 

















!
!
! ========================================
! EXE 12 ---  Relative Sea Level (RSL) ---
! ========================================
!
! RSL database 
If(trim(adjustl(option_rsla))=='y') then
	Write(2,*) " "
	Write(2,*) "#echo ------------------------------------------"
	Write(2,*) "echo"
	Write(2,*) " echo '---> The RSL database is: '", trim(adjustl(rsl_file)) 
	Write(2,*) "#echo ------------------------------------------"
ENDIF
!
! Mapping the distribution of RSL sites ... 
	If(option_rsldb=='y') then 
		Write(2,*) " "
		Write(2,*) "#echo ------------------------"
		Write(2,*) "echo"
		Write(2,*) " echo '---> rsl-s.gmt: Map of RSL sites'"
		Write(2,*) "#echo ------------------------"	
		file_gmt='rsl-s.gmt'
		call MAKE_RSLDB(RSL_FILE, RSL_DATABASE_FORMAT, NRSL, FILE_GMT)
		if(option_gmt=='y') write(2,*) "sh ", file_gmt 
		Write(2,*) "cp ", trim(adjustl(rsl_file)), " "//depot//"/rsl/rsl-sites"		
		Write(2,*) "mv rsl-s.gmt ", " "//depot//"/rsl/rsl-sites"
		if(option_gmt=='y') Write(2,*) "ps2pdf maprsl.ps"
		if(option_gmt=='y') Write(2,*) "mv maprsl.ps ", " "//depot//"/rsl/rsl-sites"
		if(option_gmt=='y') Write(2,*) "mv maprsl.pdf ", " "//depot//"/rsl/rsl-sites"
		Write(2,*) "mv lon-lat-rsl.dat ", " "//depot//"/rsl/rsl-sites"
		Write(2,*) "mv tmptitle ", " "//depot//"/rsl/rsl-sites"	
	ENDIF
!
! --- Computing synthetic RSL curves ... 
	If(option_rsl=='y') then 
		Write(2,*) " "
		Write(2,*) "#echo  ---------------------------------------------------------"
		Write(2,*) "echo"	
		Write(2,*) " echo '---> RSL.F90: Predicting RSL at the sites of database: '", trim(rsl_file)
		Write(2,*) "#echo  ---------------------------------------------------------"	
	        Write(2,*) "./shrsl.exe"
		Write(2,*) "./rsl.exe" 		
!Write(2,*) "mv shrsl.bin ", depot//"/bin"
!
! --- ... and drawing RSL figures for each site		
		If(option_rslp=='y') then 
		Write(2,*) " "
		Write(2,*) "#echo  ------------------------------------------"
		Write(2,*) "echo"		
		Write(2,*) " echo '---> rsl-curves.gmt: Postscript images of RSL data vs predictions'"
		Write(2,*) "#echo  ------------------------------------------"		
!		
		file_gmt="rsl-curves.gmt"
		CALL MAKE_RSL (NV, CODE, RSL_FILE, RSL_DATABASE_FORMAT, RUN, NINC, & 
		               NRSL, TITLICE, RESOLUTION, ITER, MODE, DEGREE, VSTRING, FILE_GMT, SHORT_VISCO_FILENAME)
		if(option_gmt=='y') Write(2,*) 'sh ', trim(adjustl(file_gmt))
		Write(2,*) "mv ", trim(adjustl(file_gmt)), " "//depot//"/rsl/rsl-curves/"
		Write(2,*) "mv tmpb* ", depot//"/rsl/rsl-curves/"
	        if(option_gmt=='y') Write(2,*) "mv rslp*.ps ", depot//"/rsl/rsl-curves/ps"
	        if(option_gmt=='y') Write(2,*) "mv rslp*.pdf ", depot//"/rsl/rsl-curves/pdf"
		Write(2,*) "/bin/rm junky*"
		Endif
!
	        Write(2,*) "mv rsld*.dat ", depot//"/rsl/rsl-curves/"
	        Write(2,*) "mv rslp*.dat ", depot//"/rsl/rsl-curves/"
!
	Write(2,*) "cp ", trim(adjustl(rsl_file)), " "//depot//"/rsl/rsl-curves/"
!
	If(option_rslsca=='n') then 
		Write(2,*) "mv scatter-data.dat ",  " "//depot//"/rsl/rsl-scplot/"
		Write(2,*) "mv scatter-pred.dat ",  " "//depot//"/rsl/rsl-scplot/"
	endif
!
	ENDIF 
!
! --- A scatterplot of all RSL data vs all observations...  
	If(option_rslsca=='y') then 
		Write(2,*) " "
		Write(2,*) "#echo  ---------------------------------"
		Write(2,*) "echo"
		Write(2,*) " echo '---> Drawing a RSL scatterplot'"
		Write(2,*) "#echo  ---------------------------------"		
!		
		file_gmt="rsl-sca.gmt"
		CALL MAKE_RSLSCA (NV, CODE, RSL_FILE, RUN, NINC, NRSL, TITLICE, RESOLUTION, & 
			          ITER, MODE, DEGREE, VSTRING, FILE_GMT, SHORT_VISCO_FILENAME)
		if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file_gmt)) 
		if(option_gmt=='y') Write(2,*) "ps2pdf scatter-plot.ps"
		if(option_gmt=='y') Write(2,*) "mv scatter-plot.ps ",  depot//"/rsl/rsl-scplot/"
		if(option_gmt=='y') Write(2,*) "mv scatter-plot.pdf ", depot//"/rsl/rsl-scplot/"
		Write(2,*) "mv scatter-data.dat ", depot//"/rsl/rsl-scplot/"
		Write(2,*) "mv scatter-pred.dat ", depot//"/rsl/rsl-scplot/"
		Write(2,*) "mv title.tmp ", 	   depot//"/rsl/rsl-scplot/"
		Write(2,*) "mv ", trim(adjustl(file_gmt)), " ", depot//"/rsl/rsl-scplot/"
	ENDIF 
!
! --- Misfit study for RSL 
	If(option_rslmf=='y') then
		Write(2,*) " "
		Write(2,*) "#echo  -------------------------------"
		Write(2,*) "echo"		
		Write(2,*) " echo '---> Misfit analysis for RSL'"	
		Write(2,*) "#echo  -------------------------------"		
		file_gmt='rsl-mis.gmt'
        	Call MAKE_RSLMIS (NV, CODE, RUN, NINC, NRSL, TITLICE, RESOLUTION, ITER, & 
			          MODE, DEGREE, VSTRING, FILE_GMT, SHORT_VISCO_FILENAME)
		if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file_gmt)) 
		if(option_gmt=='y') Write(2,*) "ps2pdf rsl-misfit.ps"
		if(option_gmt=='y') Write(2,*) "mv rsl-misfit.ps ",  depot//"/rsl/rsl-misfit/"
		if(option_gmt=='y') Write(2,*) "mv rsl-misfit.pdf ", depot//"/rsl/rsl-misfit/"
		Write(2,*) "mv rsl-mis.gmt ",    depot//"/rsl/rsl-misfit/"
		Write(2,*) "mv tmph.dat ",       depot//"/rsl/rsl-misfit/"
		Write(2,*) "mv mis.dat ",        depot//"/rsl/rsl-misfit/"
		Write(2,*) "mv gmis.dat ",       depot//"/rsl/rsl-misfit/"
	ENDIF 
!
! --- RSL table with all observations & predictions   
	If(option_rsltab=='y') then 
		Write(2,*) " "
		Write(2,*) "#echo  ---------------------------------------------"
		Write(2,*) "echo"		
		Write(2,*) " echo '---> RSL data and predictions table'"	
		Write(2,*) "#echo  ---------------------------------------------"
		Write(2,*) "mv rsl-obs-vs-predictions.dat ",  depot//"/rsl/rsl-table-obs-vs-predictions/"
	ENDIF
!
! --- <<Global RSL zones>> 
	If(option_rslz=='y') then
		Write(2,*) " "
		Write(2,*) "#echo  -------------------------------------"
		Write(2,*) "echo"					            
		Write(2,*) " echo '---> rsl-zones.gmt: Global RSL zones'"
		Write(2,*) "#echo  -------------------------------------"				    
		file1_gmt='rsl-zones.gmt'
		file2_gmt='rsl-allzones.gmt'
        	Call MAKE_RSLZONES (NV, CODE, RUN, NINC, NRSL, TITLICE, RESOLUTION, ITER, & 
			            MODE, DEGREE, VSTRING, FILE1_GMT, FILE2_GMT, SHORT_VISCO_FILENAME)
		Write(2,*) "./rslzones.exe"
		if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file1_gmt)) 
		if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file2_gmt)) 
		Write(2,*) "mv ", trim(adjustl(file1_gmt)), " "//depot//"/rsl/rsl-zones/" 
		Write(2,*) "mv ", trim(adjustl(file2_gmt)), " "//depot//"/rsl/rsl-zones/" 
		Write(2,*) "mv tmpz*", " "//depot//"/rsl/rsl-zones/" 				
		Write(2,*) "mv lonlat-*.dat ", " "//depot//"/rsl/rsl-zones/" 
		Write(2,*) "mv rslzones-*.dat ", " "//depot//"/rsl/rsl-zones/" 
		if(option_gmt=='y') Write(2,*) "mv rslzones-*.pdf ", " "//depot//"/rsl/rsl-zones/" 		
		if(option_gmt=='y') Write(2,*) "mv rslzones-*.ps ", " "//depot//"/rsl/rsl-zones/" 	
	ENDIF 
!
! --- Regional RSL contour lines... 
	If(option_rslc=='y') then
		Write(2,*) " "
		Write(2,*) "#echo  ----------------------------------"
		Write(2,*) "echo"					            			    
		Write(2,*) " echo '---> RSLC.F90: Regional RSL contour lines'"
		Write(2,*) "#echo  ----------------------------------"
	        Write(2,*) "./shrslc.exe"
		Write(2,*) "./rslc.exe" 
		file_gmt='rslc.gmt'
		Call MAKE_RSLC (TIME_BPCC, MIN_RSLC, MAX_RSLC, RSL_INT, & 
				LONMINC, LONMAXC, LATMINC, LATMAXC, & 
			        NV, CODE, RUN, NINC, NRSLC, TITLICE, RESOLUTION, & 
				ITER, MODE, DEGREE, VSTRING, NAME_OF_REGION, &
				FILE_GMT, SHORT_VISCO_FILENAME) 
		if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file_gmt)) 	
 		Write(2,*) "cp ", trim(adjustl(file_region)), " "//depot//"/rsl/rsl-contours/"
		Write(2,*) "mv lonlat_rslc.dat", " "//depot//"/rsl/rsl-contours/"	
		Write(2,*) "mv rslc-cont.dat", " "//depot//"/rsl/rsl-contours/"
		Write(2,*) "mv rslc.dat", " "//depot//"/rsl/rsl-contours/"
		Write(2,*) "/bin/rm shrslc.bin"
		Write(2,*) "mv ", trim(adjustl(file_gmt)), " "//depot//"/rsl/rsl-contours/"	
		if(option_gmt=='y') Write(2,*) "mv rslc-map.ps", " "//depot//"/rsl/rsl-contours/"
		Write(2,*) "mv rslc*.tmp", " "//depot//"/rsl/rsl-contours/"							  
		if(option_gmt=='y') Write(2,*) "mv pal_rslc.cpt", " "//depot//"/rsl/rsl-contours/"	
	ENDIF
!
!
! ==============================================
! EXE 13 --- Sea level change at tide-gauges ---
! ==============================================
!
! Tide gauge database 
If(option_tga=='y') then 
	Write(2,*) " "
	Write(2,*) "#echo -------------------------------------------------"
	Write(2,*) "echo"
	Write(2,*) " echo '---> The tide-gauge database is: '", trim(adjustl(TGAUGES_FILE)) 
	Write(2,*) "#echo -------------------------------------------------"
ENDIF
!
!
! Plotting map showing the distribution of tide-gauges...
	If(option_tgplot=='y') then 
		Write(2,*) " "
		Write(2,*) "#echo ------------------------------------------------"  ! Revised May 2011
		Write(2,*) "echo"		
		Write(2,*) " echo '---> tgauges.gmt: Plotting the distribution of tide gauges'"
		Write(2,*) "#echo ------------------------------------------------"  ! Revised May 2011	
		file_gmt="tgauges.gmt"
		Call MAKE_TGAUGES (TGAUGES_FILE, TGAUGES_DATABASE_FORMAT, FILE_GMT)
		if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file_gmt)) 	
		Write(2,*) "mv ", trim(adjustl(file_gmt)), " "//depot//"/tgauges/tgauges-sites/"	
		Write(2,*) "cp ", trim(adjustl(TGAUGES_FILE)), " "//depot//"/tgauges/tgauges-sites/"
		Write(2,*) "mv lon-lat-tgauges-*.dat ",  depot//"/tgauges/tgauges-sites/"
		Write(2,*) "mv titlege*.tmp ", depot//"/tgauges/tgauges-sites/"
		Write(2,*) "mv titleall.tmp ", depot//"/tgauges/tgauges-sites/"
		Write(2,*) "mv titletgauges.tmp ", depot//"/tgauges/tgauges-sites/"
		if(option_gmt=='y') Write(2,*) "ps2pdf map-tgauges.ps "
		if(option_gmt=='y') Write(2,*) "mv map-tgauges.ps ", depot//"/tgauges/tgauges-sites/"
		if(option_gmt=='y') Write(2,*) "mv map-tgauges.pdf ", depot//"/tgauges/tgauges-sites/"
	Endif
!
!
! --- Plotting a scatterplot of observations vs years of observations 
	If(option_tgsca=='y') then 
		Write(2,*) " "
		Write(2,*) "#echo -----------------------------------------"
		Write(2,*) "echo"
		Write(2,*) " echo '---> tgauges-scpl.gmt: Scatterplot of tide gauges trends'"
		Write(2,*) "#echo -----------------------------------------"	
		file1_gmt="tgauges-scpl.gmt"
		Call MAKE_TGAUGESSCA (TGAUGES_FILE, FILE1_GMT)
!		
! --- Scatterplots and stats 
		if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file1_gmt)) 
		if(option_gmt=='y') Write(2,*) "ps2pdf tgauges-scpl.ps"
		Write(2,*) "mv ", trim(adjustl(file1_gmt)), " "//depot//"/tgauges/tgauges-scplots/"	
		Write(2,*) "cp ", trim(adjustl(TGAUGES_FILE)), " "//depot//"/tgauges/tgauges-scplots/"
		Write(2,*) "mv tgauges-scpl.dat", " "//depot//"/tgauges/tgauges-scplots/"	
		if(option_gmt=='y') Write(2,*) "mv tgauges-scpl.ps", " "//depot//"/tgauges/tgauges-scplots/"
		if(option_gmt=='y') Write(2,*) "mv tgauges-scpl.pdf", " "//depot//"/tgauges/tgauges-scplots/"	
		if(option_gmt=='y') Write(2,*) "mv tgauges-stat.dat", " "//depot//"/tgauges/tgauges-scplots/"	
		if(option_gmt=='y') Write(2,*) "mv tgauges-scpl.tmp", " "//depot//"/tgauges/tgauges-scplots/"	
		Write(2,*) "mv tmpctitle1", " "//depot//"/tgauges/tgauges-scplots/"
		Write(2,*) "mv tmpctitle2", " "//depot//"/tgauges/tgauges-scplots/"
		Write(2,*) "mv tmpctitle3", " "//depot//"/tgauges/tgauges-scplots/"
		Write(2,*) "mv tmpctitle4", " "//depot//"/tgauges/tgauges-scplots/"
		Write(2,*) "mv tmpctitle5", " "//depot//"/tgauges/tgauges-scplots/"
		Write(2,*) "/bin/rm junk-tgauges.tmp"
!		
	Endif	
!
! --- S-dot predictions at TG stations... 
	If(option_tg=='y') then 
		Write(2,*) " "
		Write(2,*) "#echo --------------------------------------------------" 
		Write(2,*) "echo"		
		Write(2,*) " echo '---> TGAUGES.F90: Dot-S, U, and N predictions at tide gauges'"
		Write(2,*) "#echo --------------------------------------------------" 
		Write(2,*) "./shtgauges.exe"
		Write(2,*) "./tgauges.exe"					
!Write(2,*) "mv shtidegauges.bin ", depot//"/bin/"
!Write(2,*) "cp shs.bin ", depot//"/bin/"	
		Write(2,*) "mv ptidegauges.dat ", " "//depot//"/tgauges/tgauges-predictions/"			
	Endif	
!
!
! ===============================================================
! ! EXE 14 ---  Elastic rebound  ---
! ===============================================================
!
!
 IF(OPTION_REB=='y') THEN 
!
    	    Write(2,*) " " 
    	    Write(2,*) "#echo -----------------------------------------------------"
    	    Write(2,*) "echo"
    	    Write(2,*) " echo '---> ELA_REB.F90: Elastic rebound '"
    	    Write(2,*) "#echo -----------------------------------------------------"		    
!
 if(option_reb_gg=='y'.or. & 
    option_reb_gr=='y'.or. & 
    option_reb_ag=='y'.or. & 
    option_reb_ar=='y'.or. & 
    option_reb_sm=='y') Write(2,*) "./elareb.exe"  
!
! \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/   !  Greenland 
! - Regional and/or global elastic rebound for >>>> Greenland  <<<< 	
! \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/   !  Greenland 
!
 if(option_reb_gg=='y'.or.option_reb_gr=='y') then 
!	
        FILE_GMT='Greenland-elareb.gmt' 
!	
	call MAKE_ELAREB_MAPS_GREENLAND (option_reb_gg, & 
					 option_reb_gr,  & 
 	 				 resolution,      & 
					 nv,               &
					 code,              & 
					 iter,              & 
					 mode, 	           & 
					 degree,          & 
					 vstring,        & 
					 TITLICE,       & 
					 FILE_GMT)
!
	Write(2,*) "if [ -f  ./gmtdefaults4  ]" 
	Write(2,*) "then"
	Write(2,*) "/bin/rm -v ./gmtdefaults4"
	Write(2,*) "fi"        
	if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file_gmt)) 
!	Write(2,*) "cp ./gmtdefaults4 ", depot//"/elastic-rebound/Greenland/" 
        Write(2,*) "mv ", trim(adjustl(file_gmt)), " ", depot//"/elastic-rebound/Greenland/" 
!
	If(option_reb_gg=='y') then 
	     if(option_gmt=='y') Write(2,*) "ps2pdf smap_gree_glob.ps "	     
	     if(option_gmt=='y') Write(2,*) "ps2pdf umap_gree_glob.ps "	     
	     if(option_gmt=='y') Write(2,*) "ps2pdf nmap_gree_glob.ps "	
	     if(option_gmt=='y') Write(2,*) "ps2pdf smap_gree_glob_norm.ps "		          
	     Write(2,*) "mv *map_gree_glob.dat ", " ", depot//"/elastic-rebound/Greenland/"  
	     if(option_gmt=='y') Write(2,*) "mv *map_gree_glob.ps  ", " ", depot//"/elastic-rebound/Greenland/" 
	     if(option_gmt=='y') Write(2,*) "mv *map_gree_glob.pdf ", " ", depot//"/elastic-rebound/Greenland/"  
	     Write(2,*) "mv *gtmp*.dat ", " ", depot//"/elastic-rebound/Greenland/"  
	     Write(2,*) "mv smap_gree_glob_norm.dat ", " ", depot//"/elastic-rebound/Greenland/" 
	     if(option_gmt=='y') Write(2,*) "mv smap_gree_glob_norm.ps  ", " ", depot//"/elastic-rebound/Greenland/" 
	     if(option_gmt=='y') Write(2,*) "mv smap_gree_glob_norm.pdf ", " ", depot//"/elastic-rebound/Greenland/"  
	endif
	If(option_reb_gr=='y') then 
	     if(option_gmt=='y') Write(2,*) "ps2pdf smap_gree_reg.ps "	     
	     if(option_gmt=='y') Write(2,*) "ps2pdf umap_gree_reg.ps "	     
	     if(option_gmt=='y') Write(2,*) "ps2pdf nmap_gree_reg.ps "	     
	     Write(2,*) "mv *_gree_reg.dat ", " ", depot//"/elastic-rebound/Greenland/" 
	     if(option_gmt=='y') Write(2,*) "mv *_gree_reg.ps  ", " ", depot//"/elastic-rebound/Greenland/" 
	     if(option_gmt=='y') Write(2,*) "mv *_gree_reg.pdf ", " ", depot//"/elastic-rebound/Greenland/"  
	     Write(2,*) "mv *rrtmp*.dat ", " ", depot//"/elastic-rebound/Greenland/"  
	endif
!
        endif   ! <<<<<<< <<<<<<< <<<<<<< <<<<<<< <<<<<<< <<<<<<< <<<<<<< On Greenland 
!
!
! \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/   !  Antarctica  
 	if(option_reb_ag=='y'.or.option_reb_ar=='y') then 
! \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/   !  Antarctica  
!XXXXXXXXX
!XXXXXXXXX
!XXXXXXXXX
!XXXXXXXXX
!XXXXXXXXX
!XXXXXXXXX
!XXXXXXXXX
! === Revised GG May 03 2011 - for the ANT component  
        FILE_GMT='Anta-elareb.gmt' 
!	
       call MAKE_ELAREB_MAPS_ANTARCTICA  (option_reb_ag, & 
					  option_reb_ar,  & 
 	 				  resolution,      & 
					  nv,               &
					  code,              & 
					  iter,              & 
					  mode,             & 
					  degree,          & 
					  vstring,        & 
					  TITLICE,       & 
					  FILE_GMT)
!
	Write(2,*) "if [ -f  ./gmtdefaults4  ]" 
	Write(2,*) "then"
	Write(2,*) "/bin/rm -v ./gmtdefaults4"
	Write(2,*) "fi"        
	if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file_gmt)) 
!	Write(2,*) "cp ./gmtdefaults4 ", depot//"/elastic-rebound/Antarctica/" 
        Write(2,*) "mv ", trim(adjustl(file_gmt)), " ", depot//"/elastic-rebound/Antarctica/" 
!
	If(option_reb_ag=='y') then 
	     if(option_gmt=='y') Write(2,*) "ps2pdf smap_anta_glob.ps "	     
	     if(option_gmt=='y') Write(2,*) "ps2pdf umap_anta_glob.ps "	     
	     if(option_gmt=='y') Write(2,*) "ps2pdf nmap_anta_glob.ps "	
	     if(option_gmt=='y') Write(2,*) "ps2pdf smap_anta_glob_norm.ps "		          
	     Write(2,*) "mv *map_anta_glob.dat ", " ", depot//"/elastic-rebound/Antarctica/"  
	     if(option_gmt=='y') Write(2,*) "mv *map_anta_glob.ps  ", " ", depot//"/elastic-rebound/Antarctica/" 
	     if(option_gmt=='y') Write(2,*) "mv *map_anta_glob.pdf ", " ", depot//"/elastic-rebound/Antarctica/"  
	     Write(2,*) "mv *gtmp*.dat ", " ", depot//"/elastic-rebound/Antarctica/"  
	     Write(2,*) "mv smap_anta_glob_norm.dat ", " ", depot//"/elastic-rebound/Antarctica/" 
	     if(option_gmt=='y') Write(2,*) "mv smap_anta_glob_norm.ps  ", " ", depot//"/elastic-rebound/Antarctica/" 
	     if(option_gmt=='y') Write(2,*) "mv smap_anta_glob_norm.pdf ", " ", depot//"/elastic-rebound/Antarctica/"  
	endif
	If(option_reb_ar=='y') then 
	     if(option_gmt=='y') Write(2,*) "ps2pdf smap_anta_reg.ps "	     
	     if(option_gmt=='y') Write(2,*) "ps2pdf umap_anta_reg.ps "	     
	     if(option_gmt=='y') Write(2,*) "ps2pdf nmap_anta_reg.ps "	     
	     Write(2,*) "mv *_anta_reg.dat ", " ", depot//"/elastic-rebound/Antarctica/" 
	     if(option_gmt=='y') Write(2,*) "mv *_anta_reg.ps  ", " ", depot//"/elastic-rebound/Antarctica/" 
	     if(option_gmt=='y') Write(2,*) "mv *_anta_reg.pdf ", " ", depot//"/elastic-rebound/Antarctica/"  
	     Write(2,*) "mv *rrtmp*.dat ", " ", depot//"/elastic-rebound/Antarctica/"  
	endif
!
        endif   ! <<<<<<< <<<<<<< <<<<<<< <<<<<<< <<<<<<< <<<<<<< <<<<<<< On Antarctica 
!XXXXXXXXX
!XXXXXXXXX
!XXXXXXXXX
!XXXXXXXXX
!XXXXXXXXX
!XXXXXXXXX
!XXXXXXXXX
!
!
! \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/   !  Small glaciers   
 	if(option_reb_sm=='y') then 
! \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/   !  Small glaciers  
!
        FILE_GMT='Glac-elareb.gmt' 
!
	call MAKE_ELAREB_MAPS_GLAC  (option_reb_sm, & 
                                     resolution,    & 
                                     nv,            &
                                     code,          & 
				     iter,          & 
				     mode,          & 
				     degree,        & 
				     vstring,       & 
			             TITLICE,       & 
				     FILE_GMT)
!
        if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file_gmt)) 
        Write(2,*) "mv ", trim(adjustl(file_gmt)), " ", depot//"/elastic-rebound/Small_glaciers/" 
!
	If(option_reb_sm=='y') then 
	     if(option_gmt=='y') Write(2,*) "ps2pdf smap_glac_glob.ps "	     
	     if(option_gmt=='y') Write(2,*) "ps2pdf umap_glac_glob.ps "	     
	     if(option_gmt=='y') Write(2,*) "ps2pdf nmap_glac_glob.ps "	
	     if(option_gmt=='y') Write(2,*) "ps2pdf smap_glac_glob_norm.ps "		          
	     Write(2,*) "mv *map_glac_glob.dat ", " ", depot//"/elastic-rebound/Small_glaciers/"  
	     if(option_gmt=='y') Write(2,*) "mv *map_glac_glob.ps  ", " ", depot//"/elastic-rebound/Small_glaciers/" 
	     if(option_gmt=='y') Write(2,*) "mv *map_glac_glob.pdf ", " ", depot//"/elastic-rebound/Small_glaciers/"  
	     Write(2,*) "mv glc*tmp*.dat ", " ", depot//"/elastic-rebound/Small_glaciers/"  
	     Write(2,*) "mv smap_glac_glob_norm.dat ", " ", depot//"/elastic-rebound/Small_glaciers/" 
	     if(option_gmt=='y') Write(2,*) "mv smap_glac_glob_norm.ps  ", " ", depot//"/elastic-rebound/Small_glaciers/" 
	     if(option_gmt=='y') Write(2,*) "mv smap_glac_glob_norm.pdf ", " ", depot//"/elastic-rebound/Small_glaciers/"  
	endif
!
!
        endif
!		 
! 
! \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/   !  Geodetic Sites   
 	if(option_3d_reb=='y') then 
! \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/   !  Geodetic Sites   
	Write(2,*) "./geoelareb.exe" 
        Write(2,*) "mv geodetic-predictions.dat ", depot//"/elastic-rebound/Geodetic/" 
        Write(2,*) "cp ", file_3d, depot//"/elastic-rebound/Geodetic/"	
			       Endif
!
	Endif   ! <<< <<< <<< <<< On the Elastic Rebound switch
!
!
!
! =======================================================================
! ! EXE YY ---  Global maps of S, and dot-S for elastic rebound in steps
! =======================================================================
!
	if(option_ela_ms=='y') then 
		Write(2,*) " " 
		Write(2,*) "#echo -----------------------------------------------------"
		Write(2,*) "echo"
		Write(2,*) " echo '---> Global maps for elastic rebound in steps'"
		Write(2,*) "#echo -----------------------------------------------------"		
!
	        if(option_gm_ela_ms_s_var=='y'.or.&
		   option_gm_ela_ms_n_var=='y'.or.&
		   option_gm_ela_ms_u_var=='y') write(2,*) "./gmaps_ela_ms_var.exe"
	        if(option_gm_ela_ms_s_dot=='y'.or.&
		   option_gm_ela_ms_n_dot=='y'.or.&
		   option_gm_ela_ms_u_dot=='y') write(2,*) "./gmaps_ela_ms_dot.exe"
!
		file1_gmt='gmaps_ela_ms_var.gmt'  
		file2_gmt='gmaps_ela_ms_dot.gmt'  
!
        	Call MAKE_GMAPS_ELA_MS    (DEPOT,      & 
		                           TITLICE,    & 
		                           RESOLUTION, & 
					   DEGREE,     & 
					   NINC,       & 
					   CODE,       & 
					   MODE,       &
					   ITER,       &
					   FILE1_GMT,   & 
                                           FILE2_GMT,   & 
					   OPTION_GIA_CORR, & 		
		                           OPTION_GM_ELA_MS_S_VAR, & 
					   OPTION_GM_ELA_MS_N_VAR, &
					   OPTION_GM_ELA_MS_U_VAR, & 					   
					   OPTION_GM_ELA_MS_S_DOT, &  
					   OPTION_GM_ELA_MS_N_DOT, &
					   OPTION_GM_ELA_MS_U_DOT) 
!					   
! WORK IN PROGRESS HERE!!!!!!
!
                if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file1_gmt)) 
		Write(2,*) "mv ", trim(adjustl(file1_gmt)), " "//depot//"/gmaps-ela-ms/var/" 
!
                if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file2_gmt)) 
		Write(2,*) "mv ", trim(adjustl(file2_gmt)), " "//depot//"/gmaps-ela-ms/dot/" 
!
		if(option_gm_ela_ms_s_var=='y')Write(2,*) "mv svarmap*.dat ",  depot//"/gmaps-ela-ms/var/"
		if(option_gm_ela_ms_n_var=='y')Write(2,*) "mv nvarmap*.dat ",  depot//"/gmaps-ela-ms/var/"
		if(option_gm_ela_ms_u_var=='y')Write(2,*) "mv uvarmap*.dat ",  depot//"/gmaps-ela-ms/var/"
        if(option_gmt=='y') then
  		  if(option_gm_ela_ms_s_var=='y')Write(2,*) "mv svarmap*.ps  ",  depot//"/gmaps-ela-ms/var/"
		  if(option_gm_ela_ms_n_var=='y')Write(2,*) "mv nvarmap*.ps  ",  depot//"/gmaps-ela-ms/var/"
		  if(option_gm_ela_ms_u_var=='y')Write(2,*) "mv uvarmap*.ps  ",  depot//"/gmaps-ela-ms/var/"
        endif
!
		if(option_gm_ela_ms_s_dot=='y')Write(2,*) "mv sdotmap*.dat ",  depot//"/gmaps-ela-ms/dot/"
		if(option_gm_ela_ms_n_dot=='y')Write(2,*) "mv ndotmap*.dat ",  depot//"/gmaps-ela-ms/dot/"
		if(option_gm_ela_ms_u_dot=='y')Write(2,*) "mv udotmap*.dat ",  depot//"/gmaps-ela-ms/dot/"
        if(option_gmt=='y') then
		  if(option_gm_ela_ms_s_dot=='y')Write(2,*) "mv sdotmap*.ps  ",  depot//"/gmaps-ela-ms/dot/"
		  if(option_gm_ela_ms_n_dot=='y')Write(2,*) "mv ndotmap*.ps  ",  depot//"/gmaps-ela-ms/dot/"
		  if(option_gm_ela_ms_u_dot=='y')Write(2,*) "mv udotmap*.ps  ",  depot//"/gmaps-ela-ms/dot/"
        endif
!
                if(option_gm_ela_ms_s_dot=='y')Write(2,*) "mv tmpfs-sd* ", depot//"/gmaps-ela-ms/dot/"
                if(option_gm_ela_ms_n_dot=='y')Write(2,*) "mv tmpfs-nd* ", depot//"/gmaps-ela-ms/dot/"
                if(option_gm_ela_ms_u_dot=='y')Write(2,*) "mv tmpfs-ud* ", depot//"/gmaps-ela-ms/dot/"

                if(option_gm_ela_ms_s_var=='y')Write(2,*) "mv tmpfs-s* ",  depot//"/gmaps-ela-ms/var/"
                if(option_gm_ela_ms_n_var=='y')Write(2,*) "mv tmpfs-n* ",  depot//"/gmaps-ela-ms/var/"
                if(option_gm_ela_ms_u_var=='y')Write(2,*) "mv tmpfs-u* ",  depot//"/gmaps-ela-ms/var/"
!	     
        endif
	
!
!
! =======================================================================
! ! EXE YY ---  Tide gauges predictions for elastic rebound in steps
! =======================================================================
!
	if(OPTION_TG_ELA_MS=='y') then 
		Write(2,*) " " 
		Write(2,*) "#echo -------------------------------------------------------------"
		Write(2,*) "echo"
		Write(2,*) " echo '---> TG predictions (vars and rates) for elastic rebound in steps'"
		Write(2,*) "#echo -------------------------------------------------------------"		
!
                Write(2,*) "./shtgauges.exe" 
                Write(2,*) "./tgaugeselams.exe" 
!
		file_gmt='tg_plots_var.gmt'
!
        	Call MAKE_TG_PLOTS_VAR    (FILE_TG_ELA_MS,             & 
		                           ELAREB_MS_DATABASE_FORMAT,  & 
					   FILE_GMT) 
!        	Call MAKE_TG_PLOTS_DOT    (FILE_TG_ELA_MS,             & 
!		                           ELAREB_MS_DATABASE_FORMAT,  & 
!					   FILE_GMT) 


		if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file_gmt)) 
		Write(2,*) "mv ", trim(adjustl(file_gmt)), " "//depot//"/tgauges-ela-ms" 
		Write(2,*) "mv ptidegauges.dat", " "//depot//"/tgauges-ela-ms" 
		Write(2,*) "mv *.tga",    " "//depot//"/tgauges-ela-ms" 
		if(option_gmt=='y') Write(2,*) "mv *.tga.ps", " "//depot//"/tgauges-ela-ms" 
!
        endif
!
!
! ===============================================================
! ! EXE 15 ---  S, U, and N-dot at present time (Global maps) ---
! ===============================================================
!
	if(option_gm=='y') then 
		Write(2,*) " " 
		Write(2,*) "#echo -----------------------------------------------------"
		Write(2,*) "echo"
		Write(2,*) " echo '---> GMAPS.F90: Global maps of dot S, U and N at present time'"
		Write(2,*) "#echo -----------------------------------------------------"		
		file_gmt='gmaps.gmt'
		Write(2,*) "./gmaps.exe"  

		Write(2,*) "cp GMT_panoply.cpt ",  depot//"/gmaps/"

        	Call MAKE_GMAPS (TITLICE, RESOLUTION, NV, CODE, ITER, MODE,                             & 
				 DEGREE, VSTRING, OPTION_ROF, FILE_CAP, FILE_GMT, SHORT_VISCO_FILENAME, & 
                                 OPTION_DOTS,   & 
				 OPTION_DOTU,   & 
				 OPTION_DOTN,   & 
				 OPTION_DOTG,   & 
                                 OPTION_DOTW,   & 
				 OPTION_DOTFA,  & 
				 OPTION_DOTSS,  &
				 OPTION_DOTLOI, & 
				 OPTION_DOTLOO, & 
				 OPTION_DOTLOT)				 
!				 				 
		if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file_gmt)) 
		Write(2,*) "mv ", trim(adjustl(file_gmt)), " "//depot//"/gmaps" 
!
		if(option_dots== 'y')Write(2,*) "mv sdotmap.dat ",  depot//"/gmaps/"
		if(option_dotu== 'y')Write(2,*) "mv udotmap.dat ",  depot//"/gmaps/"
		if(option_dotn== 'y')Write(2,*) "mv ndotmap.dat ",  depot//"/gmaps/"
		if(option_dotg== 'y')Write(2,*) "mv gdotmap.dat ",  depot//"/gmaps/"
		if(option_dotw== 'y')Write(2,*) "mv wdotmap.dat ",  depot//"/gmaps/"
		if(option_dotfa=='y')Write(2,*) "mv fadotmap.dat ", depot//"/gmaps/"
		if(option_dotss=='y')Write(2,*) "mv ssdotmap.dat ", depot//"/gmaps/"
!
        if(option_gmt=='y') then
  		  if(option_dots== 'y')Write(2,*) "ps2pdf  sdotmap.ps"
		  if(option_dotu== 'y')Write(2,*) "ps2pdf  udotmap.ps"
		  if(option_dotn== 'y')Write(2,*) "ps2pdf  ndotmap.ps"		
		  if(option_dotg== 'y')Write(2,*) "ps2pdf  gdotmap.ps"		
		  if(option_dotw== 'y')Write(2,*) "ps2pdf  wdotmap.ps"		
		  if(option_dotfa=='y')Write(2,*) "ps2pdf fadotmap.ps"		
		  if(option_dotss=='y')Write(2,*) "ps2pdf ssdotmap.ps"		
        endif
!
        if(option_gmt=='y') then
		  if(option_dots== 'y')Write(2,*) "mv  sdotmap.ps ", depot//"/gmaps/"
		  if(option_dotu== 'y')Write(2,*) "mv  udotmap.ps ", depot//"/gmaps/"
		  if(option_dotn== 'y')Write(2,*) "mv  ndotmap.ps ", depot//"/gmaps/"
		  if(option_dotg== 'y')Write(2,*) "mv  gdotmap.ps ", depot//"/gmaps/"
		  if(option_dotw== 'y')Write(2,*) "mv  wdotmap.ps ", depot//"/gmaps/"
		  if(option_dotfa=='y')Write(2,*) "mv fadotmap.ps ", depot//"/gmaps/"
		  if(option_dotss=='y')Write(2,*) "mv ssdotmap.ps ", depot//"/gmaps/"
        endif
!
        if(option_gmt=='y') then
		  if(option_dots== 'y')Write(2,*) "mv  sdotmap.pdf ", depot//"/gmaps/"
		  if(option_dotu== 'y')Write(2,*) "mv  udotmap.pdf ", depot//"/gmaps/"
		  if(option_dotn== 'y')Write(2,*) "mv  ndotmap.pdf ", depot//"/gmaps/"
		  if(option_dotg== 'y')Write(2,*) "mv  gdotmap.pdf ", depot//"/gmaps/"
		  if(option_dotw== 'y')Write(2,*) "mv  wdotmap.pdf ", depot//"/gmaps/"
		  if(option_dotfa=='y')Write(2,*) "mv fadotmap.pdf ", depot//"/gmaps/"
		  if(option_dotss=='y')Write(2,*) "mv ssdotmap.pdf ", depot//"/gmaps/"
        endif
!
		Write(2,*) "mv tmpf* ", depot//"/gmaps/"
                Write(2,*) "mv pale.cpt ", depot//"/gmaps/"	
		 
                if(option_rof=='z')Write(2,*) "mv shore.dat ", depot//"/gmaps/"
				
	Endif
!
!
!
! ==============================================================
! EXE 16 --- S, U, and N-dot at present time (Regional maps) ---
! ==============================================================
!
	if(option_rm(0)=='y') then 
		Write(2,*) " " 
		Write(2,*) "#echo -------------------------------------------------------"
		Write(2,*) "echo"
		Write(2,*) " echo '---> RMAPS.F90: Regional maps of dot S, U and N at present time'"
		Write(2,*) "#echo -------------------------------------------------------"
		Write(2,*) "./rmaps.exe"  
		Write(2,*) " cp sdotmap.dat ", depot//"/rmaps/sun_data"
		Write(2,*) " cp udotmap.dat ", depot//"/rmaps/sun_data"
		Write(2,*) " cp ndotmap.dat ", depot//"/rmaps/sun_data"
        	Call MAKE_RMAPS (TITLICE, RESOLUTION, NV, CODE, ITER, MODE, DEGREE, VSTRING, & 
		                 OPTION_RM, SHORT_VISCO_FILENAME)
			if(option_rm(1)=='y') then 
			Write(2,*) " echo '     - Italy'"
				if(option_gmt=='y') Write(2,*) "sh italy.gmt"  
				if(option_gmt=='y') Write(2,*) "ps2pdf sdot-italy.ps"	
				if(option_gmt=='y') Write(2,*) "ps2pdf udot-italy.ps"	
				if(option_gmt=='y') Write(2,*) "ps2pdf ndot-italy.ps"
				Write(2,*) " mv italy.gmt ", depot//"/rmaps/Italy"
				Write(2,*) "mv *dot*ital* ", depot//"/rmaps/Italy"
				Write(2,*) "mv *pale*ital* ",depot//"/rmaps/Italy"		
			Endif			
			if(option_rm(2)=='y') then 
			Write(2,*) " echo '     - Mediterranean'"
				if(option_gmt=='y') Write(2,*) "sh mediterranean.gmt"  	
				if(option_gmt=='y') Write(2,*) "ps2pdf sdot-mediterranean.ps"	
				if(option_gmt=='y') Write(2,*) "ps2pdf udot-mediterranean.ps"	
				if(option_gmt=='y') Write(2,*) "ps2pdf ndot-mediterranean.ps"
				Write(2,*) " mv mediterranean.gmt ", depot//"/rmaps/Mediterranean"
				Write(2,*) "mv *pale*medi* ",        depot//"/rmaps/Mediterranean"
				Write(2,*) "mv *dot*medi* ",         depot//"/rmaps/Mediterranean"			
			Endif			
			if(option_rm(3)=='y') then 
			Write(2,*) " echo '     - Europe'"
				if(option_gmt=='y') Write(2,*) "sh europe.gmt"  
				if(option_gmt=='y') Write(2,*) "ps2pdf sdot-europe.ps"	
				if(option_gmt=='y') Write(2,*) "ps2pdf udot-europe.ps"	
				if(option_gmt=='y') Write(2,*) "ps2pdf ndot-europe.ps"
				Write(2,*) " mv europe.gmt ", depot//"/rmaps/Europe"
				Write(2,*) "mv *pale*euro* ", depot//"/rmaps/Europe"
				Write(2,*) "mv *dot*euro* ",  depot//"/rmaps/Europe"
			Endif							
			if(option_rm(4)=='y') then 
			Write(2,*) " echo '     - Fennoscandia'"
				if(option_gmt=='y') Write(2,*) "sh fennoscandia.gmt"  	
				if(option_gmt=='y') Write(2,*) "ps2pdf sdot-fennoscandia.ps"	
				if(option_gmt=='y') Write(2,*) "ps2pdf udot-fennoscandia.ps"	
				if(option_gmt=='y') Write(2,*) "ps2pdf ndot-fennoscandia.ps"
				Write(2,*) " mv fennoscandia.gmt ", depot//"/rmaps/Fennoscandia"
				Write(2,*) "mv *pale*fenn* ",       depot//"/rmaps/Fennoscandia"
				Write(2,*) "mv *dot*fenn* ",        depot//"/rmaps/Fennoscandia"
			Endif					
			if(option_rm(5)=='y') then 
			Write(2,*) " echo '     - Greenland'"
				if(option_gmt=='y') Write(2,*) "sh greenland.gmt"  	
				if(option_gmt=='y') Write(2,*) "ps2pdf sdot-greenland.ps"	
				if(option_gmt=='y') Write(2,*) "ps2pdf udot-greenland.ps"	
				if(option_gmt=='y') Write(2,*) "ps2pdf ndot-greenland.ps"
				Write(2,*) " mv greenland.gmt ", depot//"/rmaps/Greenland"
				Write(2,*) "mv *pale*gree* ",    depot//"/rmaps/Greenland"
				Write(2,*) "mv *dot*gree* ",     depot//"/rmaps/Greenland"
			Endif					
			if(option_rm(6)=='y') then 
			Write(2,*) " echo '     - North America'"
				if(option_gmt=='y') Write(2,*) "sh north-america.gmt"  	
				if(option_gmt=='y') Write(2,*) "ps2pdf sdot-namerica.ps"	
				if(option_gmt=='y') Write(2,*) "ps2pdf udot-namerica.ps"	
				if(option_gmt=='y') Write(2,*) "ps2pdf ndot-namerica.ps"
				Write(2,*) " mv north-america.gmt ", depot//"/rmaps/North_America"
			        Write(2,*) "mv *pale*name* ",        depot//"/rmaps/North_America"
				Write(2,*) "mv *dot*name* ",         depot//"/rmaps/North_America"
			Endif
			if(option_rm(7)=='y') then 
			Write(2,*) " echo '     - Antarctica'"
				if(option_gmt=='y') Write(2,*) "sh antarctica.gmt"  	
				if(option_gmt=='y') Write(2,*) "ps2pdf sdot-antarctica.ps"	
				if(option_gmt=='y') Write(2,*) "ps2pdf udot-antarctica.ps"	
				if(option_gmt=='y') Write(2,*) "ps2pdf ndot-antarctica.ps"
				Write(2,*) " mv antarctica.gmt ", depot//"/rmaps/Antarctica"
				Write(2,*) "mv *pale*anta* ",     depot//"/rmaps/Antarctica"
				Write(2,*) "mv *dot*anta* ",      depot//"/rmaps/Antarctica"
			Endif	
	Endif
!
!
!
! ======================================================================
! EXE 17 --- 3D velocity & S and N-dot at present time at specific sites 
! ======================================================================
!
	if(option_3d=='y') then 
		Write(2,*) " " 
		Write(2,*) "#echo -----------------------------------------------------------"
		Write(2,*) "echo"
		Write(2,*) " echo '---> GEO.F90: 3D velocity and S and N-dot today at specific sites'"
		Write(2,*) "#echo -----------------------------------------------------------"
		Write(2,*) "./geo.exe"  
		Write(2,*) "cp ", file_3d, depot//"/geod/sites" 		
		Write(2,*) "mv geodetic-predictions.dat ", depot//"/geod/sites" 		
	Endif 
!
!
!
! ========================================================
! EXE 18 --- 3D velocity at present time at points on maps  
! ========================================================
!
	if(option_3d_regions=='y') then 
		Write(2,*) " " 
		Write(2,*) "#echo -----------------------------------------------------------"
		Write(2,*) "echo"
		Write(2,*) " echo '---> GEO_MAPS.F90: 3D velocity and U-dot today at points on maps'"
		Write(2,*) "#echo -----------------------------------------------------------"
		Write(2,*) "./geomaps.exe"  
		
		
		
!Write(2,*) "cp ", file_3d, depot//"/geod/sites" 		
!Write(2,*) "mv geodetic-predictions.dat ", depot//"/geod/sites" 		
	Endif 
!
!
!
!
! =========================================================================
! EXE 19 ---  Rate of variation of Stokes coefficients at present time  ---
! =========================================================================
!
	if(option_st=='y') then 
		Write(2,*) " " 
		Write(2,*) "#echo ----------------------------------------------------------------"
		Write(2,*) "echo"
		Write(2,*) " echo '---> STOKES.F90: Present-time rate of change of Stokes coefficients'"
		Write(2,*) "#echo ----------------------------------------------------------------"		
		file_gmt='stokes.gmt'
		Write(2,*) "./stokes.exe"  					   
                Call MAKE_STOKES (RUN, NV, CODE, TITLICE, RESOLUTION, ITER, & 
			          MODE, DEGREE, VSTRING, SHORT_VISCO_FILENAME, FILE_GMT)
		if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file_gmt)) 
		Write(2,*) "mv ", trim(adjustl(file_gmt)), " "//depot//"/stokes/"  
		if(option_gmt=='y') Write(2,*) "ps2pdf stokes.ps"
		if(option_gmt=='y') Write(2,*) "mv stokes.ps ", depot//"/stokes/"
		if(option_gmt=='y') Write(2,*) "mv stokes.pdf ", depot//"/stokes/"
		Write(2,*) "mv stokes.dat ", depot//"/stokes/"
		if(option_gmt=='y') Write(2,*) "mv cosine.tmp ", depot//"/stokes/"
		if(option_gmt=='y') Write(2,*) "mv sine.tmp ", depot//"/stokes/"
		Write(2,*) "mv stokes*.tmp ", depot//"/stokes/"	
		Write(2,*) "mv title_stokes.tmp ", depot//"/stokes/"		
	Endif
!
!
!
Write(2,*) " "
Write(2,*) " "
Write(2,*) "echo ''"
Write(2,*) "echo '---------------------------------'"
Write(2,*) "echo ' >>> 2. Cleaning the directory...'"
Write(2,*) "echo '---------------------------------'"
!
Write(2,*) "if [ -f  ./gmtdefaults4  ]" 
Write(2,*) "then"
Write(2,*) "/bin/rm -v ./gmtdefaults4"
Write(2,*) "fi"
!
write(2,*) "/bin/rm *.mod"
write(2,*) "/bin/rm *.mod" 
write(2,*) "/bin/rm ./ALMA/main_module.mod" 
write(2,*) "/bin/rm ./ALMA/common.mod" 
!
Write(2,*) "/bin/rm -v *brok*.bin"
!
Write(2,*) "/bin/rm -v *.exe"
Write(2,*) "/bin/rm -v *.o"
!
Write(2,*) "/bin/rm -v *.grd"
Write(2,*) "/bin/rm -v *jun*"  
!
If(option_nm=='y') Write(2,*) "mv common*.mod ", depot//"/MODULES/"
If(option_nm=='y') Write(2,*) "mv strata.mod ", depot//"/MODULES/"
IF(OPTION_pw=='y') Write(2,*) "mv alma-logfile.dat ",  depot//"/log/" 
IF(OPTION_pw=='y') Write(2,*) "/bin/rm main_module.mod "
Write(2,*) "cp shtools.mod ", depot//"/MODULES/"
!
If(option_nm=='y') Write(2,*) "mv taboo.log ",    depot//"/log/"
Write(2,*) "mv selen.log ",    depot//"/log/"
Write(2,*) "cp config.dat ",   depot//"/log/"
Write(2,*) "cp config-MOD.f90 ",   depot//"/log/"
Write(2,*) "cp data.inc ",     depot//"/log/"
Write(2,*) "cp selen.sh ",     depot//"/log/"
!
Write(2,*) "mv pxa.dat ",      depot//"/px"
Write(2,*) "mv weta.dat ",     depot//"/px/"
Write(2,*) "mv drya.dat ",     depot//"/px/"
Write(2,*) "mv px-lat.dat ",   depot//"/px/"
Write(2,*) "mv px.dat ",       depot//"/px/"
Write(2,*) "mv px-table.dat ", depot//"/px/"
!
If    (option_nm=='y')then
! 
	Write(2,*) "mv ebu.dat ", depot//"/Love-Numbers-by-TABOO"
	Write(2,*) "mv ebv.dat ", depot//"/Love-Numbers-by-TABOO"
	Write(2,*) "mv ebs.dat ", depot//"/Love-Numbers-by-TABOO"
	Write(2,*) "mv ebn.dat ", depot//"/Love-Numbers-by-TABOO"
elseif(option_pw=='y')then 
! 
	Write(2,*) "mv ebu.dat ", depot//"/Love-Numbers-by-ALMA"
	Write(2,*) "mv ebv.dat ", depot//"/Love-Numbers-by-ALMA"
	Write(2,*) "mv ebs.dat ", depot//"/Love-Numbers-by-ALMA/"
	Write(2,*) "mv ebn.dat ", depot//"/Love-Numbers-by-ALMA/"
Endif
!
Write(2,*) "mv ",  trim(adjustl(ice_file)), &
           "   "//depot//"/"//trim(adjustl(titlice))//"/original/"	   
!
If    (option_topo=='y')then
        Write(2,*) "mv ptopo-*.dat ",    " "//depot//"/paleo-topography/data"
        Write(2,*) "mv imask-*.dat ",    " "//depot//"/paleo-topography/data"
        Write(2,*) "mv oc-*.dat ",       " "//depot//"/paleo-topography/more-data"			
        Write(2,*) "mv shof-*.dat ",       " "//depot//"/paleo-topography/more-data"	
        Write(2,*) "mv px-table-*.dat ", " "//depot//"/paleo-topography/more-data"	
Endif
!
!
Write(2,*) "echo ''"
Write(2,*) "echo '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -'"
Write(2,*) "echo '     SELEN, a Sea levEL EquatioN solver, Version 3.2     '"
Write(2,*) "echo '                                                         '" 
Write(2,*) "echo '      Web page: http://fcolleoni.free.fr/SELEN.html      '"	    
Write(2,*) "echo '   Send comments, requests of help and suggestions to:   '" 
Write(2,*) "echo '                <giorgio.spada@gmail.com>                '"
Write(2,*) "echo '                            -                            '"
Write(2,*) "echo '                    Copyright(C) 2008                    '"    
Write(2,*) "echo '     Giorgio Spada, Florence Colleoni & Paolo Stocchi    '"
Write(2,*) "echo '                          * * *                          '"
Write(2,*) "echo '     This programs comes with  ABSOLUTELY NO WARRANTY    '"
Write(2,*) "echo ' This is free software and you are welcome to distribute '"
Write(2,*) "echo '              it under certain conditions.               '"
Write(2,*) "echo '    For details, visit  <http://www.gnu.org/licenses/>   '"
Write(2,*) "echo '                  or edit file COPYING                   '"
Write(2,*) "echo '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -'"
Write(2,*) "echo ''"
!
!
Write(2,*) "echo ''"
Write(2,*) "echo ' >>> Outputs for this run are available in directory: '", trim(adjustl(depot))
!
!Write(2,*) "echo ''"
!
!
! --- Closing "selen.sh"
   close(2)
!
! --- Closing "config.dat"
   close(1)
!
!
!
!
! =========================================================
!
! 	Part #3 : preparing the include file 'data.inc'
!
! =========================================================
!
! 
! >>>>>> A time stamp on 'data.inc'
!
  call DATE_AND_TIME (date,timc)      
!
open(3,file='data.inc',status='unknown')
!
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!", " File <<data.inc>>, created by <<config.f90>> on ", & 
             date(1:4), '.', date(5:6), '.', date(7:8), ' ', & 
	     timc(1:2), '.', timc(3:4), '.', timc(5:6) 
Write(3,*)"!                  				                                "
Write(3,*)"! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
Write(3,*)"! Copyright (C) 2008 Giorgio Spada, Florence Colleoni, and Paolo Stocchi     "
Write(3,*)"! 								                "								
Write(3,*)"! This file is part of SELEN.                                                " 
Write(3,*)"!       							                " 
Write(3,*)"! SELEN is free software: you can redistribute it and/or modify it under the "
Write(3,*)"! terms of the GNU General Public License as published by the Free Software  " 
Write(3,*)"! Foundation, either version 3 of the License, or at your option) any later  " 
Write(3,*)"! version. 									"
Write(3,*)"! 									        "
Write(3,*)"! SELEN is distributed in the /hope that it will be useful, but WITHOUT ANY  "
Write(3,*)"! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS  "
Write(3,*)"! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more      "
Write(3,*)"! details. 									"
Write(3,*)"!   										"
Write(3,*)"! You should have received a copy of the GNU General Public License along    "
Write(3,*)"! with SELEN.  If not, see <http://www.gnu.org/licenses/>.                   "
Write(3,*)"! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
Write(3,*) "!"
Write(3,*) "!------------------- 1- General settings --------------------------"  
Write(3,*) "!"
Write(3,*) "! --- Repository name"
Write(3,*) "CHARACTER*4, PARAMETER :: RUN=",trim(adjustl(runp))
Write(3,*) "! --- Tegmark resolution"
Write(3,*) "INTEGER, PARAMETER :: RES=", resolution
Write(3,*) "! --- Number of pixels"
Write(3,*) "INTEGER, PARAMETER :: NP=2*RES*(RES-1)*20+12"  
Write(3,*) "! --- Maximum harmonic degree"
Write(3,*) "INTEGER, PARAMETER :: LMAX=", degree   
!
Write(3,*) "! ********************************************** "
Write(3,*) "! --- Dealing with the degree 1 Love numbers --- "
Write(3,*) "! ********** [Updated GS Jan 29 2011] ********** "
Write(3,*) "! ********************************************** "
!
If (option_deg1=='y')    Write(3,*) "INTEGER, PARAMETER :: DEG1= 1"
If (option_deg1=='n')    Write(3,*) "INTEGER, PARAMETER :: DEG1= 0"
!
!If (option_deg1=='y') then 
If (option_rframe=='CM') Write(3,*) "CHARACTER*2, PARAMETER :: RFRAME= 'CM'" 
If (option_rframe=='CE') Write(3,*) "CHARACTER*2, PARAMETER :: RFRAME= 'CE'" 
!Endif
!
Write(3,*) "! --- Jmax index"
Write(3,*) "INTEGER, PARAMETER :: JMAX=(LMAX+1)*(LMAX+2)/2" 
Write(3,*) "! --- Pi"
Write(3,*) "REAL*8,  PARAMETER :: PI=3.14159265358979323840" 
Write(3,*) "! --- Water and ice densities (kg/m^3)"
!
Write(3,*) "REAL*4,  PARAMETER :: RHOW= ", RHO_WATER
Write(3,*) "REAL*4,  PARAMETER :: RHOI= ", RHO_ICE 
!
! -"Reference radius" and "reference gravity" (SI units)
Write(3,*) "REAL*4,  PARAMETER :: RAD_REF=6.371E6"
Write(3,*) "REAL*4,  PARAMETER :: GRA_REF=9.81E0 " 
!
Write(3,*) "! --- Type of numerical derivative used in SELEN"
If(option_der=='1') Write(3,*) "INTEGER, PARAMETER :: IDER= 1"
If(option_der=='2') Write(3,*) "INTEGER, PARAMETER :: IDER= 2"
!
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!------------------- 2- Ice models settings -----------------------"  
Write(3,*) "!"
Write(3,*) "! --- Ice sheet model file"
Write(3,*) "CHARACTER*",len_ice,", PARAMETER :: ICE_MODEL=", ice_filename
Write(3,*) "! --- Number of time steps"
Write(3,*) "INTEGER, PARAMETER :: NN=", NINC		
Write(3,*) "! --- Number of ice elements"
write(3,*) "INTEGER, PARAMETER :: NEL=", nice
Write(3,*) "! --- Time increment, ka"
Write(3,*) "REAL*4, PARAMETER :: DELTA=", DELTA   
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!"
!
If(option_lnap=='y') then
Write(3,*) "!------------------- 3- (ELASTIC - external Love numbers) settings ------"  
Write(3,*) "!"
Write(3,*) "! --- Number of viscous layers and fake code" 
Write(3,*) "INTEGER, PARAMETER :: CDE=", code
Endif
!
If    (option_nm=='y') then
Write(3,*) "!------------------- 3- Normal modes (TABOO) settings -------------"  
Write(3,*) "!"
Write(3,*) "! --- Earth model (see TABOO user guide)"
Write(3,*) "INTEGER, PARAMETER :: NV=", nv, ", CDE=", code
Write(3,*) "! --- Viscosity & Litho model file"
write(3,*) "CHARACTER*",LEN_VISCO,", PARAMETER :: VISCO_MODEL= &"
write(3,*)  trim(adjustl(visco_filename))
!
If(OPTION_TLOVE=='y') then 
Write(3,*) "INTEGER, PARAMETER :: TLOVE= 1"  
else
Write(3,*) "INTEGER, PARAMETER :: TLOVE= 0"  
endif
!
!If(OPTION_TLOVE_EXT=='y') then                       ! In PROGRESS === ! In PROGRESS === 
!Write(3,*) "INTEGER, PARAMETER :: TLOVE= 1"  
!else
!Write(3,*) "INTEGER, PARAMETER :: TLOVE= 0"  
!endif                                                ! In PROGRESS === ! In PROGRESS === 
!
elseif (option_pw=='y') then 
Write(3,*) "!------------------- 3- (ALMA - Option 2) settings -------------"  
Write(3,*) "!"
Write(3,*) "! --- Number of viscous layers and fake code" 
Write(3,*) "INTEGER, PARAMETER :: NV=", nv, ", CDE=", code
Write(3,*) "! --- Viscosity & Litho model file"
write(3,*) "CHARACTER*",LEN_VISCO,", PARAMETER :: VISCO_MODEL= & "
write(3,*)  trim(adjustl(visco_filename))
elseif (option_pwa=='y') then 
Write(3,*) "!------------------- 3- (ALMA - pre computed) settings -------------"  
Write(3,*) "!"
Write(3,*) "! --- Number of viscous layers and fake code" 
Write(3,*) "INTEGER, PARAMETER :: NV=", nv, ", CDE=", code
Write(3,*) "! --- Viscosity & Litho model file"
write(3,*) "CHARACTER*",LEN_VISCO,", PARAMETER :: VISCO_MODEL= & "
write(3,*)  trim(adjustl(visco_filename))
Endif
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!------------------- 4- Sea level Equation settings ---------------"  
Write(3,*) "!"
Write(3,*) "! --- Number of iterations"
Write(3,*) "INTEGER, PARAMETER :: SMAX=", iter_rev
Write(3,*) "! --- Mode of solution"
Write(3,*) "INTEGER, PARAMETER :: IMODE=", mode
If(option_topo=='y') then
Write(3,*) "! --- Topo file (for evolving shorelines)"
Write(3,*) "CHARACTER*",len_file_topo, ", PARAMETER :: TOPO_FILE=", "'"//trim(adjustl(file_topo))//"'" 
Write(3,*) "! --- Pixelized Topo file (for evolving shorelines)"
Write(3,*) "CHARACTER*",len_file_pxtopo, ", PARAMETER :: PXTOPO_FILE=", "'"//trim(adjustl(file_pxtopo))//"'" 
Endif 
!
!
If(option_rsl=='y'.or.option_rslp=='y')then 
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!------------------- 5- RSL analysis at site of a database --------"
Write(3,*) "!"
Write(3,*) "! --- RSL database"
Write(3,*) "CHARACTER*", len_rsl, ", PARAMETER :: RSL_DATABASE=", rsl_database
Write(3,*) "! --- Format of the RSL database"
Write(3,*) "CHARACTER* 1, PARAMETER :: RSL_DATABASE_FORMAT = ", & 
		"'"//trim(adjustl(RSL_DATABASE_FORMAT))//"'"
Write(3,*) "! --- Number of RSL sites in database"
Write(3,*) "INTEGER, PARAMETER :: NRSL=", nrsl
Endif
!
If(option_rslc=='y')then 
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!------------------- 6- Regional RSL analysis ---------------------"
Write(3,*) "!"
Write(3,*) "! --- Virtual RSL sites lon-lat file"
Write(3,*) "CHARACTER*", LEN_RSLC, ", PARAMETER :: RSLC_FILE=", RSLC_LONLAT_FILE
Write(3,*) "! --- Number of virtual RSL sites"
Write(3,*) "INTEGER, PARAMETER :: NRSLC=", nrslc
Write(3,*) "! --- Time of analysis (ka)"
Write(3,*) "REAL, PARAMETER :: TIME_BPC=", time_bpc
Endif
!
If(option_tg=='y'.or.option_tgplot=='y')then 
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!------------------- 7- Predictions at tide-gauges ----------------"
Write(3,*) "!"
Write(3,*) "! --- Tide-gauge database"
Write(3,*) "CHARACTER*", LEN_TGAUGES, ", PARAMETER :: TGAUGES_DATABASE=", TGAUGES_DATABASE	
!					      
Write(3,*) "! --- Tide-gauge database format"
Write(3,*) "CHARACTER* 1, PARAMETER :: TG_DATABASE_FORMAT = ", & 
		"'"//trim(adjustl(TGAUGES_DATABASE_FORMAT))//"'"    
!
Write(3,*) "! --- Number of tide gauges"
Write(3,*) "INTEGER, PARAMETER :: NTIDEGAUGES=", NTIDEGAUGES
Endif 
!
If(option_reb=='y') then 
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!--------------- 8- Elastic rebound maps ----------------"
Write(3,*) "!"
!
! -------- New as for June 30, 2012- 
If(option_xygrid=='1')     Write(3,*) "INTEGER, PARAMETER :: ELAREB_GRID=1"
If(option_xygrid=='2')     Write(3,*) "INTEGER, PARAMETER :: ELAREB_GRID=2"
! -------- New as for June 30, 2012- 
!
! Greenland, regional 
If  (option_reb_gr=='y') Write(3,*) "INTEGER, PARAMETER :: ELAREB_GR = 1"
If  (option_reb_gr=='n') Write(3,*) "INTEGER, PARAMETER :: ELAREB_GR = 0"
! Greenland, global 
If  (option_reb_gg=='y') Write(3,*) "INTEGER, PARAMETER :: ELAREB_GG = 1"
If  (option_reb_gg=='n') Write(3,*) "INTEGER, PARAMETER :: ELAREB_GG = 0"
! Antarctica, regional
If  (option_reb_ar=='y') Write(3,*) "INTEGER, PARAMETER :: ELAREB_AR = 1"
If  (option_reb_ar=='n') Write(3,*) "INTEGER, PARAMETER :: ELAREB_AR = 0"
! Antarctica, global
If  (option_reb_ag=='y') Write(3,*) "INTEGER, PARAMETER :: ELAREB_AG = 1"
If  (option_reb_ag=='n') Write(3,*) "INTEGER, PARAMETER :: ELAREB_AG = 0"
! Small glaciers, global
If  (option_reb_sm=='y') Write(3,*) "INTEGER, PARAMETER :: ELAREB_SM = 1"
If  (option_reb_sm=='n') Write(3,*) "INTEGER, PARAMETER :: ELAREB_SM = 0"
!
! Regional resolutions
Write(3,*) "REAL*4, PARAMETER :: RES_REG_GRC = ",trim(res_reg_grc)
Write(3,*) "REAL*4, PARAMETER :: RES_REG_ANC = ",trim(res_reg_anc)
!
Endif 
!
!
If(option_gm=='y') then 
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!------------------- 9- Global maps ----------------"
Write(3,*) "!"
!
! -------- New as for July 15, 2012- 
If(option_gmaps_xygrid=='1')     Write(3,*) "INTEGER, PARAMETER :: GMAPS_GRID=1"
If(option_gmaps_xygrid=='2')     Write(3,*) "INTEGER, PARAMETER :: GMAPS_GRID=2"
! -------- New as for July 15, 2012- 
!
If  (option_dots=='y') Write(3,*) "INTEGER, PARAMETER :: GLOB_I_S= 1"
If  (option_dots=='n') Write(3,*) "INTEGER, PARAMETER :: GLOB_I_S= 0"
If  (option_dotu=='y') Write(3,*) "INTEGER, PARAMETER :: GLOB_I_U= 1"
If  (option_dotu=='n') Write(3,*) "INTEGER, PARAMETER :: GLOB_I_U= 0"
If  (option_dotn=='y') Write(3,*) "INTEGER, PARAMETER :: GLOB_I_N= 1"
If  (option_dotn=='n') Write(3,*) "INTEGER, PARAMETER :: GLOB_I_N= 0"
If  (option_dotg=='y') Write(3,*) "INTEGER, PARAMETER :: GLOB_I_G= 1"
If  (option_dotg=='n') Write(3,*) "INTEGER, PARAMETER :: GLOB_I_G= 0"
If  (option_dotw=='y') Write(3,*) "INTEGER, PARAMETER :: GLOB_I_W= 1"
If  (option_dotw=='n') Write(3,*) "INTEGER, PARAMETER :: GLOB_I_W= 0"
!
If (option_dotfa=='y') Write(3,*) "INTEGER, PARAMETER :: GLOB_I_FA= 1"
If (option_dotfa=='n') Write(3,*) "INTEGER, PARAMETER :: GLOB_I_FA= 0"
If (option_dotss=='y') Write(3,*) "INTEGER, PARAMETER :: GLOB_I_SS= 1"
If (option_dotss=='n') Write(3,*) "INTEGER, PARAMETER :: GLOB_I_SS= 0"
!
If(option_dotloi=='y') Write(3,*) "INTEGER, PARAMETER :: GLOB_I_LOI= 1"
If(option_dotloi=='n') Write(3,*) "INTEGER, PARAMETER :: GLOB_I_LOI= 0"
If(option_dotloo=='y') Write(3,*) "INTEGER, PARAMETER :: GLOB_I_LOO= 1"
If(option_dotloo=='n') Write(3,*) "INTEGER, PARAMETER :: GLOB_I_LOO= 0"
If(option_dotlot=='y') Write(3,*) "INTEGER, PARAMETER :: GLOB_I_LOT= 1"
If(option_dotlot=='n') Write(3,*) "INTEGER, PARAMETER :: GLOB_I_LOT= 0"
!
Endif 
!


If(OPTION_ELA_MS=='y') then 
!
!
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!----------------- 9-bis Global maps and/or TG predictions for elastic steps ----------------"
Write(3,*) "!"
IF(OPTION_GIA_CORR=='y') Write(3,*) "INTEGER, PARAMETER :: GIA_CORR= 1"
IF(OPTION_GIA_CORR=='n') Write(3,*) "INTEGER, PARAMETER :: GIA_CORR= 0"
Write(3,*) "!"
If(OPTION_GM_ELA_MS_S_VAR=='y') Write(3,*) "INTEGER, PARAMETER :: GLOB_ELA_S_VAR= 1"
If(OPTION_GM_ELA_MS_S_VAR=='n') Write(3,*) "INTEGER, PARAMETER :: GLOB_ELA_S_VAR= 0"
If(OPTION_GM_ELA_MS_N_VAR=='y') Write(3,*) "INTEGER, PARAMETER :: GLOB_ELA_N_VAR= 1"
If(OPTION_GM_ELA_MS_N_VAR=='n') Write(3,*) "INTEGER, PARAMETER :: GLOB_ELA_N_VAR= 0"
If(OPTION_GM_ELA_MS_U_VAR=='y') Write(3,*) "INTEGER, PARAMETER :: GLOB_ELA_U_VAR= 1"
If(OPTION_GM_ELA_MS_U_VAR=='n') Write(3,*) "INTEGER, PARAMETER :: GLOB_ELA_U_VAR= 0"
If(OPTION_GM_ELA_MS_S_DOT=='y') Write(3,*) "INTEGER, PARAMETER :: GLOB_ELA_S_DOT= 1"
If(OPTION_GM_ELA_MS_S_DOT=='n') Write(3,*) "INTEGER, PARAMETER :: GLOB_ELA_S_DOT= 0"
If(OPTION_GM_ELA_MS_N_DOT=='y') Write(3,*) "INTEGER, PARAMETER :: GLOB_ELA_N_DOT= 1"
If(OPTION_GM_ELA_MS_N_DOT=='n') Write(3,*) "INTEGER, PARAMETER :: GLOB_ELA_N_DOT= 0"
If(OPTION_GM_ELA_MS_U_DOT=='y') Write(3,*) "INTEGER, PARAMETER :: GLOB_ELA_U_DOT= 1"
If(OPTION_GM_ELA_MS_U_DOT=='n') Write(3,*) "INTEGER, PARAMETER :: GLOB_ELA_U_DOT= 0"
!
IF(OPTION_TG_ELA_MS=='y') then 
!
Write(3,*) "!"
Write(3,*) "! --- Tide-gauge database"
Write(3,*) "CHARACTER*", LEN_FILE_TG_ELA_MS, ", PARAMETER :: TGAUGES_DATABASE=", DATABASE_TG_ELA_MS	
!					      
Write(3,*) "! --- Tide-gauge database format"
Write(3,*) "CHARACTER* 1, PARAMETER :: TG_DATABASE_FORMAT = ", & 
		"'"//trim(adjustl(ELAREB_MS_DATABASE_FORMAT))//"'"    
!
Write(3,*) "! --- Number of tide gauges"
Write(3,*) "INTEGER, PARAMETER :: NTIDEGAUGES=", NGEOD
!
ENDIF
!
Endif 







!
!
If(option_3d=='y')then 
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!------------------- 10- Predictions at geodetic sites ----------------"
Write(3,*) "!"
Write(3,*) "! "
Write(3,*) "! --- Database of geodetic points"
Write(3,*) "CHARACTER*", LEN_GEOD, ", PARAMETER :: GEODETIC_DATABASE=", & 
						   GEO_DATABASE
Write(3,*) "! --- Number of points"
Write(3,*) "INTEGER, PARAMETER :: NGEOD=", NGEOD
Endif 
!
If(option_3d_reb=='y')then 
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!--------- 10- Predictions at geodetic sites [Elastic rebound] ---------"
Write(3,*) "!"
Write(3,*) "! --- Database of geodetic points"
Write(3,*) "CHARACTER*", LEN_GEOD, ", PARAMETER :: GEODETIC_DATABASE=", &    
						   GEO_DATABASE_ELAREB
Write(3,*) "! --- Number of points"
Write(3,*) "INTEGER, PARAMETER :: NGEOD=", NGEOD
!
Write(3,*) "! --- Tide-gauge database format"
Write(3,*) "CHARACTER* 1, PARAMETER :: ELAREB_DATABASE_FORMAT = ", & 
		"'"//trim(adjustl(ELAREB_DATABASE_FORMAT))//"'"    
Endif 
!
If(option_3d_regions=='y')then 
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!----------------- 11- Regional 3D geodetic predictions --------------"
Write(3,*) "!"
Write(3,*) "! --- File with description of regions"
Write(3,*) "CHARACTER*", LEN_3D_REGIONS, ", PARAMETER :: TRED_REGIONS_DATABASE=", &
							 TRED_REGIONS_DATABASE
Write(3,*) "! --- Number of regions"
Write(3,*) "INTEGER, PARAMETER :: N_3D_REGIONS=", N_3D_REGIONS
Write(3,*) "! --- Names of regions"
Write(3,*) "CHARACTER*20 TRED_REGIONS_NAME(N_3D_REGIONS)"
do j=1, N_3D_REGIONS
        Write(3,*) "DATA TRED_REGIONS_NAME(", j, ")", "/"//trim(adjustl(TRED_REGIONS_NAME(j)))//"/" 
enddo
Endif
!
!
If(option_st=='y')then 
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!------------------ 12- Stokes coefficients (SC) settings ---------"
Write(3,*) "!"
Write(3,*) "! --- Min. degree for the SC"
Write(3,*) "INTEGER, PARAMETER :: STMIN=", degree_st_min
Write(3,*) "! --- Max. degree for the SC"
Write(3,*) "INTEGER, PARAMETER :: STMAX=", degree_st_max
endif 
!
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!"
!
! --- Closing "data.inc"
   close(3)
!
! --- Closing "selen.log"
   close(88)
!
STOP
END
!
!
!
!
!
!
 SUBROUTINE MAKE_OFDVMAP (DEGREE, FILE_GMT) 
 IMPLICIT NONE
 INTEGER NRESOLUTION 
 CHARACTER*3  DEGREE
 CHARACTER*20 FILE_GMT, PSFILENAME  
 CHARACTER*80 R_OPTION, B_OPTION, G_OPTION, & 
              W_OPTION, H_OPTION, S_OPTION, & 
	      J_OPTION, XY_OPTION
!
!
! # Revised GS June 2008 for upgrade to "Selen 2.6" 
!
!
  open(9,file=file_gmt,status='unknown')
!
  Write(9,*) "gmtset PAPER_MEDIA A4+"
  Write(9,*) "gmtset LABEL_FONT_SIZE 26p"
  Write(9,*) "gmtset ANOT_FONT_SIZE 16p"
  Write(9,*) " "
!
!
! ======= Settings for various options ======
!
! --- Projection ---
	   J_OPTION =  "-JX16l/16l"
!
! --- Range option ---
           R_OPTION = "-R0.5/300/1e-4/1e-1"
!
! --- Basemap option ---
	   B_OPTION = "-Bf3a2g3:'Harmonic degree':/f3a1g3p:'OF degree variance':WSne"
!
! --- Size of symbol --- 
	   S_OPTION = "-Ss0.45"
!
! --- Color of symbol --- 
	   G_OPTION = "-G0/0/255"
!
! --- Style of interpolating line --- 
	   W_OPTION = "-W8/255/0/0ta" 
!
! --- X-Y offset 
	   XY_OPTION = "-X4 -Y4" 
!
! --- Postscript output filename 
 	   psfilename = "ofdv.ps"
!
!
! ======= Writing the GMT script ======
!

! - Psbasemap 	
  Write(9,*) "psbasemap -U/-2/-2/'SELEN 3.2' -K", " ", & 
  	     trim(adjustl(XY_OPTION)), " ", & 
	     trim(adjustl(R_OPTION)),  " ", &    
	     trim(adjustl(B_OPTION)),  " ", &    
	     trim(adjustl(J_OPTION)),  " ", &
	     "> ", trim(adjustl(psfilename))  
!	   
! - A filter 
  Write(9,*) "awk '{print $1, $3}' ofdv.dat > predof.dat"
!
!
! - Plotting the interpolating line  
  Write(9,*) "psxy predof.dat -O -K -H3 -B -R -JX ", & 
  	     trim(adjustl(W_OPTION)), " ", & 
	     ">> ", trim(adjustl(psfilename))  
!
! - Plotting symbols 
  Write(9,*) "psxy ofdv.dat -O -K -H3 -B -R -JX ", & 
  	     trim(adjustl(S_OPTION)), " ", & 
  	     trim(adjustl(G_OPTION)), " ", & 
	     ">> ", trim(adjustl(psfilename))  
!
! - Plotting a string with harmonic degree	     
  Write(9,*) "echo '20 5e-2 18 0 0 BL MAX DEGREE=", trim(adjustl(DEGREE)), "'", & 
  	     " | pstext -N ", trim(adjustl(J_OPTION)), " ", & 
	                      trim(adjustl(R_OPTION)), " ", & 
			      " -G0 -W255 -O -K ", & 
	     		      ">> ", trim(adjustl(psfilename))  
!
! - Plotting a string with "power" of the interpolating line
  Write(9,*) "pstext -O dv-power.tmp -N -JX -R -G0 -W255 >> ", trim(adjustl(psfilename)) 
!
!
! - Cleaning 
  Write(9,*) "/bin/rm predof.dat"
!
 close(9)
!
 End Subroutine Make_ofdvmap
!
!
!
!
!
!
 SUBROUTINE MAKE_RSLC   & 
 	   (TIMEBP,     & 
	    MINRSL,     &
	    MAXRSL,     &
	    RSLINT,     &
	    LONMINC,    &
	    LONMAXC,    &
	    LATMINC,    &
	    LATMAXC,    & 
	    NV,         &
	    CODE,       &
	    RUN,        &
	    NINC,       &
	    NRSLC,      & 
	    TITLICE,    &
	    RESOLUTION, & 
	    ITER,       &
	    MODE,       &
	    DEGREE,     & 
	    VSTRING,    &
	    NAME_REGION,&
	    FILE_GMT,   &
	    SHORT)      
!
!     ----------------------------------------------------------------
! --- A GMT script for plotting RSL CONTOUR LINES - GS January 29 2008
!     ** Inspired by the work of Paolo Stocchi in the Mediterranean ** 
!     ---------------------------------------------------------------- 
!     Revided July 2006 for version 2.6 
!     Revised on April 2010 by GS for version 3.1 ALMA & g95 
!
 IMPLICIT NONE
 REAL*4 LAT_TIME, LON_TIME, LAT_STRING, MID_LON 
 REAL*4 FLAT_MIN, FLAT_MAX, FLON_MIN, FLON_MAX 
 REAL*4 RSL 
 INTEGER I, J, NN, NRSLC 
 CHARACTER*50 NAME_REGION
 CHARACTER*100 TITLE
 CHARACTER*10 TIMEBP, MINRSL, MAXRSL, RSLINT, LONMINC, LONMAXC, LATMINC, LATMAXC  
 CHARACTER*80 T_OPTION, R_OPTION, B_OPTION, G_OPTION, & 
              W_OPTION, C_OPTION, S_OPTION, H_OPTION, & 
	      D_OPTION, J_OPTION, W1_OPTION, W2_OPTION, ANNOT 
 CHARACTER*10 DEGREE, RESOLUTION, TITLICE
 CHARACTER*10 PROJ, COPT, GSPAC 
 CHARACTER*2  LABEL 
 CHARACTER*3  NINC, NV, CODE
 CHARACTER*1  ITER, MODE 
 CHARACTER*20 DATE, TIMC
 CHARACTER*20 FILE_GMT
 CHARACTER*30 VSTRING
 CHARACTER*4  RUN
 CHARACTER*100 SHORT
!	
  open(9,file=file_gmt,status='unknown')
  Write(9,*) "gmtset PAPER_MEDIA A4+"
  Write(9,*) "gmtset LABEL_FONT_SIZE 24p"
  Write(9,*) "gmtset ANOT_FONT_SIZE 10p"
  Write(9,*) "gmtset FRAME_WIDTH 0.15c"
  Write(9,*) " "
!
!
!
! ======= Settings for various options ======
!
! --- Range option ---
           R_OPTION = "-R"//trim(adjustl(lonminc))& 
                     //"/"//trim(adjustl(lonmaxc))& 
		     //"/"//trim(adjustl(latminc))&
		     //"/"//trim(adjustl(latmaxc))  
!
! --- Cpt option ---
           T_OPTION = "-T"//trim(adjustl(minrsl))&
                     //"/"//trim(adjustl(maxrsl))&
		     //"/"//trim(adjustl(RSLINT)) 
!
! --- Resolution of coastlines ---
           D_OPTION = "-Di"		     
!
! --- Shades for oceans and lands ---
	   S_OPTION = "-S0/60/255" 
           G_OPTION = "-G150/220/150" 
!
! --- Number of header lines in "rslc.dat" ---
           H_OPTION = "-H9"	   	
! 
! --- C option: a thick second level contour
	   C_OPTION="-C10"	
!
! --- Grid spacing for "Surface" ---
	   GSPAC="-I0.1"
!
! --- Projection ---
	   J_OPTION="-JM24"	 
!
! --- Annotation interval ---
           ANNOT='4'
!
! --- B Option ---
           title='"'//trim(adjustl(name_region))//'"'
	   B_OPTION="-Ba"//trim(adjustl(ANNOT))//"/a"&
	                 //trim(adjustl(ANNOT))//"WSEn"&
			 //":."//trim(adjustl(title))//":"
!
! --- Pen attributes for first level of contour lines ---	
           W1_OPTION = "-W3/255/0/0"
!
! --- Pen attributes for second level of contour lines ---	
           W2_OPTION = "-W10/255"
!
!
! ======= Writing the GMT script ======
!	   
! --- Surface ---  
    	   Write(9,*) "surface ", & 
	               trim(adjustl(GSPAC)), " ", & 
		       trim(adjustl(H_OPTION)), & 
		       " rslc.dat ", & 
		       trim(adjustl(R_OPTION)), & 
		       " -Gg.grd"
!
! --- A no-green palette ---
           Write(9,*) "makecpt -Cno_green ", & 
	   	       trim(adjustl(GSPAC)), " ", & 
		       trim(adjustl(T_OPTION)), & 
		       " > pal_rslc.cpt"
!
! --- Psbasemap ---
	   Write(9,*) "psbasemap -Y4 -X4 -R  ", & 
	   	       trim(adjustl(J_OPTION)), " ", & 
		       trim(adjustl(B_OPTION)), & 
		       " -P -K > rslc-map.ps"
!
! --- Coastlines ---
	   Write(9,*) "pscoast ", &
	   	       trim(adjustl(R_OPTION)), " ", &  
		       trim(adjustl(J_OPTION)), " ", & 
		       trim(adjustl(D_OPTION)), " ", & 
       	               " -B ", & 
		       trim(adjustl(S_OPTION)), " ", & 
		       trim(adjustl(G_OPTION)), " ", & 
		       " -O -K >> rslc-map.ps"		 
!
! --- First contour ---
	   Write(9,*) "grdcontour -Cpal_rslc.cpt ", & 
	               trim(adjustl(W1_OPTION)), & 
		       " -R  -G4/10 g.grd -JM -B -O -K >> rslc-map.ps"		       
!
! --- Second contour ---
	   Write(9,*) "grdcontour -U/0.5/0.5/'SELEN 3.2' ", & 
	               trim(adjustl(C_OPTION)), " ", & 
	               trim(adjustl(W2_OPTION)), & 
		       " -R -G4/4 -A1f10 g.grd -JM -B -O -K >> rslc-map.ps"
!
! --- Placing some text strings ---
           call CHAR10_2_FLOAT(lonminc, flon_min)
           call CHAR10_2_FLOAT(lonmaxc, flon_max)
 	   call CHAR10_2_FLOAT(latminc, flat_min)
 	   call CHAR10_2_FLOAT(latmaxc, flat_max)	  
	   mid_lon=flon_min+(flon_max-flon_min)/2. 	  
	   lat_string=flat_max+(flat_max-flat_min)/10.
	   lat_string=48.80 
	   lon_time=flon_max-(flon_max-flon_min)/30.
           lat_time=flat_min+(flat_max-flat_min)/10.  	   	   		       		      
!
	   If(code=='-1') then 
 	   open(4,file='rslc1.tmp',status='unknown') 
 	   Write(4,*) mid_lon, lat_string, "14 0 2 BC -Ice model: ", trim(adjustl(TITLICE)), &  	 
 	   " -LMAX=", trim(adjustl(DEGREE)),  " -RES=",  trim(adjustl(RESOLUTION)),& 	  
 	   " -ALMA rheology:", trim(adjustl(SHORT)), & 
 	   " -MODE=", trim(adjustl(MODE)),    " -ITER=", trim(adjustl(ITER)) 	   
 	   close(4) 
	   Else 
 	   open(4,file='rslc1.tmp',status='unknown') 
 	   Write(4,*) mid_lon, lat_string, "14 0 2 BC -Ice model: ", trim(adjustl(TITLICE)), & 
 	   " -Viscosity: ", trim(adjustl(VSTRING)),&
 	   " -LMAX=", trim(adjustl(DEGREE)),  " -RES=",  trim(adjustl(RESOLUTION)),& 	  
 	   " -NV=",   trim(adjustl(NV)),      " -CODE=", trim(adjustl(CODE)),& 
 	   " -MODE=", trim(adjustl(MODE)),    " -ITER=", trim(adjustl(ITER)) 	   
 	   close(4) 	   
	   Endif
!
 	   open(4,file='rslc2.tmp',status='unknown') 
	   Write(4,*) lon_time, lat_time, "16 0 0 BR ", trim(adjustl(timebp)), " ka"
	   close(4)	   	   
!
 	   Write(9,*) "pstext rslc1.tmp -N -JM -R -G0 -O -K >> rslc-map.ps" 
 	   Write(9,*) "pstext rslc2.tmp -N -JM -R -G0 -W255 -O >> rslc-map.ps" 

!
 End Subroutine MAKE_RSLC
!
!
!
!
!
!
	SUBROUTINE SCAN_REGION (FILEIN, FILEOUT, N, TIME, TIMEC, & 
		                LONMINC, LONMAXC, LATMINC, LATMAXC, & 
				MIN_RSLC, MAX_RSLC, CINT, & 
				NAME_OF_REGION)
	IMPLICIT NONE
!	
! # Scans "FILEIN" to retrieve information upon the number of points and 
!   time for the rsl contour analysis     === GS January 28, 2007 ===
!
	CHARACTER*10 LONMINC, LONMAXC, LATMINC, LATMAXC, TIMEC, MIN_RSLC, MAX_RSLC, CINT 
	CHARACTER*30 FILEIN, FILEOUT
	CHARACTER*50 NAME_OF_REGION
        CHARACTER*100 SS(10) 
        CHARACTER*200 LINE
	INTEGER, PARAMETER :: LARGE_INTEGER = 100
	INTEGER, PARAMETER :: MAXP = 29999
        INTEGER J, N, IERR, NOUT
	INTEGER ILAT, ILON, I_LON_MAX, I_LAT_MAX
	REAL*4 TIME, LON_LAT_INC, LON, LAT, LON_MIN, LON_MAX, LAT_MIN, LAT_MAX
	REAL*4, PARAMETER :: FRAME_WIDTH = 2. 
!
!
!
! ==> 1) Reading bounds from the input file
!
	IERR=0
!
	open(70,file=filein,status='old')
!
		do j=1, large_integer 
!
	        read(70,'(a200)',end=1) line 
!
		if(line(1:2)=='10') then 
				      CALL scan_string(line, 1, ss, nout)
				      name_of_region=trim(adjustl(ss(1)))
				 endif	
		if(line(1:2)=='20') then 
				      CALL scan_string(line, 1, ss, nout)
				      timec=trim(adjustl(ss(1)))
				      CALL CHAR100_2_REAL(ss(1), time)
				 endif
		if(line(1:2)=='30') then 
				      CALL scan_string(line, 2, ss, nout)		
			              CALL CHAR100_2_REAL(ss(1), lon_min)
			              CALL CHAR100_2_REAL(ss(2), lon_max)	
				      If(lon_max<lon_min) IERR=1  	
				      LONMINC=trim(adjustl(ss(1))) 
				      LONMAXC=trim(adjustl(ss(2))) 	
            			 endif				    
		if(line(1:2)=='40') then 
				      CALL scan_string(line, 2, ss, nout)
				      CALL CHAR100_2_REAL(ss(1), lat_min)
			              CALL CHAR100_2_REAL(ss(2), lat_max)
				      If(lat_max<lat_min) IERR=1  
				      LATMINC=trim(adjustl(ss(1))) 
				      LATMAXC=trim(adjustl(ss(2))) 				 			      
				 endif
		if(line(1:2)=='50') then 
				      CALL scan_string(line, 1, ss, nout)
				      CALL CHAR100_2_REAL(ss(1), lon_lat_inc)
				      If(lon_lat_inc<=0.)IERR=1			
				 endif 
		if(line(1:2)=='60') then 
				      CALL scan_string(line, 3, ss, nout)
				      MIN_RSLC=trim(adjustl(ss(1)))
				      MAX_RSLC=trim(adjustl(ss(2)))
	 			      CINT=trim(adjustl(ss(3)))
				 endif 								 				 		  
		enddo
1       continue
        close(70)
!
	If(IERR==1) then 
        	Write(*, *) "File ", trim(adjustl(filein)), " is badly configured"
        	Write(88,*) "File ", trim(adjustl(filein)), " is badly configured"
        	Call Stop_Config 
	Endif	
	
!
! ==> 2) Computing lon-lat of points and writing on the outpout file... 
!
! --- Counts the number of points in the User-supplied file "filein" 

	
! NEW code as of November 2009 
	i_lon_max=aint((lon_max-lon_min)/lon_lat_inc) + 1
	i_lat_max=aint((lat_max-lat_min)/lon_lat_inc) + 1
	n=0 
	do 2 ilon=1, i_lon_max 
		do 2 ilat=1, i_lat_max 
		n=n+1	
2	continue 
        Write(88,*) "The number of points in file ", trim(adjustl(filein)), " is", n
!	
!	n=0
!	do 2 lon = lon_min, lon_max, lon_lat_inc 
!		do 2 lat = lat_min, lat_max, lon_lat_inc 
!		n=n+1	
!2       continue
!        Write(88,*) "The number of points in file ", trim(adjustl(filein)), " is", n
!	
! --- Adds a frame to avoid border effects using GMT  
	lon_min=lon_min - frame_width 
	lon_max=lon_max + frame_width 
	lat_min=lat_min - frame_width
	lat_max=lat_max + frame_width
!
! NEW code as of November 2009 
	i_lon_max=aint((lon_max-lon_min)/lon_lat_inc) + 1
	i_lat_max=aint((lat_max-lat_min)/lon_lat_inc) + 1
	n=0 
	do 3 ilon=1, i_lon_max 
		do 3 ilat=1, i_lat_max 
		n=n+1	
3	continue 
        Write(88,*) "After framing (2 degrees NSEW), the number of points is ", n	
!
!	n=0
!	do 3 lon = lon_min, lon_max, lon_lat_inc 
!		do 3 lat = lat_min, lat_max, lon_lat_inc 
!		n=n+1	
!3       continue
!        Write(88,*) "After framing (2 degrees NSEW), the number of points is ", n	
!
	If(n>=maxp) then 
        	Write(*, *) "The total number of points in file ", trim(adjustl(filein)), " is", n
		Write(*, *) "this exceeds the maximum allowed (i. e.,", maxp, ")" 
        	Write(88,*) "The total number of points in file ", trim(adjustl(filein)), " is", n 
		Write(88,*) "this exceeds the maximum allowed (i. e.,", maxp, ")" 
        	Call Stop_Config 
	Endif	
!
! --- Composes the final set of points and reports on "fileout"
!     Negative lon values are transformed in the range [0:360] 
!
        open(71,file=fileout,status='unknown') 	
	i_lon_max=aint((lon_max-lon_min)/lon_lat_inc) + 1
	i_lat_max=aint((lat_max-lat_min)/lon_lat_inc) + 1
!
	do 4 ilon=1, i_lon_max 
		lon=lon_min+(ilon-1)*lon_lat_inc		
!
		do 4 ilat=1, i_lat_max	
		lat=lat_min+(ilat-1)*lon_lat_inc							
!
			if(lon>=0.) then 
				write(71,*) lon, lat 
				else 
				write(71,*) lon+360., lat 		
			Endif		
4       continue
	close(71) 
!
!        open(71,file=fileout,status='unknown') 	
!	 do 4 lon = lon_min, lon_max, lon_lat_inc 
!		do 4 lat = lat_min, lat_max, lon_lat_inc 				
!			if(lon>=0.) then 
!				write(71,*) lon, lat 
!				else 
!				write(71,*) lon+360., lat 		
!			Endif		
!4       continue
!	 close(71) 
! 	
	End subroutine scan_region
!
!
!
!
!
!
   	SUBROUTINE MAKE_ESLPLOT_PM (TITLE, FILE_GMT)
!
        IMPLICIT NONE
!
!--- Creates a GMT script for plotting the Equivalent Sea Level as a function
!    of time, for the ice model with name "titlice".  == GS Urbino 8/12/07 ==
!    RE-DESIGNED by on MARCH 02, 2012 for the "pm" ice models...
!
	CHARACTER*20 FILE_GMT 
	CHARACTER*10 TITLE
	INTEGER K
	
	CHARACTER*30, PARAMETER:: FILE_IN1="esl-pm.dat"
	CHARACTER*30, PARAMETER:: FILE_IN2="esl-pm-thin.dat"
	CHARACTER*30, PARAMETER:: FILE_OUT="esl-pm.ps"
	CHARACTER*80, PARAMETER:: B_OPTION="-Ba10f5:'time (years)':/a0.1f0.1:'Equivalent Sea Level (ESL, m)':SWen"
	CHARACTER*30, PARAMETER:: R_OPTION="-R-5/205/-0.1/1"
	CHARACTER*30, PARAMETER:: J_OPTION="-JX17.33/7.4164"
        CHARACTER*80, PARAMETER :: XY_OPTION1="-X4 -Y4"
        CHARACTER*80, PARAMETER :: XY_OPTION2="-X0 -Y0"
        CHARACTER*80, PARAMETER :: W_OPTION1=" -W4  "
        CHARACTER*80, PARAMETER :: W_OPTION2=" -W1  "  


!
!
        open(9,file=file_gmt,status='unknown')
        Write(9,*) " gmtset PAPER_MEDIA A4+"
        Write(9,*) " gmtset LABEL_FONT_SIZE 16p"
        Write(9,*) " gmtset ANOT_FONT_SIZE 12p"
	Write(9,*) " gmtset FRAME_WIDTH 0.01c"
	Write(9,*) " gmtset TICK_LENGTH -0.2c"
	Write(9,*) " gmtset FRAME_PEN = 0.5p"
        Write(9,*) " "
!
        Write(9,*) " "
        Write(9,*) "FILE_IN1=",  '"'//trim(adjustl(file_in1))//'"'          
        Write(9,*) "FILE_IN2=",  '"'//trim(adjustl(file_in2))//'"'          
        Write(9,*) "FILE_OUT=",  '"'//trim(adjustl(file_out))//'"'       
        Write(9,*) " "
        Write(9,*) "psxy ", trim(adjustl(XY_OPTION1)), " ", trim(adjustl(R_OPTION)),  " ", "$FILE_IN1 ", & 
	                    trim(adjustl(W_OPTION1)),  " ", trim(adjustl(J_OPTION)), " ", & 
			    trim(adjustl(B_OPTION)),    " -K > $FILE_OUT "
        Write(9,*) " "
        Write(9,*) "psxy ", trim(adjustl(XY_OPTION2)), " ", trim(adjustl(R_OPTION)),  " ", "$FILE_IN2 ", & 
	                    trim(adjustl(W_OPTION2)),  " ", trim(adjustl(J_OPTION)), " ", & 
			    trim(adjustl(B_OPTION)),    " -O >> $FILE_OUT "
!
         close(9) 
!
   	 End Subroutine make_eslplot_pm
!
!
!
!
!
!
!
!
!
!
!
!
   	SUBROUTINE MAKE_ESLPLOT (NINC, TITLICE, SHOF_FILE, FILE_GMT)
!
        IMPLICIT NONE
!
!--- Creates a GMT script for plotting the Equivalent Sea Level as a function
!    of time, for the ice model with name "titlice".  == GS Urbino 8/12/07 ==
!
!    Revised July 4 2008 for v. 2.6 
!    Also revised July 8
!    *** Revised GS & GG May 2011 - Implementation of the Antarctic model "ANTA"
!    *** Revised GS & GG June 2011 - Implementation of the "ANTA" ice model...
! :: Revised GS (RECTANGULAR ICE for the benchmark SLE II ) <<<<== IN PROGRESS 
! // Revised GS on february 2012  for the "GRD" ice eimplementation    
!!   Jul 1, 2013; The "FLO" ice for coupling with her models 

        CHARACTER*100 B_OPTION
	CHARACTER*40 Y_LABEL, R_OPTION, ESL_STRING
	CHARACTER*20 FILE_GMT, SHOF_FILE
	CHARACTER*20 AWKSTRING
	CHARACTER*10 TITLICE
	CHARACTER*3  NINC
	REAL VMAX, VMIN, ESL
	INTEGER K
!
!
        open(9,file=file_gmt,status='unknown')
        Write(9,*) "gmtset PAPER_MEDIA A4+"
        Write(9,*) "gmtset LABEL_FONT_SIZE 24p"
        Write(9,*) "gmtset ANOT_FONT_SIZE 12p"
        Write(9,*) " "
!
!
! --- R- option 
        if(titlice=='ICE1')     R_OPTION="-R-2/23/0/220"
        if(titlice=='ICE3G')    R_OPTION="-R-2/23/0/220"
        if(titlice=='IMED1')    R_OPTION="-R-2/23/0/220"
        if(titlice=='DISK')     R_OPTION="-R-2/23/0/220"
        if(titlice=='ICE5G')    R_OPTION="-R-2/23/0/220" 
        if(titlice=='ICE5G26')  R_OPTION="-R-2/28/0/220" 
        if(titlice=='IJ05MOD')  R_OPTION="-R-2/23/0/220" 
        if(titlice=='ANU05')    R_OPTION="-R-2/32/0/220" 
        if(titlice=='ICAP')     R_OPTION="-R-2/23/0/220"
        if(titlice=='ALPS')     R_OPTION="-R-2/23/0/2" 
        if(titlice=='FLOR')     R_OPTION="-R-2/44/0/220"  
        if(titlice=='GREEN')    R_OPTION="-R-2/23/0/1" 
        if(titlice=='GLAC')     R_OPTION="-R-2/23/0/1" 	
        if(titlice=='ANTA')     R_OPTION="-R-2/23/0/1" 
        if(titlice=='TRANS')    R_OPTION="-R-2/50/0/10" 
        if(titlice=='RETT')     R_OPTION="-R-2/50/0/10" 
        if(titlice=='GRD')      R_OPTION="-R-2/11/-1/1"
        if(titlice=='BENO')      R_OPTION="-R-2/201/-1/1"  
!
! --- y label 
	if(titlice=='ICE1')     y_label='ICE1 ESL (m)'
	if(titlice=='ICE3G')    y_label='ICE3G ESL (m)'
	if(titlice=='IMED1')    y_label='IMED1 ESL (m)'
	if(titlice=='DISK')     y_label='DISK ESL (m)'
	if(titlice=='IJ05MOD')  y_label='IJ05MOD ESL (m)'
	if(titlice=='ICE5G')    y_label='ICE5G ESL (m)'
	if(titlice=='ICE5G26')  y_label='ICE5G ESL (m)'
        if(titlice=='ANU05')    y_label='ANU ESL (m)'
	if(titlice=='ALPS')     y_label='ALPS ESL (m)'
	if(titlice=='FLOR')     y_label='FLO ESL (m)'
        if(titlice=='ICAP')     y_label='ICAP ESL (m)'
        if(titlice=='GREEN')    y_label='GREEN ESL (m)'
        if(titlice=='ANTA')     y_label='ANTA ESL (m)'
        if(titlice=='GLAC')     y_label='GLAC ESL (m)'
        if(titlice=='TRANS')    y_label='TRANS ESL (m)'
        if(titlice=='RETT')     y_label='RECT ESL (m)'
        if(titlice=='GRD')      y_label='GRD ESL (m)' 
        if(titlice=='BENO')      y_label='LEGOS ESL (m)' 

!
! --- B option 
	                    b_option=" -Ba2f1:'time (ka)':/a50f10WSen:"//"'"//trim(adjustl(y_label))//"':"

        if(titlice=='GRD')  b_option=" -Ba2f1:'time (*10 years)':/a0.1f0.1WSen:"//"'"//trim(adjustl(y_label))//"':"
	
! --- A psbasemap 
        Write(9,*) "psbasemap -U'SELEN 3.2' -X6 -Y10 ", " ", & 
	            trim(adjustl(b_option)), " ", & 
		    trim(adjustl(R_OPTION)), " -JX14/9 -K > plot.ps"
 
! --- Draw thick and thin lines 
	if(ninc=='30') awkstring="'{print 30.-$1, $2}'"
	if(ninc=='26') awkstring="'{print 26.-$1, $2}'"
	if(ninc=='21') awkstring="'{print 21.-$1, $2}'"
	if(ninc=='21') awkstring="'{print 24.-$1, $2}'"
	if(ninc=='18') awkstring="'{print 18.-$1, $2}'"
	if(ninc=='1' ) awkstring="'{print  1.-$1, $2}'"
	if(ninc=='8' ) awkstring="'{print  8.-$1, $2}'"
	if(ninc=='50') awkstring="'{print 50.-$1, $2}'"
	if(ninc=='41') awkstring="'{print 41.-$1, $2}'"  
	if(ninc=='20') awkstring="'{print 20.-$1, $2}'"
	if(ninc=='10') awkstring="'{print 10.-$1, $2}'"
	if(ninc=='200')awkstring="'{print 200.-$1, $2}'"


!	
	Write(9,*) "awk "//trim(adjustl(awkstring))//" esl-thin.dat > esl.tmp"	
	Write(9,*) "psxy esl.tmp -M -B -R -JX -W2/0 -O -K >> plot.ps"
!
	Write(9,*) "awk "//trim(adjustl(awkstring))//" esl.dat > esl.tmp"	
	Write(9,*) "psxy esl.tmp -M -B -R -JX -W6/0 -O -K >> plot.ps"
!
! --- Inset with total Equivalent Sea level 
	 Write(9,*) "pstext -N esl-tot.dat -JX -R -O >> plot.ps"
!
! --- The end 
         Write(9,*) "mv plot.ps esl.ps"
!
         close(9) 
!
   	 End Subroutine make_eslplot 
!
!
!
!
!
!
	SUBROUTINE OF_PIXELIZATION_REAL
        IMPLICIT NONE
	CHARACTER*20, PARAMETER :: FILENAME='px.gmt'
!				        
!--- REALISTIC OCEAN FUNCTION (according to GMT topography database)   
!
!--- Creates a GMT script named "px.gmt", which is used to separate wet from dry pixels. 
!    The wet/dry distribution of pixels *is* sensitive to the GMT options -D (coastlines 
!    resolution) and -A (spatial filter), so different choices may  give (slightly) 
!    different distributions, affecting all the computations in SELEN. 
!				        
!    [*** GS & FC Bologna 8/12/07 ***]
!    Last change: GS February 21, 2008
!    Name changed on July 19 for v. 2.6 
!
        open(13,file=filename,status='unknown')
!Write(13,*)"echo '     - ", trim(adjustl(filename))//":", " gmt-selecting pixels in <<pxa.dat>>'"
	Write(13,*)"#"
	Write(13,*)"# GMT script for sepating wet from dry pixels: px.gmt"
	Write(13,*)"#               REALISTIC OCEAN FUNCTION             " 
	Write(13,*)"# ====== Made by GS and FC on December 8, 2007 ======"	
	Write(13,*)"#" 

!
! Until June 2014 (VB December 2011)
!Write(13,*)"gmtselect pxa.dat -H4 -Di -R0/360/-90/90  -A0 -JQ180/200 -Nk/s/k/s/k > weta.dat"
!Write(13,*)"gmtselect pxa.dat -H4 -Di -R              -A0 -JQ        -Ns/k/s/k/s > drya.dat"

! Tests for the new ice model 
Write(13,*)"gmtselect pxa.dat -H4 -Dc -R0/360/-90/90  -A0 -JQ180/200 -Nk/s/s/s/s > weta.dat"
Write(13,*)"gmtselect pxa.dat -H4 -Dc -R              -A0 -JQ        -Ns/k/k/k/k > drya.dat"



! To reproduce the previous GJI computations- GS June 5, 2012 [AND install the shorelines...]
! Write(13,*)"gmtselect pxa.dat -H4 -Df -R0/360/-90/90  -A0 -JQ180/200 -Nk/s/k/s/k > weta.dat"
! Write(13,*)"gmtselect pxa.dat -H4 -Df -R              -A0 -JQ        -Ns/k/s/k/s > drya.dat"

!	Write(13,*)"gmtselect pxa.dat -H4 -Di -R0/360/-90/90  -A0 -JQ0/16 -Nk/s/s/s/s > weta.dat"
!	Write(13,*)"gmtselect pxa.dat -H4 -Di -R              -A0 -JQ     -Ns/k/k/k/k > drya.dat"

	close(13) 		
	End Subroutine OF_PIXELIZATION_REAL
!
!
!
!
!
!
	SUBROUTINE OF_PIXELIZATION_ZONAL(FILECAP,AMPC)
        IMPLICIT NONE
	CHARACTER*20, PARAMETER :: FILENAME='px.gmt'
	CHARACTER*20 FILECAP
	CHARACTER*10 AMPC 
	REAL*4 AMP, LATCAP 
!
!				        
!--- "ZONAL" OCEAN FUNCTION    
!
!--- Creates a GMT script named "px.gmt", which is used to separate wet from dry pixels. 
!    The wet/dry distribution of pixels *is* sensitive to the GMT options -D (coastlines 
!    resolution) and -A (spatial filter), so different choices may  give (slightly) 
!    different distributions, affecting all the computations in SELEN. 
!				        
!    Created by GS on July 19, 2008 for v. 2.6
!
!
!
!#------ Creating a 'cap' multi-segment file where a 
!        circular polygon of radius (colatitude) AMP
!
	Call CHAR10_2_REAL(ampc, amp)
	If(amp.le.0.or.amp.ge.180.)then 
          Write(88,*) "The size of the Zonal Ocean seems to be out of bounds..." 
          Write(*, *) "The size of the Zonal Ocean seems to be out of bounds..." 
	  call Stop_Config	
	Endif				
!	
	Open(12,file=filecap,status='unknown')
	  latcap=90.-amp 
	  Write(12,'(a1)')">"
	  Write(12,'(a6)')"0   90"  
 	  Write(12,'(a3,1x,f14.8)')"0  ", latcap  
 	  Write(12,'(a3,1x,f14.8)')"90 ", latcap  
 	  Write(12,'(a3,1x,f14.8)')"180", latcap  
 	  Write(12,'(a3,1x,f14.8)')"270", latcap  
 	  Write(12,'(a3,1x,f14.8)')"360", latcap  
	  Write(12,'(a6)')"360 90"
	Close(12)  
!
!
!#------ Creating a GMT script for separating wet from dry pixels... 
!
        Open(13,file=filename,status='unknown')
          Write(13,*)"echo '     - ", trim(adjustl(filename))//":", " gmt-selecting pixels in <<pxa.dat>>'"
	  Write(13,*)"#"
	  Write(13,*)"#GMT script for sepating wet from dry pixels: px.gmt"
	  Write(13,*)"#               *ZONAL* OCEAN FUNCTION              " 
	  Write(13,*)"# ============ Made by GS July 19, 2008 ============"
	  Write(13,*)"#" 
	  Write(13,*)"gmtselect -F"//trim(adjustl(filecap)), " pxa.dat -H4 -JQ180/200 -If > weta.dat"
	  Write(13,*)"gmtselect -F"//trim(adjustl(filecap)), " pxa.dat -H4 -JQ            > drya.dat"
	Close(13) 
!		
	End Subroutine OF_PIXELIZATION_ZONAL
!
!
!
!   	
!
!
        SUBROUTINE MAKE_RAINBOW_PALETTE 
        IMPLICIT NONE
!
!--- Creates a "rainbow palette", named "ice-pal.cpt" suitable
!    for plotting both the original and the reconstructed ice 
!    distribution...     *** GS and FC Bologna 7/12/07 ***
!
        open(13,file='ice-pal.cpt',status='unknown') 
	Write(13, '(a1)') "#"
	Write(13,'(a74)') "# A rainbow palette by Florence and Giorgio as of 7 Dec. 2007: ice-pal.cpt"
	Write(13, '(a1)') "#"	
	Write(13,'(a58)') "-500   255     0       255     0       255     0       255"
	Write(13,'(a58)') "0      164     0       255     500     164     0       255"	
	Write(13,'(a58)') "500    73      0       255     1000    73      0       255"
	Write(13,'(a58)') "1000   0      18       255     1500    0      18       255"
	Write(13,'(a58)') "1500   0     109       255     2000    0     109       255"
	Write(13,'(a58)') "2000   0     201       255     2500    0     201       255"
	Write(13,'(a58)') "2500   0     255       219     3000    0     255       219"
	Write(13,'(a58)') "3000   0     255       127     3500    0     255       127"
	Write(13,'(a57)') "3500   0     255        36     4000    0     255       36"
	Write(13,'(a56)') "4000   54    255         0     4500    54    255       0"
	Write(13,'(a56)') "4500   146   255         0     5000    146   255       0"
	Write(13,'(a56)') "5000   237   255         0     5500    237   255       0"
	Write(13,'(a56)') "5500   255   182         0     6000    255   182       0"
	Write(13,'(a56)') "6000   255   91          0     6500    255   91        0"
	Write(13,'(a17)') "B  255    0   255"
	Write(13,'(a17)') "F  255    0     0"
	Write(13,'(a17)') "N  128  128   128"	    
        close(13) 	
!
 	END SUBROUTINE MAKE_RAINBOW_PALETTE 
!
!
!
!
!
!
!
 SUBROUTINE MAKE_TGAUGESSCA (DATAFILE, FILE1_GMT)
 IMPLICIT NONE
!
! --- This routine perfoms simple statistical analyses on tide-gauge data 
!     and creates a two-frames scatterplot of tide gauges trends vs years 
!     of data 
!
!     *** Last modified GS December 8, 2007 *** 
!
!     === Re touched on July 2008 for SLEN 2.6 ===
! 
!
 CHARACTER*30 DATAFILE
 CHARACTER*20 FILE1_GMT, FILE2_GMT 
 INTEGER I, J, K, N, P, NH, IMIN, IMAX, PMIN, PMAX, IEARS_MAX, IEARS_MIN
 INTEGER, PARAMETER :: ILARGE=10000 ! A large number 
 INTEGER IEARS(ILARGE) 
 CHARACTER*100 RIGA		   ! A row 
 CHARACTER*8 PCODE                 ! Tide gauge code  
 CHARACTER*3 GLOSS_CODE            ! GLOSS code 
 CHARACTER*3 YEARS                 ! Number of years 
 CHARACTER*11 RANGE_YEARS	   ! Range of years 
 CHARACTER*9 DATUM		   ! Trend 
 CHARACTER*3 PLUS_MINUS 	   ! +/- 
 CHARACTER*5 ERROR		   ! Error 
 CHARACTER*7 STDV		   ! Std. deviation of residues 						       
 CHARACTER*2 LAT1, LAT2 	   ! Latitide degrees and minutes      
 CHARACTER*3 LON1		   ! Longitude degrees 
 CHARACTER*2 LON2		   ! Longitude minutes 
 CHARACTER*30 NAME(ILARGE)	   ! Name of the station 
 CHARACTER*1 WLAT, WLON 	   ! Where is the station (e, w, s, n) 
 CHARACTER*10 YMINC, YMAXC, DELTAC, IMINC, IMAXC, PMINC, PMAXC 
 REAL*8 LAT(ILARGE), NLAT, DLAT 	   ! Latitude  
 REAL*8 LON(ILARGE), NLON, DLON 	   ! Longitude 
 REAL*8 TREND(ILARGE), D_TREND(ILARGE)
 REAL*8 AVE(ILARGE), DAVE(ILARGE), SUMP(ILARGE)	
 REAL*8 AVEW, DAVEW, SUMPW	
 REAL*8 RANGE, DELTA
 REAL*8 AVE_RANGE, AVE_MAX, AVE_MIN, AVE_DELTA
 REAL*8 TREND_MAX, TREND_MIN 
 CHARACTER*100  B_OPTION, R_OPTION  
 CHARACTER*80 STRING1, STRING2
	  
!
!
! ===============================================
! ====== Analysis of tide gauges database =======
! ===============================================
!
! --- Open and reads the tide-gauges data file for getting info...
!
! --- Counting the header lines (beginning with '#'), and the data lines
       open(10,file=DATAFILE,status='unknown')
       nh=0 ; n=0
       do i=1, ilarge 
       		read(10,'(a100)',end=1) riga
       		if(riga(1:1)=='#') nh=nh+1
       		if(riga(1:1)/='#') n=n+1
       enddo 
 1     close(10)  
!
! --- Opening the tide-gauges data file... 
       open(10,file=DATAFILE,status='unknown')
!
! --- Reading again the header lines  
       do i=1, nh ; read(10,'(a100)',end=1) riga ; enddo   
!
! --- Loop on the tide-gauge stations 
       do 5 i=1, n
! 
! --- Reading one line   
       	     read  (10,200) pcode, gloss_code, years, range_years, & 
   		     	    datum , plus_minus, error, stdv, & 
			    lat1, lat2, wlat, lon1, lon2, wlon, name(i)  
! --- Reading format 
200    format(a8,1x,a3,1x,a3,2x,a11,1x,a8,1x, & 
              a3,1x,a5,1x,a7,3x,a2,1x,a2,1x,a1,1x,& 
              a3,1x,a2,1x,a1,3x,a30)			      
!
! --- Converts the range in floating-point format 
	     open(50,file='junk-tgauges.tmp',status='unknown') ; write(50,'(a3)') years ; close(50) 
	     open(50,file='junk-tgauges.tmp',status='unknown') ; read (50,*)      iears(i) ; close(50)
!
! --- Converts the trend in floating-point format 
	     open(50,file='junk-tgauges.tmp',status='unknown') ; write(50,'(a9)') datum ; close(50) 
	     open(50,file='junk-tgauges.tmp',status='unknown') ; read (50,*)      trend(i) ; close(50)	     
!
! --- Converts the error on trend in floating-point format 
	     open(50,file='junk-tgauges.tmp',status='unknown') ; write(50,'(a5)') error ; close(50) 
	     open(50,file='junk-tgauges.tmp',status='unknown') ; read (50,*)    d_trend(i) ; close(50)	     	     	  
!
! --- Extracts lon and lat of the station from the texst strings... 
             call  find_lon_lat (lat1, lat2, wlat, lon1, lon2, wlon, lon(i), lat(i)) 
!
! --- End of Loop
 5 continue 	
!
      close(10)		
!
!
! ==========================================
! --- Filter of data for stats and plots ---
! ==========================================
!
! --- Target file for stats
!
        open(73,file='tgauges-stat.dat',status='unknown')
	Write(73,*) " "    
!
! >>> #1 Creating a scatterplot file  
!
  	Open(28,file='tgauges-scpl.dat',status='unknown') 
!
	do i=1, n 
		write(28,*) iears(i), trend(i), d_trend(i), lon(i), lat(i), name(i)  
	enddo 
	Close(28) 
!
! >>> #2 Weighted average of the whole set of data  
!
	avew  = 0. 
	davew = 0. 
	sumpw = 0.
	do i=1, n      
		avew  = avew  + trend(i)   *float(iears(i))
		davew = davew + d_trend(i) *float(iears(i)) 
		sumpw = sumpw +             float(iears(i))
	enddo
	avew  = avew/sumpw 
	davew = davew/sumpw
	Write(73,*) "Average weighted rate of sea level change: "
	Write(73,'(f5.2,a1,a3,a1,f5.2,a6)') & 
		    avew, " ", plus_minus, " ", davew, ' mm/yr'
	Write(73,*) " "
!
! >>> #3 Average over subsets of increasing number of years 
!
	iears_max=-99999
	iears_min=+99999
	do i=1, n 
		if(iears(i)>=iears_max)iears_max=iears(i) 
		if(iears(i)<=iears_min)iears_min=iears(i) 		
	enddo
	Write(73,*) 'Max number of years of data = ', iears_max 
	Write(73,*) 'Min number of years of data = ', iears_min 
	Write(73,*) " "	
!
!
! -------------------------------
! --- GMT scripts for Figures ---
! -------------------------------
!	
! === Figure 1: A two-frame figure with scattered tide-gauges data
!		with and w/o error bars, with the weighted average 
!		based on all data. Useful for presentations. 
!
! --- Determines min and max trend of tide-gauge data
	trend_max=-99999
	trend_min=+99999
	do i=1, n 
		if(trend(i)>=trend_max)then
					trend_max=trend(i) 
					imax=i 
		Endif
		if(trend(i)<=trend_min)then 
					trend_min=trend(i) 
					imin=i
	        Endif			
	enddo
	Write(73,*) 'Max trend (mm/yr) = ', trend_max 
	Write(73,*) 'Min trend (mm/yr) = ', trend_min
	Write(73,*) " "	
	close(73) 
!
	trend_max=trend_max+d_trend(imax)
	trend_min=trend_min-d_trend(imin)	
!
	range = trend_max - trend_min 
!	
	if(range>=0  .and.range<50)  delta=10  
	if(range>=50 .and.range<100) delta=10	
	if(range>=100.and.range<200) delta=20
	if(range>200) 		  delta=50 
!
   	OPEN  (10,FILE='junk.dat',STATUS='unknown')
	WRITE (10,'(i10)') int(trend_min) ; WRITE (10,'(i10)') int(trend_max)
	WRITE (10,'(i10)') int(delta) 
	WRITE (10,'(i10)') int(iears_min) ; WRITE (10,'(i10)') int(iears_max)	
	CLOSE(10)

   	OPEN  (10,FILE='junk.dat',STATUS='unknown')
	READ (10,'(a10)') yminc ; READ (10,'(a10)') ymaxc
	READ (10,'(a10)') deltac
	READ (10,'(a10)') iminc ; READ (10,'(a10)') imaxc 	
	CLOSE(10)
	iminc='0'
!
! --- Target fiel for GMT script 
 	 open (19,file=file1_gmt,status='unknown')
	 Write(19,*) "gmtset PAPER_MEDIA A4+" 
	 Write(19,*) "gmtset HEADER_FONT_SIZE 24p"
	 Write(19,*) "gmtset FRAME_WIDTH 0.1c"
	 Write(19,*) " "
!
! --- -R and -B options 
	 r_option = "-R"//trim(adjustl(iminc))//"/"//trim(adjustl(imaxc))//"/"//trim(adjustl(yminc))//"/"//trim(adjustl(ymaxc))
	 b_option = "-Ba30f10:'Lenght of record (yr)':/a"//trim(adjustl(deltac))//"f"//trim(adjustl(deltac))//"WSen:'Rate of SLC, mm/yr':"
!
! --- Basemap for the top frame 
	 Write(19,*) "psbasemap -X4 -Y17 ", " "//trim(adjustl(b_option)), & 	 		                    
	 	 			            " "//trim(adjustl(r_option)), " -P -JX14/8  -K >  plot.ps" 
!
! --- A title for top frame
	 open (8,file='tmpctitle1',status='unknown')
	 Write(8,*) iears_max/2, trend_max+trend_max/6, " 16 0 1 BC Scatterplot from file: ", trim(adjustl(datafile)) ; close(8)  
	 Write(19,*) "pstext -N tmpctitle1", " -JX -R -O -K >> plot.ps"
!
! --- Filters the useful data	 
         Write(19,*) "awk '{print $1, $2, $3}' tgauges-scpl.dat > tgauges-scpl.tmp"
!	 
! --- Top frame: data with errorbars  
	 Write(19,*) "psxy tgauges-scpl.tmp ", "-Ey0.25/2 -B -R -JX -O -K >> plot.ps" 
!
! --- Basemap for the bottom frame 
	 Write(19,*) "psbasemap -U'SELEN 3.2' -X0 -Y-12 ", " "//trim(adjustl(b_option)), & 
	 		                        " "//trim(adjustl(r_option)), " -P -JX14/8  -O -K >>  plot.ps" 	
! --- A title for bottom frame
	 open (8,file='tmpctitle2',status='unknown')
	 Write(8,*) iears_max/2, trend_max+trend_max/6, " 16 0 1 BC Scatterplot from file: ", trim(adjustl(datafile)) ; close(8)  
	 Write(19,*) "pstext -N tmpctitle2", " -JX -R -O -K >> plot.ps"
!	 
! --- Bottom frame: data w/o errorbars, only dots   
	 Write(19,*) "psxy tgauges-scpl.tmp ", "-Sc0.1 -G0 -B -R -JX -O -K >> plot.ps" 
!	 
! --- A note about the average rate, in the bottom frame
	open (8,file='tmpctitle3',status='unknown')
	Write(8,*) "Average (weighted) rate:"
	Write(8,'(f5.2,a1,a3,a1,f5.2,a6)') avew, " ", plus_minus, " ", davew, ' mm/yr'
	close(8) 
	open (8,file='tmpctitle3',status='unknown')
	read(8,'(a80)') string1 
	read(8,'(a80)') string2		
	close(8) 
	open (8,file='tmpctitle4',status='unknown')
	Write(8,*) iears_max-iears_max/4, trend_max-1*delta/1,   " 12 0 1 BC ", trim(adjustl(string1)) ; close(8)
	Write(19,*) "pstext -N tmpctitle4", " -JX -R -O -K >> plot.ps" 
	open (8,file='tmpctitle5',status='unknown')
	Write(8,*) iears_max-iears_max/4, trend_max-3*delta/2, " 14 0 2 BC ", trim(adjustl(string2)) ; close(8)	
	Write(19,*) "pstext -N tmpctitle5", " -G0 -JX -R -O >> plot.ps" 	
	Write(19,*) "mv plot.ps tgauges-scpl.ps"
!
! --- Closing 
	close(19) 
!
 END SUBROUTINE MAKE_TGAUGESSCA
!
!
!
!
!
!
 SUBROUTINE MAKE_TGAUGES (DATAFILE, FORMAT_OPTION, FILE_GMT)
 IMPLICIT NONE
 CHARACTER*30 DATAFILE
 CHARACTER*20 FILE_GMT
 CHARACTER*1  FORMAT_OPTION
 INTEGER I, J, K, N, NH, NALL, IEARS, NGE30, NGE60, NGE90 
!
! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +   
! --- This routine creates a GMT script for plotting the spatial distribution 
!     of the tide gauges in database DATAFILE. Last modified GS Nov. 07 2007 
! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +   
! Retouched on April 2010 by GS
! Retouched on May 2011 by GS and GG 
!
 INTEGER, PARAMETER :: MAXN=100000 ! A large number  
 CHARACTER*100 RIGA		   ! A row 
 CHARACTER*8 PCODE                 ! Tide gauge code  
 CHARACTER*3 GLOSS_CODE            ! GLOSS code 
 CHARACTER*3 YEARS                 ! Number of years 
 CHARACTER*11 RANGE_YEARS	   ! Range of years 
 CHARACTER*9 DATUM		   ! Trend 
 CHARACTER*3 PLUS_MINUS 	   ! +/- 
 CHARACTER*5 ERROR		   ! Error 
 CHARACTER*7 STDV		   ! Std. deviation of residues 						       
 CHARACTER*2 LAT1, LAT2 	   ! Latitide degrees and minutes      
 CHARACTER*3 LON1		   ! Longitude degrees 
 CHARACTER*2 LON2		   ! Longitude minutes 
 CHARACTER*30 NAME		   ! Name of the station 
 CHARACTER*1 WLAT, WLON 	   ! Where is the station (e, w, s, n)  
 REAL*8 LAT, NLAT, DLAT 	   ! Latitude  
 REAL*8 LON, NLON, DLON 	   ! Longitude  		  
!
!
!
! ====== Part 1: Analysis of tide gauges database...
!
! --- Target files for lon-lat of tide gauge sites: 
!     	Unit 39: all stations, 
!     	Unit 40: stations with >= 30 years of data, 
!     	Unit 41: stations with >= 60 years of data, 
!     	Unit 42: stations with >= 90 years of data...
! 
       open(39,file='lon-lat-tgauges-all.dat', status='unknown')
       open(40,file='lon-lat-tgauges-ge30.dat',status='unknown')
       open(41,file='lon-lat-tgauges-ge60.dat',status='unknown')
       open(42,file='lon-lat-tgauges-ge90.dat',status='unknown')
!
! --- An header 
	Write(39,*) " lon-lat (deg) & name of tide gauges sites - all data from <<rlr-trends.txt>> "
	Write(40,*) " lon-lat (deg) & name of tide gauges sites - only stations with >= 30 yrs of data"
	Write(41,*) " lon-lat (deg) & name of tide gauges sites - only stations with >= 60 yrs of data"
	Write(42,*) " lon-lat (deg) & name of tide gauges sites - only stations with >= 90 yrs of data"	
!
! --- Open and reads the tide-gauges data file for getting info...
!
! --- Counting the header lines (beginning with '#'), and the data lines
       open(10,file=DATAFILE,status='unknown')
       nh=0 ; n=0
       do i=1, maxn 
       		read(10,'(a100)',end=1) riga
       		if(riga(1:1)=='#') nh=nh+1
       		if(riga(1:1)/='#') n=n+1
       enddo 
 1     close(10)  
!
! --- Counters
	nall=n ; nge30=0 ; nge60=0 ; nge90=0
!
! --- Opening the tide gauge data file... 
       open(10,file=DATAFILE,status='unknown')
!
! --- Reading again the header lines  
       do i=1, nh ; read(10,'(a100)',end=1) riga ; enddo  
!
! *************************************************************************************************
! *************************************************************************************************
!
       If(FORMAT_OPTION=='0')       THEN  
!
!
! --- Loop on the tide gauges stations 
       do 5 i=1, n
! 
! --- Reading one line   
       	     read  (10,200) pcode, gloss_code, years, range_years, & 
   		     	    datum , plus_minus, error, stdv, & 
			    lat1, lat2, wlat, lon1, lon2, wlon, name      
!
! --- Converts the range in floating-point format 
	     open(50,file='junk-tgauges.tmp',status='unknown') ; write(50,'(a3)') years ; close(50) 
	     open(50,file='junk-tgauges.tmp',status='unknown') ; read (50,*)      iears ; close(50)	     
!
! --- Extracts lon and lat of the station from the texst strings... 
             call  find_lon_lat (lat1, lat2, wlat, lon1, lon2, wlon, lon, lat) 
!
! --- Reports lon-lat on lon-lat files... 
             write(39,*) lon, lat, name 
	     if(iears>=30)  then 
	     			nge30 = nge30 + 1 
				write(40,*) lon, lat, name 
				endif			    
	     if(iears>=60)  then
	     			nge60 = nge60 + 1 
				write(41,*) lon, lat, name 
				endif 		    
	     if(iears>=90)  then 					
	     			nge90 = nge90 + 1 
				write(42,*) lon, lat, name 
	     endif				     			
!
! --- End of Loop
  5     continue 	
!
!
	ELSEIF(FORMAT_OPTION=='1')   THEN
!
        Do 6 i=1, n 
!
! --- Reading one line      		     	    
        Read(10,119)  j, j, lon, lat, j, j, iears 
!
! --- End of Loop
  6     continue 	  
!
!
        ENDIF
!
! *************************************************************************************************
! *************************************************************************************************
!
      close(10) ; do i=39, 42 ; close(i) ; enddo
!
! --- Reading formats 
200    format(a8,1x,a3,1x,a3,2x,a11,1x,a8,1x, & 
              a3,1x,a5,1x,a7,3x,a2,1x,a2,1x,a1,1x,& 
              a3,1x,a2,1x,a1,3x,a30)
!
119  FORMAT (2(i5, 1x), 2(f10.4, 1x), 3(i5, 1x), F7.2, 1x, F7.2, 4x, A30)
!
!
! ====== Part 2: Creates a GMT script for plotting a map of tide gauges sites  
!
! ------ Target GMT file for a map of RSL sites 
!
 	open(29,file=file_gmt,status='unknown')
!
	Write(29,*) "gmtset PAPER_MEDIA A4+" 
	Write(29,*) "gmtset HEADER_FONT_SIZE 24p"
	Write(29,*) "gmtset FRAME_WIDTH 0.1c"
	Write(29,*) "gmtset ANOT_FONT_SIZE 12p"
!
! ------ Four frames (Mercator projection)
!
! --- Frame with all data 
	Write(29,*) "psbasemap -X5 -Y11 -Ba180/a80f80Wsen -R0/360/-80/80  -JM9 -K > map-tgauges.ps"
 	Write(29,*) "pscoast -G0/120/0 -S0/0/220 -B -R -O -K  -JM -Dc  -A10000 >> map-tgauges.ps"
 	Write(29,*) "psxy -H1 lon-lat-tgauges-all.dat -B -R -JM -Sc0.08 -G220 -O -K >> map-tgauges.ps"
	open(18,file='titleall.tmp',status='unknown') 
	Write(18,*) " 180 84.4 16 0 1 BC all tide gauge stations"
	Write(18,*) " 180 82.2 12 0 2 BC N=", nall ; close(18) 	
	Write(29,*) "pstext -N titleall.tmp -G0 -JM -R -O -K>> map-tgauges.ps"
	Write(29,*) ""
!
! --- Title for all frames 
	open(18,file='titletgauges.tmp',status='unknown') 
	Write(18,*) "0 87.4 24 0 3 BL Distribution of tide gauge stations as of 1/22/07" ; close(18) 		
	Write(29,*) "pstext -N titletgauges.tmp -G0 -JM -R -O -K >> map-tgauges.ps"
	Write(29,*) ""
!
! --- Frame with data with age >= 30 yrs  
	Write(29,*) "psbasemap -X10 -Y0 -Ba180/a80f80wsEn -R -JM -O -K >> map-tgauges.ps"
 	Write(29,*) "pscoast -G120 -S0/0/220 -B -R -O -K  -JM -Dc  -A10000 >> map-tgauges.ps"
 	Write(29,*) "psxy -H1 lon-lat-tgauges-ge30.dat -B -R -JM -Sc0.08 -G220 -O -K >> map-tgauges.ps"
	open(18,file='titlege30.tmp',status='unknown') 
 	Write(18,*) "180 84.4 16 0 1 BC more than 30 years of data"
 	Write(18,*) "180 82.2 12 0 2 BC N=", nge30  ; close(18) 
 	Write(29,*) "pstext -N titlege30.tmp -G0 -JM -R -O -K >> map-tgauges.ps"
	Write(29,*) ""
!
! --- Frame with data with age >= 60 yrs  
 	Write(29,*) "psbasemap -X-10 -Y-9 -Ba180/a80f80WSen -R -JM -O -K -U"//"'SELEN 3.2'", " >> map-tgauges.ps"
 	Write(29,*) "pscoast -G120 -S0/0/220 -B -R -O -K  -JM -Dc  -A10000 >> map-tgauges.ps"
 	Write(29,*) "psxy -H1 lon-lat-tgauges-ge60.dat -B -R -JM -Sc0.08 -G220 -O -K >> map-tgauges.ps"
	open(18,file='titlege60.tmp',status='unknown') 
 	Write(18,*) "180 84.4 16 0 1 BC more than 60 years of data"
 	Write(18,*) "180 82.2 12 0 2 BC N=", nge60 ; close(18)  
 	Write(29,*) "pstext -N titlege60.tmp -G0 -JM -R  -O -K >> map-tgauges.ps"
	Write(29,*) ""
!
! --- Frame with data with age >= 90 yrs 
	Write(29,*) "psbasemap -X10 -Y0 -Ba180/a80f80wSEn -R -JM -O -K >> map-tgauges.ps"
 	Write(29,*) "pscoast -G120 -S0/0/220 -B -R -O -K  -JM -Dc  -A10000 >> map-tgauges.ps"
 	Write(29,*) "psxy -H1 lon-lat-tgauges-ge90.dat -B -R -JM -Sc0.08 -G220 -O -K >> map-tgauges.ps"
	open(18,file='titlege90.tmp',status='unknown')  	
 	Write(18,*) "180 84.4 16 0 1 BC more than 90 years of data"
	Write(18,*) "180 82.2 12 0 2 BC N=", nge90 ; close(18)  !
	Write(29,*) "pstext -N titlege90.tmp -G0 -JM -R -O >> map-tgauges.ps"
!
	close(29) 	
!	
 END SUBROUTINE MAKE_TGAUGES
!
!
!
!
!
!
 SUBROUTINE MAKE_RMAPS (TITLICE, &
 			RESOLUTION, & 
			NV, & 
			CODE, & 
			ITER, & 
			MODE, & 
			DEGREE, & 
			VSTRING, & 
			OPT, & 
			SHORT_VISCO)
!			
! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + 
! Prepares a GMT script for *** REGIONAL maps *** of dot-S, U & N at present time 
! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + 
! --- Last change GS March 30 08 --- 
! Also revised April 2010 for ALMA coupling 
!
!
 IMPLICIT NONE
 INTEGER, PARAMETER :: NRG=10 
 CHARACTER*1  OPT(0:NRG)
 CHARACTER*80 TITRE
 CHARACTER*30 NAMEIN, NAMEOUT, NAMEF   
 CHARACTER*10 TITLICE, RESOLUTION, DEGREE
 CHARACTER*30 VSTRING
 CHARACTER*20 FILE1_GMT
 CHARACTER*20 R_OPTION, T_OPTION
 CHARACTER*3  NV, CODE
 CHARACTER*1  ITER, MODE 
 CHARACTER*33 TABLE_NAME, TABLE_NAME_SU, TABLE_NAME_N, SCRIPT_NAME  
 CHARACTER*33 TMP_FILE, DATA_FILE, PS_FILE 
 CHARACTER*99 TMP_TITLE
 CHARACTER*100 SHORT_VISCO
 INTEGER K, U, US 
!
!
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This subroutine prepares a GMT script for plotting dot S, N & U 
! across specific regions. Regions available (as of March 2008): 
!  
! -1) Italy  
! -2) Mediterranean 
! -3) Europe 
! -4) Fennoscandia 
! -5) Greenland 
! -6) North America  
! -7) Antarctica 
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
!
!
  If(opt(1)=='y') then 
!
! ################      ################     ################     ################
! Region 1 - ITALY      Region 1 - ITALY     Region 1 - ITALY	  Region 1 - ITALY
! ################      ################     ################	  ################
!
  table_name = "pale-italy.cpt" ; script_name = "italy.gmt"
!   
  Call color_tables (table_name)  
!
  us=11 
  Open (us,file=script_name,status='unknown') 
  Write(us,*)"gmtset PAPER_MEDIA A4+"
  Write(us,*)"gmtset LABEL_FONT_SIZE 24p"
  Write(us,*)"gmtset ANOT_FONT_SIZE 16p" 
  Write(us,*)"gmtset FRAME_WIDTH 0.1c"
  Write(us,*)" " 
!
  do 100 u=12, 14 
! 
  if    (u==12) then 
  	    tmp_file = "sdot-italy.tmp" ; tmp_title= "Rate of sea level change today - ITALY"
	    data_file= "sdotmap.dat"    ; ps_file  = "sdot-italy.ps"  
  elseif(u==13) then
  	    tmp_file = "udot-italy.tmp" ; tmp_title= "Rate of vertical uplift today - ITALY"
	    data_file= "udotmap.dat"    ; ps_file  = "udot-italy.ps"
  elseif(u==14) then
  	    tmp_file = "ndot-italy.tmp" ; tmp_title= "Rate of geoid height variation today - ITALY"
	    data_file= "ndotmap.dat"    ; ps_file  = "ndot-italy.ps"
  Endif 
!
!
! ..............................................................................
  Open (u,file=tmp_file,status='unknown') 
!
  Write(u,*)"13 +49.0 24 0 3 BC ", trim(adjustl(tmp_title))
!
  If(CODE=='-1') then 
!
  Write(u,*)"13 +48.0 16 0 0 BC ", & 
  		               " -Ice model: ", trim(adjustl(titlice)), & 
  			       " -ALMA rheology: ", trim(adjustl(short_visco))
  Write(u,*)"13 +30.0 16 0 0 BC ", & 
  	                       " -LMAX=", trim(adjustl(DEGREE)), & 
			       " -RES=",  trim(adjustl(RESOLUTION)), & 
			       " -NV=",   trim(adjustl(NV)), & 
			       " -MODE=", trim(adjustl(MODE)), & 
			       " -ITER=", trim(adjustl(ITER))	
  Else
!
  Write(u,*)"13 +48.0 16 0 0 BC ", & 
  		               " -Ice model: ", trim(adjustl(titlice)), & 
  			       " -Viscosity profile: ", trim(adjustl(vstring))
  Write(u,*)"13 +30.0 16 0 0 BC ", & 
  	                       " -LMAX=", trim(adjustl(DEGREE)), & 
			       " -RES=",  trim(adjustl(RESOLUTION)), & 
			       " -NV=",   trim(adjustl(NV)), & 
			       " -CODE=", trim(adjustl(CODE)), & 
			       " -MODE=", trim(adjustl(MODE)), & 
			       " -ITER=", trim(adjustl(ITER))	
  Endif
!
  Close(u)			       		       
! ..............................................................................
!
!
  Write(us,*)"psbasemap -P -X4 -Y10 -Ba2f2WSEn -R5/21/35/47 -JJ13/14 -K > ", & 
  	      trim(adjustl(ps_file))
  Write(us,*)"pscontour -I -JJ -R -O -K ", trim(adjustl(data_file)), " -C"//trim(adjustl(table_name)), " >> ", & 
              trim(adjustl(ps_file)) 
  Write(us,*) "pscontour -G16 -W1/255 -A -JJ -R -O -K ", trim(adjustl(data_file)), & 
              " -C"//trim(adjustl(table_name)), " >> ",  trim(adjustl(ps_file))	    
  Write(us,*)"pscoast -R -JJ -Di -B -W2/0 -A1000 -O -K >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)"pstext -N -R -JJ -B -G0 -O -K ", trim(adjustl(tmp_file)), " >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)"psscale -E -U/0.5/0.5/'SELEN 3.2' ", "-C"//trim(adjustl(table_name)), " -Bf0.1a0.5/:mm/yr: -D7/-2/10/1h -O >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)" "
!
  100 Continue 
!
  close(us)   
!
  ENDIF
!
! ----------------      ----------------     ----------------     ----------------
! Region 1 - ITALY      Region 1 - ITALY     Region 1 - ITALY	  Region 1 - ITALY
! ----------------      ----------------     ----------------	  ----------------
!
  If(opt(2)=='y') then 
!
! ########################     ########################    ########################
! Region 2 - MEDITERRANEAN     Region 2 - MEDITERRANEAN    Region 2 - MEDITERRANEAN
! ########################     ########################    ########################
!
  table_name = "pale-mediterranean.cpt" ; script_name = "mediterranean.gmt"
!   
  Call color_tables (table_name)  
!
  us=11 
  Open (us,file=script_name,status='unknown') 
  Write(us,*)"gmtset PAPER_MEDIA A4+"
  Write(us,*)"gmtset LABEL_FONT_SIZE 24p"
  Write(us,*)"gmtset ANOT_FONT_SIZE 16p" 
  Write(us,*)"gmtset FRAME_WIDTH 0.1c"
  Write(us,*)" " 
!
  do 200 u=12, 14 
! 
  if    (u==12) then 
  	    tmp_file = "sdot-mediterranean.tmp" ; tmp_title= "Rate of sea level change today - MEDITERRANEAN"
	    data_file= "sdotmap.dat"    ; ps_file  = "sdot-mediterranean.ps"  
  elseif(u==13) then
  	    tmp_file = "udot-mediterranean.tmp" ; tmp_title= "Rate of vertical uplift today - MEDITERRANEAN"
	    data_file= "udotmap.dat"    ; ps_file  = "udot-mediterranean.ps"
  elseif(u==14) then
  	    tmp_file = "ndot-mediterranean.tmp" ; tmp_title= "Rate of geoid height variation today - MEDITERRANEAN"
	    data_file= "ndotmap.dat"    ; ps_file  = "ndot-mediterranean.ps"
  Endif 
!
! ..............................................................................
  Open (u,file=tmp_file,status='unknown') 
!
  Write(u,*)"20 +53.0   22 0 3 BC ", trim(adjustl(tmp_title))
!
  If(CODE=='-1') then 
!
  Write(u,*)"20 +51.5 16 0 0 BC", & 
  		               " -Ice model: ", trim(adjustl(titlice)), & 
  			       " -ALMA rheology: ", trim(adjustl(short_visco))
  Write(u,*)"20 +21.0 16 0 0 BC", & 
  	                       " -LMAX=", trim(adjustl(DEGREE)), & 
			       " -RES=",  trim(adjustl(RESOLUTION)), & 
			       " -NV=",   trim(adjustl(NV)), & 
			       " -MODE=", trim(adjustl(MODE)), & 
			       " -ITER=", trim(adjustl(ITER))	
  Else
!
  Write(u,*)"20 +51.5 16 0 0 BC", & 
  		               " -Ice model: ", trim(adjustl(titlice)), & 
  			       " -Viscosity profile: ", trim(adjustl(vstring))
  Write(u,*)"20 +21.0 16 0 0 BC", & 
  	                       " -LMAX=", trim(adjustl(DEGREE)), & 
			       " -RES=",  trim(adjustl(RESOLUTION)), & 
			       " -NV=",   trim(adjustl(NV)), & 
			       " -CODE=", trim(adjustl(CODE)), & 
			       " -MODE=", trim(adjustl(MODE)), & 
			       " -ITER=", trim(adjustl(ITER))	
  Endif
!
  Close(u)			       		       
! ..............................................................................
!
!  
  Write(us,*)"psbasemap -X3 -Y6 -Ba8f4WSEn -R-8/48/30/50 -JJ18/24 -K > ", & 
  	      trim(adjustl(ps_file))
  Write(us,*)"pscontour -I -JJ -R -O -K ", trim(adjustl(data_file)), " -C"//trim(adjustl(table_name)), " >> ", & 
              trim(adjustl(ps_file)) 
  Write(us,*)"pscontour -G16 -W1/255 -A -JJ -R -O -K ", trim(adjustl(data_file)), " -C"//trim(adjustl(table_name)), " >> ", & 
              trim(adjustl(ps_file))	     	      
  Write(us,*)"pscoast -R -JJ -Di -B -W2/0 -A1000 -O -K >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)"pstext -N -R -JJ -B -G0 -O -K ", trim(adjustl(tmp_file)), " >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)"psscale -E -U/0.5/0.5/'SELEN 3.2' ", "-C"//trim(adjustl(table_name)), " -Bf0.1a0.5/:mm/yr: -D12/-1.5/10/1h -O >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)" "
!
  200 Continue 
!
  close(us)  
!
  ENDIF  
!
! ------------------------     ------------------------    ------------------------
! Region 2 - MEDITERRANEAN     Region 2 - MEDITERRANEAN    Region 2 - MEDITERRANEAN
! ------------------------     ------------------------    ------------------------
!
  If(opt(3)=='y') then 
!
! #################    #################    #################    #################
! Region 3 - EUROPE    Region 3 - EUROPE    Region 3 - EUROPE	 Region 3 - EUROPE
! #################    #################    #################	 #################
!
  table_name = "pale-europe.cpt" ; script_name = "europe.gmt"
!   
  Call color_tables (table_name)  
!
  us=11 
  Open (us,file=script_name,status='unknown') 
  Write(us,*)"gmtset PAPER_MEDIA A4+"
  Write(us,*)"gmtset LABEL_FONT_SIZE 24p"
  Write(us,*)"gmtset ANOT_FONT_SIZE 16p" 
  Write(us,*)"gmtset FRAME_WIDTH 0.1c"
  Write(us,*)" " 
!
  do 300 u=12, 14 
! 
  if    (u==12) then 
  	    tmp_file = "sdot-europe.tmp" ; tmp_title= "Rate of sea level change today - EUROPE"
	    data_file= "sdotmap.dat"    ; ps_file  = "sdot-europe.ps"  
  elseif(u==13) then
  	    tmp_file = "udot-europe.tmp" ; tmp_title= "Rate of vertical uplift today - EUROPE"
	    data_file= "udotmap.dat"    ; ps_file  = "udot-europe.ps"
  elseif(u==14) then
  	    tmp_file = "ndot-europe.tmp" ; tmp_title= "Rate of geoid height variation today - EUROPE"
	    data_file= "ndotmap.dat"    ; ps_file  = "ndot-europe.ps"
  Endif 
!
!
! ..............................................................................
  Open (u,file=tmp_file,status='unknown') 
!
  Write(u,*)"20 +80.0 26 0 3 BC ", trim(adjustl(tmp_title))
!
  If(CODE=='-1') then 
!
  Write(u,*)"20 +77.0 16 0 0 BC ", & 
  		               " -Ice model: ", trim(adjustl(titlice)), & 
  			       " -ALMA rheology: ", trim(adjustl(short_visco))
  Write(u,*)"20 -2.0 16 0 0 BC ", & 
  	                       " -LMAX=", trim(adjustl(DEGREE)), & 
			       " -RES=",  trim(adjustl(RESOLUTION)), & 
			       " -NV=",   trim(adjustl(NV)), & 
			       " -MODE=", trim(adjustl(MODE)), & 
			       " -ITER=", trim(adjustl(ITER))	
  Else
!
  Write(u,*)"20 +77.0 16 0 0 BC ", & 
  		               " -Ice model: ", trim(adjustl(titlice)), & 
  			       " -Viscosity profile: ", trim(adjustl(vstring))
  Write(u,*)"20 -2.0 16 0 0 BC ", & 
  	                       " -LMAX=", trim(adjustl(DEGREE)), & 
			       " -RES=",  trim(adjustl(RESOLUTION)), & 
			       " -NV=",   trim(adjustl(NV)), & 
			       " -CODE=", trim(adjustl(CODE)), & 
			       " -MODE=", trim(adjustl(MODE)), & 
			       " -ITER=", trim(adjustl(ITER))	
  Endif
!
  Close(u)			       		       
! ..............................................................................
!
!
!  
  Write(us,*)"psbasemap -X3 -Y5 -Ba20f10/a10f5WSEn -R-40/80/25/75 -JJ20/18 -K > ", & 
  	      trim(adjustl(ps_file))
  Write(us,*)"pscontour -I -JJ -R -O -K ", trim(adjustl(data_file)), " -C"//trim(adjustl(table_name)), " >> ", & 
              trim(adjustl(ps_file)) 
  Write(us,*)"pscontour -G16 -W1/255 -A -JJ -R -O -K ", trim(adjustl(data_file)), " -C"//trim(adjustl(table_name)), " >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)"pscoast -R -JJ -Di -B -W2/0 -A1000 -O -K >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)"pstext -N -R -JJ -B -G0 -O -K ", trim(adjustl(tmp_file)), " >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)"psscale -E -U/0.5/0.5/'SELEN 3.2' ", "-C"//trim(adjustl(table_name)), " -Bf1a4/:mm/yr: -D9/-1.5/10/1h -O >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)" "
!
  300 Continue 
!
  close(us)  
!
  ENDIF
!
! -----------------    -----------------    -----------------    -----------------
! Region 3 - EUROPE    Region 3 - EUROPE    Region 3 - EUROPE	 Region 3 - EUROPE
! -----------------    -----------------    -----------------	 -----------------
!
  If(opt(4)=='y') then 
!
! #######################  #######################  #######################
! Region 4 - FENNOSCANDIA  Region 4 - FENNOSCANDIA  Region 4 - FENNOSCANDIA
! #######################  #######################  #######################
!
  table_name_su = "pale-fennoscandia-su.cpt" 
  Call color_tables (table_name_su)  
!
  table_name_n = "pale-fennoscandia-n.cpt" 
  Call color_tables (table_name_n) 
!
  script_name = "fennoscandia.gmt"
!
  us=11 
  Open (us,file=script_name,status='unknown') 
  Write(us,*)"gmtset PAPER_MEDIA A4+"
  Write(us,*)"gmtset LABEL_FONT_SIZE 24p"
  Write(us,*)"gmtset ANOT_FONT_SIZE 16p" 
  Write(us,*)"gmtset FRAME_WIDTH 0.1c"
  Write(us,*)" " 
!
  do 400 u=12, 14 
! 
  if    (u==12) then 
  	    tmp_file = "sdot-fennoscandia.tmp" ; tmp_title= "Rate of sea level change today - FENNOSCANDIA"
	    data_file= "sdotmap.dat"    ; ps_file  = "sdot-fennoscandia.ps"  
  elseif(u==13) then
  	    tmp_file = "udot-fennoscandia.tmp" ; tmp_title= "Rate of vertical uplift today - FENNOSCANDIA"
	    data_file= "udotmap.dat"    ; ps_file  = "udot-fennoscandia.ps"
  elseif(u==14) then
  	    tmp_file = "ndot-fennoscandia.tmp" ; tmp_title= "Rate of geoid height variation today - FENNOSCANDIA"
	    data_file= "ndotmap.dat"    ; ps_file  = "ndot-fennoscandia.ps"
  Endif 
!
!
! ..............................................................................
  Open (u,file=tmp_file,status='unknown') 
!
  Write(u,*)"30 +77.0 24 0 3 BC ", trim(adjustl(tmp_title))
!
  If(CODE=='-1') then 
!
  Write(u,*)"30 +76.0 16 0 0 BC ", & 
  		               " -Ice model: ", trim(adjustl(titlice)), & 
  			       " -ALMA rheology: ", trim(adjustl(short_visco))
  Write(u,*)"30 +40.0 16 0 0 BC ", & 
  	                       " -LMAX=", trim(adjustl(DEGREE)), & 
			       " -RES=",  trim(adjustl(RESOLUTION)), & 
			       " -NV=",   trim(adjustl(NV)), & 
			       " -MODE=", trim(adjustl(MODE)), & 
			       " -ITER=", trim(adjustl(ITER))	
  Else
!
  Write(u,*)"30 +76.0 16 0 0 BC ", & 
  		               " -Ice model: ", trim(adjustl(titlice)), & 
  			       " -Viscosity profile: ", trim(adjustl(vstring))
  Write(u,*)"30 +40.0 16 0 0 BC ", & 
  	                       " -LMAX=", trim(adjustl(DEGREE)), & 
			       " -RES=",  trim(adjustl(RESOLUTION)), & 
			       " -NV=",   trim(adjustl(NV)), & 
			       " -CODE=", trim(adjustl(CODE)), & 
			       " -MODE=", trim(adjustl(MODE)), & 
			       " -ITER=", trim(adjustl(ITER))	
  Endif
!
  Close(u)			       		       
! ..............................................................................

!
!  
  if(u==12.or.u==13) table_name=table_name_su
  if(u==14) 	     table_name=table_name_n
!
  Write(us,*)"psbasemap -X3 -Y5 -Ba10f5WSEn -R0/60/50/75 -JJ26/20 -K  > ", & 
  	      trim(adjustl(ps_file))
  Write(us,*)"pscontour -I -JJ -R -O -K ", trim(adjustl(data_file)), " -C"//trim(adjustl(table_name)), " >> ", & 
              trim(adjustl(ps_file)) 
  Write(us,*)"pscontour -G16 -W1/255 -A -JJ -R -O -K ", trim(adjustl(data_file)), " -C"//trim(adjustl(table_name)), " >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)"pscoast -R -JJ -Di -B -W2/0 -A1000 -O -K >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)"pstext -N -R -JJ -B -G0 -O -K ", trim(adjustl(tmp_file)), " >> ", & 
              trim(adjustl(ps_file))
  if(u==12.or.u==13) &   
  Write(us,*)"psscale -E -U/0.5/0.5/'SELEN 3.2' ", "-C"//trim(adjustl(table_name)), " -Bf1a5/:mm/yr: -D10/-1.5/10/1h -O >> ", & 
              trim(adjustl(ps_file))
  if(u==14) &   
  Write(us,*)"psscale -E -U/0.5/0.5/'SELEN 3.2' ", "-C"//trim(adjustl(table_name)), " -Bf2a1/:mm/yr: -D10/-1.5/10/1h -O >> ", & 
              trim(adjustl(ps_file))	      	      	      	      
  Write(us,*)" "
!
  400 Continue 
!
  close(us) 
!
  ENDIF
!
! -----------------------   -----------------------   -----------------------
! Region 4 - Fennoscandia   Region 4 - Fennoscandia   Region 4 - Fennoscandia
! -----------------------   -----------------------   -----------------------
!
  If(opt(5)=='y') then 
!
! ####################       ####################       ####################
! Region 5 - GREENLAND       Region 5 - GREENLAND       Region 5 - GREENLAND
! ####################       ####################       ####################
!
  table_name_su = "pale-greenland-su.cpt" 
  Call color_tables (table_name_su)  
!
  table_name_n = "pale-greenland-n.cpt" 
  Call color_tables (table_name_n) 
!
  script_name = "greenland.gmt"
!
  us=11 
  Open (us,file=script_name,status='unknown') 
  Write(us,*)"gmtset PAPER_MEDIA A4+"
  Write(us,*)"gmtset LABEL_FONT_SIZE 24p"
  Write(us,*)"gmtset ANOT_FONT_SIZE 16p" 
  Write(us,*)"gmtset FRAME_WIDTH 0.1c"
  Write(us,*)" " 
!
  do 500 u=12, 14 
! 
  if    (u==12) then 
  	    tmp_file = "sdot-greenland.tmp" ; tmp_title= "Rate of sea level change today - GREENLAND"
	    data_file= "sdotmap.dat"    ; ps_file  = "sdot-greenland.ps"  
  elseif(u==13) then
  	    tmp_file = "udot-greenland.tmp" ; tmp_title= "Rate of vertical uplift today - GREENLAND"
	    data_file= "udotmap.dat"    ; ps_file  = "udot-greenland.ps"
  elseif(u==14) then
  	    tmp_file = "ndot-greenland.tmp" ; tmp_title= "Rate of geoid height variation today - GREENLAND"
	    data_file= "ndotmap.dat"    ; ps_file  = "ndot-greenland.ps"
  Endif 
!
! ..............................................................................
  Open (u,file=tmp_file,status='unknown') 
!
  Write(u,*)"320  91   22 0 3 BC ", trim(adjustl(tmp_title))
!
  If(CODE=='-1') then 
!
  Write(u,*)"320  89.5 16 0 0 BC ", & 
  		               " -Ice model: ", trim(adjustl(titlice)), & 
  			       " -ALMA rheology: ", trim(adjustl(short_visco))
  Write(u,*)"320  30   16 0 0 BC ", & 
  	                       " -LMAX=", trim(adjustl(DEGREE)), & 
			       " -RES=",  trim(adjustl(RESOLUTION)), & 
			       " -NV=",   trim(adjustl(NV)), & 
			       " -MODE=", trim(adjustl(MODE)), & 
			       " -ITER=", trim(adjustl(ITER))	
  Else
!
  Write(u,*)"320  89.5 16 0 0 BC ", & 
  		               " -Ice model: ", trim(adjustl(titlice)), & 
  			       " -Viscosity profile: ", trim(adjustl(vstring))
  Write(u,*)"320  30   16 0 0 BC ", & 
  	                       " -LMAX=", trim(adjustl(DEGREE)), & 
			       " -RES=",  trim(adjustl(RESOLUTION)), & 
			       " -NV=",   trim(adjustl(NV)), & 
			       " -CODE=", trim(adjustl(CODE)), & 
			       " -MODE=", trim(adjustl(MODE)), & 
			       " -ITER=", trim(adjustl(ITER))	
  Endif
!
  Close(u)			       		       
! ..............................................................................
! 
! 
  if(u==12.or.u==13) table_name=table_name_su
  if(u==14) 	     table_name=table_name_n
!
  Write(us,*)"psbasemap -P -X3 -Y9 -Ba10f10/f10a10WSEn -R-80/0/50/87.5 -JJ-40/16 -K  > ", & 
  	      trim(adjustl(ps_file))
  Write(us,*)"pscontour -I -JJ -R -O -K ", trim(adjustl(data_file)), " -C"//trim(adjustl(table_name)), " >> ", & 
              trim(adjustl(ps_file)) 
  Write(us,*)"pscontour -G16 -W1/255 -A -JJ -R -O -K ", trim(adjustl(data_file)), " -C"//trim(adjustl(table_name)), " >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)"pscoast -R -JJ -Di -B -W2/0 -A1000 -O -K >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)"pstext -N -R -JJ -B -G0 -O -K ", trim(adjustl(tmp_file)), " >> ", & 
              trim(adjustl(ps_file))
  if(u==12.or.u==13) &   
  Write(us,*)"psscale -U/0.5/13/'SELEN 3.2' -E ", "-C"//trim(adjustl(table_name)), " -Bf1a5/:mm/yr: -D8/-1.5/10/1h -O >> ", & 
              trim(adjustl(ps_file))
  if(u==14) &   
  Write(us,*)"psscale -U/0.5/13/'SELEN 3.2' -E ", "-C"//trim(adjustl(table_name)), " -Bf.1a.5/:mm/yr: -D8/-1.5/10/1h -O >> ", & 
              trim(adjustl(ps_file))	      	      	      	      
  Write(us,*)" "
!
  500 Continue 
!
  close(us) 
!
  ENDIF
!
! --------------------       --------------------       --------------------
! Region 5 - GREENLAND       Region 5 - GREENLAND       Region 5 - GREENLAND
! --------------------       --------------------       --------------------
!
  If(opt(6)=='y') then   
!
! ########################    ########################    ######################## 
! Region 6 - North AMERICA    Region 6 - North AMERICA    Region 6 - North AMERICA 
! ########################    ########################    ######################## 
!
  table_name_su = "pale-namerica-su.cpt" 
  Call color_tables (table_name_su)  
!
  table_name_n = "pale-namerica-n.cpt" 
  Call color_tables (table_name_n) 
!
  script_name = "north-america.gmt"
!
  us=11 
  Open (us,file=script_name,status='unknown') 
  Write(us,*)"gmtset PAPER_MEDIA A4+"
  Write(us,*)"gmtset LABEL_FONT_SIZE 24p"
  Write(us,*)"gmtset ANOT_FONT_SIZE 16p" 
  Write(us,*)"gmtset FRAME_WIDTH 0.1c"
  Write(us,*)" " 
!
  do 600 u=12, 14 
! 
  if    (u==12) then 
  	    tmp_file = "sdot-namerica.tmp" ; tmp_title= "Rate of sea level change today - N. AMERICA"
	    data_file= "sdotmap.dat"    ; ps_file  = "sdot-namerica.ps"  
  elseif(u==13) then
  	    tmp_file = "udot-namerica.tmp" ; tmp_title= "Rate of vertical uplift today - N. AMERICA"
	    data_file= "udotmap.dat"    ; ps_file  = "udot-namerica.ps"
  elseif(u==14) then
  	    tmp_file = "ndot-namerica.tmp" ; tmp_title= "Rate of geoid height variation today - N. AMERICA"
	    data_file= "ndotmap.dat"    ; ps_file  = "ndot-namerica.ps"
  Endif 
!
! ..............................................................................
  Open (u,file=tmp_file,status='unknown') 
!
  Write(u,*)"270  92 20 0 3 BC ", trim(adjustl(tmp_title))
!
  If(CODE=='-1') then 
!
  Write(u,*)"270  90 16 0 0 BC ", & 
  		               " -Ice model: ", trim(adjustl(titlice)), & 
  			       " -ALMA rheology: ", trim(adjustl(short_visco))
  Write(u,*)"270 -30 16 0 0 BC ", & 
  	                       " -LMAX=", trim(adjustl(DEGREE)), & 
			       " -RES=",  trim(adjustl(RESOLUTION)), & 
			       " -NV=",   trim(adjustl(NV)), & 
			       " -MODE=", trim(adjustl(MODE)), & 
			       " -ITER=", trim(adjustl(ITER))	
  Else
!
  Write(u,*)"270  90 16 0 0 BC ", & 
  		               " -Ice model: ", trim(adjustl(titlice)), & 
  			       " -Viscosity profile: ", trim(adjustl(vstring))
  Write(u,*)"270 -30 16 0 0 BC ", & 
  	                       " -LMAX=", trim(adjustl(DEGREE)), & 
			       " -RES=",  trim(adjustl(RESOLUTION)), & 
			       " -NV=",   trim(adjustl(NV)), & 
			       " -CODE=", trim(adjustl(CODE)), & 
			       " -MODE=", trim(adjustl(MODE)), & 
			       " -ITER=", trim(adjustl(ITER))	
  Endif
!
  Close(u)			       		       
! ..............................................................................
!
!  
  if(u==12.or.u==13) table_name=table_name_su
  if(u==14) 	     table_name=table_name_n
!
  Write(us,*)"psbasemap -P -X3 -Y9 -Ba20f10/f10a10WSEn -R200/340/10/87.5 -JJ-90/16 -K  > ", & 
  	      trim(adjustl(ps_file))
  Write(us,*)"pscontour -I -JJ -R -O -K ", trim(adjustl(data_file)), " -C"//trim(adjustl(table_name)), " >> ", & 
              trim(adjustl(ps_file)) 
  Write(us,*)"pscontour -G16 -W1/255 -A -JJ -R -O -K ", trim(adjustl(data_file)), " -C"//trim(adjustl(table_name)), " >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)"pscoast -R -JJ -Di -B -W2/0 -A1000 -O -K >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)"pstext -N -R -JJ -B -G0 -O -K ", trim(adjustl(tmp_file)), " >> ", & 
              trim(adjustl(ps_file))
  if(u==12.or.u==13) &   
  Write(us,*)"psscale -U/0.5/0.5/'SELEN 3.2' -E ", "-C"//trim(adjustl(table_name)), " -Bf1a4/:mm/yr: -D8/-1.5/10/1h -O >> ", & 
              trim(adjustl(ps_file))
  if(u==14) &   
  Write(us,*)"psscale -U/0.5/0.5/'SELEN 3.2' -E ", "-C"//trim(adjustl(table_name)), " -Bf1a2/:mm/yr: -D8/-1.5/10/1h -O >> ", & 
              trim(adjustl(ps_file))	      	      	      	      
  Write(us,*)" "
!
  600 Continue 
!
  close(us) 
!
  ENDIF
!
! ------------------------    ------------------------    ------------------------ 
! Region 6 - North AMERICA    Region 6 - North AMERICA    Region 6 - North AMERICA 
! ------------------------    ------------------------    ------------------------ 
!
  If(opt(7)=='y') then   
!
! #####################    #####################    #####################
! Region 7 - Antarctica    Region 7 - Antarctica    Region 7 - Antarctica
! #####################    #####################    #####################
!
  table_name_su = "pale-antarctica-su.cpt" 
  Call color_tables (table_name_su)  
!
  table_name_n = "pale-antarctica-n.cpt" 
  Call color_tables (table_name_n) 
!
  script_name = "antarctica.gmt"
!
  us=11 
  Open (us,file=script_name,status='unknown') 
  Write(us,*)"gmtset PAPER_MEDIA A4+"
  Write(us,*)"gmtset LABEL_FONT_SIZE 24p"
  Write(us,*)"gmtset ANOT_FONT_SIZE 16p" 
  Write(us,*)"gmtset FRAME_WIDTH 0.1c"
  Write(us,*)" " 
!
  do 700 u=12, 14 
! 
  if    (u==12) then 
  	    tmp_file = "sdot-antarctica.tmp" ; tmp_title= "Rate of sea level change today - ANTARCTICA"
	    data_file= "sdotmap.dat"    ; ps_file  = "sdot-antarctica.ps"  
  elseif(u==13) then
  	    tmp_file = "udot-antarctica.tmp" ; tmp_title= "Rate of vertical uplift today - ANTARCTICA"
	    data_file= "udotmap.dat"    ; ps_file  = "udot-antarctica.ps"
  elseif(u==14) then
  	    tmp_file = "ndot-antarctica.tmp" ; tmp_title= "Rate of geoid height variation today - ANTARCTICA"
	    data_file= "ndotmap.dat"    ; ps_file  = "ndot-antarctica.ps"
  Endif 
!
! ..............................................................................
  Open (u,file=tmp_file,status='unknown') 
!
  Write(u,*)"0 +18.0 22 0 3 BC ", trim(adjustl(tmp_title))
!
  If(CODE=='-1') then 
!
  Write(u,*)"0 +10.0 16 0 0 BC ", & 
  		               " -Ice model: ", trim(adjustl(titlice)), & 
  			       " -ALMA rheology: ", trim(adjustl(short_visco))
  Write(u,*)"0 -115.0 16 0 0 BC ", & 
  	                       " -LMAX=", trim(adjustl(DEGREE)), & 
			       " -RES=",  trim(adjustl(RESOLUTION)), & 
			       " -NV=",   trim(adjustl(NV)), & 
			       " -MODE=", trim(adjustl(MODE)), & 
			       " -ITER=", trim(adjustl(ITER))	
  Else
!
  Write(u,*)"0 +10.0 16 0 0 BC ", & 
  		               " -Ice model: ", trim(adjustl(titlice)), & 
  			       " -Viscosity profile: ", trim(adjustl(vstring))
  Write(u,*)"0 -115.0 16 0 0 BC ", & 
  	                       " -LMAX=", trim(adjustl(DEGREE)), & 
			       " -RES=",  trim(adjustl(RESOLUTION)), & 
			       " -NV=",   trim(adjustl(NV)), & 
			       " -CODE=", trim(adjustl(CODE)), & 
			       " -MODE=", trim(adjustl(MODE)), & 
			       " -ITER=", trim(adjustl(ITER))	
  Endif
!
  Close(u)			       		       
! ..............................................................................
!
!
  if(u==12.or.u==13) table_name=table_name_su
  if(u==14) 	     table_name=table_name_n
!
  Write(us,*)"pscontour -P -X4 -Y8 -I -JE0/-90/12 -R0/360/-90/-52 -K ", trim(adjustl(data_file)), & 
  	     " -C"//trim(adjustl(table_name)), " > ", trim(adjustl(ps_file)) 
  Write(us,*)"pscontour -G15 -W1/255 -A -JE -R -O -K ", trim(adjustl(data_file)), " -C"//trim(adjustl(table_name)), " >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)"pscoast -R -JE -Di -Ba30f30g30 -W2/255/0/0 -A0 -O -K >> ", & 
              trim(adjustl(ps_file))
  if(u==12.or.u==13) &   
  Write(us,*)"psscale -E -X-6 -Y-7 ", "-C"//trim(adjustl(table_name)), " -Bf1a2/:mm/yr: -D12/5/8/1h -O -K >> ", & 
              trim(adjustl(ps_file))
  if(u==14) &   
  Write(us,*)"psscale -E -X-6 -Y-7 ", "-C"//trim(adjustl(table_name)), " -Bf0.5a1/:mm/yr: -D12/5/8/1h -O -K >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)"pstext -X6 -N -R -JE -G0 -O -U/4/20.5/'SELEN 3.2' ", trim(adjustl(tmp_file)), " >> ", & 
              trim(adjustl(ps_file))	      	      	      	      
  Write(us,*)" "
!
  700 Continue 
!
  close(us) 
!
  ENDIF
!
! ---------------------    ---------------------    ---------------------
! Region 7 - Antarctica    Region 7 - Antarctica    Region 7 - Antarctica
! ---------------------    ---------------------    ---------------------
!
 END SUBROUTINE MAKE_RMAPS
!
!
!
!
!
!
 SUBROUTINE COLOR_TABLES (TABLE_NAME)
 IMPLICIT NONE
 CHARACTER*33 TABLE_NAME
 CHARACTER*5 FMT1, FMT2
!
! # Draws color tables for regional maps- 
!
 if(table_name=='pale-italy.cpt'          .or.& 
    table_name=='pale-mediterranean.cpt')     & 
    then 
!
 	open(9,file=table_name,status='unknown') 
	fmt1='(a59)'
	fmt2='(a27)'
	Write(9,fmt1) '#----------------------------------------------------------'
 	Write(9,fmt1) '#    cpt created by: makecpt -Cno_green -T-1.5/1.5/0.1     '
	Write(9,fmt1) '#----------------------------------------------------------'	
 	Write(9,fmt1) '-1.5    32      96      255     -1.4    32      96      255'
 	Write(9,fmt1) '-1.4    32      96      255     -1.3    32      96      255'
 	Write(9,fmt1) '-1.3    32      159     255     -1.2    32      159     255'
 	Write(9,fmt1) '-1.2    32      159     255     -1.1    32      159     255'
 	Write(9,fmt1) '-1.1    32      191     255     -1      32      191     255'
 	Write(9,fmt1) '-1      32      191     255     -0.9    32      191     255'
 	Write(9,fmt1) '-0.9    0       207     255     -0.8    0       207     255'
 	Write(9,fmt1) '-0.8    0       207     255     -0.7    0       207     255'
 	Write(9,fmt1) '-0.7    42      255     255     -0.6    42      255     255'
 	Write(9,fmt1) '-0.6    42      255     255     -0.5    42      255     255'
 	Write(9,fmt1) '-0.5    85      255     255     -0.4    85      255     255'
 	Write(9,fmt1) '-0.4    85      255     255     -0.3    85      255     255'
 	Write(9,fmt1) '-0.3    127     255     255     -0.2    127     255     255'
 	Write(9,fmt1) '-0.2    127     255     255     -0.1    127     255     255'
 	Write(9,fmt1) '-0.1    170     255     255     0       170     255     255'
 	Write(9,fmt1) '0       255     255     84      0.1     255     255     84 '
 	Write(9,fmt1) '0.1     255     255     84      0.2     255     255     84 '
 	Write(9,fmt1) '0.2     255     240     0       0.3     255     240     0  '
 	Write(9,fmt1) '0.3     255     240     0       0.4     255     240     0  '
 	Write(9,fmt1) '0.4     255     191     0       0.5     255     191     0  '
	Write(9,fmt1) '0.5     255     191     0       0.6     255     191     0  '
	Write(9,fmt1) '0.6     255     168     0       0.7     255     168     0  '
	Write(9,fmt1) '0.7     255     168     0       0.8     255     168     0  '
	Write(9,fmt1) '0.8     255     138     0       0.9     255     138     0  '
	Write(9,fmt1) '0.9     255     138     0       1       255     138     0  '
	Write(9,fmt1) '1       255     112     0       1.1     255     112     0  '
	Write(9,fmt1) '1.1     255     112     0       1.2     255     112     0  '
	Write(9,fmt1) '1.2     255     77      0       1.3     255     77      0  '
	Write(9,fmt1) '1.3     255     77      0       1.4     255     77      0  '
	Write(9,fmt1) '1.4     255     0       0       1.5     255     0       0  '
 	Write(9,fmt2) 'B       32      96      255'
 	Write(9,fmt2) 'F       255     0         0'
 	Write(9,fmt2) 'N       128     128     128'
 	close(9)
!
 	Elseif(table_name=='pale-europe.cpt')  then
! 
 	open(9,file=table_name,status='unknown') 
	fmt1='(a59)'
	fmt2='(a27)'
	Write(9,fmt1) '#----------------------------------------------------------'
 	Write(9,fmt1) '#      cpt created by: makecpt -Cno_green -T-8/8/1         '
	Write(9,fmt1) '#----------------------------------------------------------'	
 	Write(9,fmt1) '-8      32      96      255     -7      32      96      255'
 	Write(9,fmt1) '-7      32      159     255     -6      32      159     255'
 	Write(9,fmt1) '-6      32      191     255     -5      32      191     255'
 	Write(9,fmt1) '-5      0       207     255     -4      0       207     255'
 	Write(9,fmt1) '-4      42      255     255     -3      42      255     255'
 	Write(9,fmt1) '-3      85      255     255     -2      85      255     255'
 	Write(9,fmt1) '-2      127     255     255     -1      127     255     255'
 	Write(9,fmt1) '-1      170     255     255     0       170     255     255'
 	Write(9,fmt1) '0       255     255     84      1       255     255     84 '
 	Write(9,fmt1) '1       255     240     0       2       255     240     0  '
 	Write(9,fmt1) '2       255     191     0       3       255     191     0  '
 	Write(9,fmt1) '3       255     168     0       4       255     168     0  '
 	Write(9,fmt1) '4       255     138     0       5       255     138     0  '
 	Write(9,fmt1) '5       255     112     0       6       255     112     0  '
 	Write(9,fmt1) '6       255     77      0       7       255     77      0  '
 	Write(9,fmt1) '7       255     0       0       8       255     0       0  '
 	Write(9,fmt2) 'B       32      96      255'
 	Write(9,fmt2) 'F       255     0         0'
 	Write(9,fmt2) 'N       128     128     128'
 	close(9) 
! 
 	Elseif(table_name=='pale-fennoscandia-su.cpt'.or.    & 
	       table_name=='pale-greenland-su.cpt'.or. &  
	       table_name=='pale-antarctica-su.cpt')       then
!
 	open(9,file=table_name,status='unknown') 
	fmt1='(a59)'
	fmt2='(a27)'
	Write(9,fmt1) '#----------------------------------------------------------'
 	Write(9,fmt1) '#      cpt created by: makecpt -Cno_green -T-10/10/1       '
	Write(9,fmt1) '#----------------------------------------------------------'	
 	Write(9,fmt1) '-10     32      96      255     -9      32      96      255'
 	Write(9,fmt1) '-9      32      96      255     -8      32      96      255'
 	Write(9,fmt1) '-8      32      159     255     -7      32      159     255'
 	Write(9,fmt1) '-7      32      191     255     -6      32      191     255'
 	Write(9,fmt1) '-6      0       207     255     -5      0       207     255'
 	Write(9,fmt1) '-5      42      255     255     -4      42      255     255'
 	Write(9,fmt1) '-4      42      255     255     -3      42      255     255'
 	Write(9,fmt1) '-3      85      255     255     -2      85      255     255'
 	Write(9,fmt1) '-2      127     255     255     -1      127     255     255'
 	Write(9,fmt1) '-1      170     255     255     0       170     255     255'
 	Write(9,fmt1) '0       255     255     84      1       255     255     84 '
 	Write(9,fmt1) '1       255     255     84      2       255     255     84 '
 	Write(9,fmt1) '2       255     240     0       3       255     240     0  '
 	Write(9,fmt1) '3       255     191     0       4       255     191     0  '
 	Write(9,fmt1) '4       255     168     0       5       255     168     0  '
 	Write(9,fmt1) '5       255     138     0       6       255     138     0  '
 	Write(9,fmt1) '6       255     138     0       7       255     138     0  '
 	Write(9,fmt1) '7       255     112     0       8       255     112     0  '
 	Write(9,fmt1) '8       255     77      0       9       255     77      0  '
 	Write(9,fmt1) '9       255     0       0       10      255     0       0  '
 	Write(9,fmt2) 'B       32      96      255'
 	Write(9,fmt2) 'F       255     0         0'
 	Write(9,fmt2) 'N       128     128     128'
 	close(9) 
!
 	Elseif(table_name=='pale-fennoscandia-n.cpt'.or.   & 
	       table_name=='pale-antarctica-n.cpt')    then
!
 	open(9,file=table_name,status='unknown') 
	fmt1='(a59)'
	fmt2='(a27)'
	Write(9,fmt1) '#----------------------------------------------------------'
 	Write(9,fmt1) '#      cpt created by: makecpt -Cno_green -T-2/2/0.2       '
	Write(9,fmt1) '#----------------------------------------------------------'	
 	Write(9,fmt1) '-2      32      96      255     -1.8    32      96      255'
 	Write(9,fmt1) '-1.8    32      96      255     -1.6    32      96      255'
 	Write(9,fmt1) '-1.6    32      159     255     -1.4    32      159     255'
 	Write(9,fmt1) '-1.4    32      191     255     -1.2    32      191     255'
 	Write(9,fmt1) '-1.2    0       207     255     -1      0       207     255'
 	Write(9,fmt1) '-1      42      255     255     -0.8    42      255     255'
 	Write(9,fmt1) '-0.8    42      255     255     -0.6    42      255     255'
 	Write(9,fmt1) '-0.6    85      255     255     -0.4    85      255     255'
 	Write(9,fmt1) '-0.4    127     255     255     -0.2    127     255     255'
 	Write(9,fmt1) '-0.2    170     255     255     0       170     255     255'
 	Write(9,fmt1) '0       255     255     84      0.2     255     255     84 '
 	Write(9,fmt1) '0.2     255     255     84      0.4     255     255     84 '
 	Write(9,fmt1) '0.4     255     240     0       0.6     255     240     0  '
 	Write(9,fmt1) '0.6     255     191     0       0.8     255     191     0  '
 	Write(9,fmt1) '0.8     255     168     0       1       255     168     0  '
 	Write(9,fmt1) '1       255     138     0       1.2     255     138     0  '
 	Write(9,fmt1) '1.2     255     138     0       1.4     255     138     0  '
 	Write(9,fmt1) '1.4     255     112     0       1.6     255     112     0  '
 	Write(9,fmt1) '1.6     255     77      0       1.8     255     77      0  '
 	Write(9,fmt1) '1.8     255     0       0       2       255     0       0  '
 	Write(9,fmt2) 'B       32      96      255'
 	Write(9,fmt2) 'F       255     0         0'
 	Write(9,fmt2) 'N       128     128     128'
 	close(9) 
!
 	Elseif(table_name=='pale-greenland-n.cpt')  then
!
 	open(9,file=table_name,status='unknown') 
	fmt1='(a59)'
	fmt2='(a27)'
	Write(9,fmt1) '#----------------------------------------------------------'
 	Write(9,fmt1) '#      cpt created by: makecpt -Cno_green -T-1/1/0.2       '
	Write(9,fmt1) '#----------------------------------------------------------'	
 	Write(9,fmt1) '-1      32      96      255     -0.8    32      96      255'
 	Write(9,fmt1) '-0.8    32      159     255     -0.6    32      159     255'
 	Write(9,fmt1) '-0.6    0       207     255     -0.4    0       207     255'
 	Write(9,fmt1) '-0.4    42      255     255     -0.2    42      255     255'
 	Write(9,fmt1) '-0.2    127     255     255     0       127     255     255'
 	Write(9,fmt1) '0       255     255     84      0.2     255     255     84 '
 	Write(9,fmt1) '0.2     255     240     0       0.4     255     240     0  '
 	Write(9,fmt1) '0.4     255     168     0       0.6     255     168     0  '
 	Write(9,fmt1) '0.6     255     138     0       0.8     255     138     0  '
 	Write(9,fmt1) '0.8     255     77      0       1       255     77      0  '
 	Write(9,fmt2) 'B       32      96      255'
 	Write(9,fmt2) 'F       255     0         0'
 	Write(9,fmt2) 'N       128     128     128'
 	close(9) 
!
 	Elseif(table_name=='pale-namerica-su.cpt')  then
!
 	open(9,file=table_name,status='unknown') 
	fmt1='(a59)'
	fmt2='(a27)'
	Write(9,fmt1) '#----------------------------------------------------------'
 	Write(9,fmt1) '#      cpt created by: makecpt -Cno_green -T-24/24/2       '
	Write(9,fmt1) '#----------------------------------------------------------'	
 	Write(9,fmt1) '-24     32      96      255     -22     32      96      255'
 	Write(9,fmt1) '-22     32      96      255     -20     32      96      255'
 	Write(9,fmt1) '-20     32      159     255     -18     32      159     255'
 	Write(9,fmt1) '-18     32      191     255     -16     32      191     255'
 	Write(9,fmt1) '-16     32      191     255     -14     32      191     255'
 	Write(9,fmt1) '-14     0       207     255     -12     0       207     255'
 	Write(9,fmt1) '-12     42      255     255     -10     42      255     255'
 	Write(9,fmt1) '-10     42      255     255     -8      42      255     255'
 	Write(9,fmt1) '-8      85      255     255     -6      85      255     255'
 	Write(9,fmt1) '-6      127     255     255     -4      127     255     255'
 	Write(9,fmt1) '-4      127     255     255     -2      127     255     255'
 	Write(9,fmt1) '-2      170     255     255     0       170     255     255'
 	Write(9,fmt1) '0       255     255     84      2       255     255     84 '
 	Write(9,fmt1) '2       255     255     84      4       255     255     84 '
 	Write(9,fmt1) '4       255     240     0       6       255     240     0  '
 	Write(9,fmt1) '6       255     191     0       8       255     191     0  '
 	Write(9,fmt1) '8       255     191     0       10      255     191     0  '
 	Write(9,fmt1) '10      255     168     0       12      255     168     0  '
 	Write(9,fmt1) '12      255     138     0       14      255     138     0  '
 	Write(9,fmt1) '14      255     138     0       16      255     138     0  '
	Write(9,fmt1) '16      255     112     0       18      255     112     0  '
	Write(9,fmt1) '18      255     77      0       20      255     77      0  '
	Write(9,fmt1) '20      255     77      0       22      255     77      0  '
	Write(9,fmt1) '22      255     0       0       24      255     0       0  ' 
 	Write(9,fmt2) 'B       32      96      255'
 	Write(9,fmt2) 'F       255     0         0'
 	Write(9,fmt2) 'N       128     128     128'
 	close(9) 
!
 	Elseif(table_name=='pale-namerica-n.cpt')  then
!
 	open(9,file=table_name,status='unknown') 
	fmt1='(a59)'
	fmt2='(a27)'
	Write(9,fmt1) '#----------------------------------------------------------'
 	Write(9,fmt1) '#      cpt created by: makecpt -Cno_green -T-10/10/1       '
	Write(9,fmt1) '#----------------------------------------------------------'	
 	Write(9,fmt1) '-10     32      96      255     -9      32      96      255'
 	Write(9,fmt1) '-9      32      96      255     -8      32      96      255'
 	Write(9,fmt1) '-8      32      159     255     -7      32      159     255'
 	Write(9,fmt1) '-7      32      191     255     -6      32      191     255'
 	Write(9,fmt1) '-6      0       207     255     -5      0       207     255'
 	Write(9,fmt1) '-5      42      255     255     -4      42      255     255'
 	Write(9,fmt1) '-4      42      255     255     -3      42      255     255'
 	Write(9,fmt1) '-3      85      255     255     -2      85      255     255'
 	Write(9,fmt1) '-2      127     255     255     -1      127     255     255'
 	Write(9,fmt1) '-1      170     255     255     0       170     255     255'
 	Write(9,fmt1) '0       255     255     84      1       255     255     84 '
 	Write(9,fmt1) '1       255     255     84      2       255     255     84 '
 	Write(9,fmt1) '2       255     240     0       3       255     240     0  '
 	Write(9,fmt1) '3       255     191     0       4       255     191     0  '
 	Write(9,fmt1) '4       255     168     0       5       255     168     0  '
 	Write(9,fmt1) '5       255     138     0       6       255     138     0  '
 	Write(9,fmt1) '6       255     138     0       7       255     138     0  '
 	Write(9,fmt1) '7       255     112     0       8       255     112     0  '
 	Write(9,fmt1) '8       255     77      0       9       255     77      0  '
 	Write(9,fmt1) '9       255     0       0       10      255     0       0  ' 
 	Write(9,fmt2) 'B       32      96      255'
 	Write(9,fmt2) 'F       255     0         0'
 	Write(9,fmt2) 'N       128     128     128'
 	close(9) 
!
 Endif
!
 END SUBROUTINE COLOR_TABLES
!
!
!
!
!
!
 SUBROUTINE MAKE_ELAREB_MAPS_GLAC      (option_reb_sm, & 
 	  			        resolution,    & 
	   				nv,            & 
	  				code,          & 
	  				iter,          & 
	  				mode,          &
	  				degree,        & 
	  				vstring,       & 
					TITLICE,       & 
					FILE_GMT) 
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! Prepares a GMT script for *** GLOBAL maps *** showing variations of S, U & N in 
! response to melting of "Small Glaciers" spread over the globe *** GS Aug 2010 ** 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
 IMPLICIT NONE
 CHARACTER*20  FILE_GMT
 CHARACTER*1   OPTION_REB_SM
 CHARACTER*30  NAMEIN, NAMEOUT, NAMEF   
 CHARACTER*80  TITRE
 INTEGER I 

 CHARACTER*100 SHORT_VISCO 
 CHARACTER*10  TITLICE, RESOLUTION, DEGREE
 CHARACTER*30  VSTRING
 CHARACTER*20  R_OPTION, T_OPTION, J_OPTION, B_OPTION
 CHARACTER*3   NV, CODE
 CHARACTER*1   OPTION_ROF
 CHARACTER*1   ITER, MODE 

 INTEGER K
!
!
!
!
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\ 
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
! A  file for global analysis 
!
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\ 
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\ 
!
 OPEN (9,FILE=FILE_GMT,STATUS='UNKNOWN')
! 
 Write(9,*)"gmtset PAPER_MEDIA A4+"
 Write(9,*)"gmtset LABEL_FONT_SIZE 24p"
 Write(9,*)"gmtset ANOT_FONT_SIZE 16p"
 Write(9,*)"gmtset FRAME_WIDTH 0.1c"
 Write(9,*)" "
!
!
!




!
! ///////////////////////////////////////// ---------- Small glaciers
!        Global map   
! ///////////////////////////////////////// ---------- Small glaciers
! 
 IF(OPTION_REB_SM=='y') THEN 
! 
  T_OPTION="-T-0.4/0.4/0.05"         ! Range of the palette 
  R_OPTION="-R0/360/-80/80"     ! Range of plot  
  B_OPTION="-B1f0.1a0.2/:mm/yr:"       ! B option for "psscale", S, U, and N 
! 
 Write(9,*) "# "
 Write(9,*) "# "
 Write(9,*) "# "
 Write(9,*) "# GLOBAL MAP - GLOBAL MAP - GLOBAL MAP - GLOBAL MAP - GLOBAL MAP "
 Write(9,*) "# GLOBAL MAP - GLOBAL MAP - GLOBAL MAP - GLOBAL MAP - GLOBAL MAP "
 Write(9,*) "# GLOBAL MAP - GLOBAL MAP - GLOBAL MAP - GLOBAL MAP - GLOBAL MAP "
 Write(9,*) "# "
!
 DO 10 I=1, 3 
!
 if(i==1)namein ="smap_glac_glob.dat" 
 if(i==2)namein ="umap_glac_glob.dat" 
 if(i==3)namein ="nmap_glac_glob.dat" 
 if(i==1)nameout="smap_glac_glob.ps" 
 if(i==2)nameout="umap_glac_glob.ps" 
 if(i==3)nameout="nmap_glac_glob.ps" 
 if(i==1)namef="glc_tmpfs.dat"
 if(i==2)namef="glc_tmpfu.dat" 
 if(i==3)namef="glc_tmpfn.dat"  
 if(i==1)TITRE="Rate of relative SL change (\dot S)"
 if(i==2)TITRE="Rate of vertical uplift (\dot U)"
 if(i==3)TITRE="Rate of absolute SL change (\dot N)"
!
 Write(9,*) " "
 if(i==1)Write(9,*)"# ---- Global map of S at present time for ELASTIC REBOUND small glaciers ----" 		      	   
 if(i==2)Write(9,*)"# ---- Global map of U at present time for ELASTIC REBOUND small glaciers ----" 		      	   
 if(i==3)Write(9,*)"# ---- Global map of N at present time for ELASTIC REBOUND small glaciers ----" 		      	   
!
 if(i==1)Write(9,*) "echo", "    - S for small glaciers - global" 
 if(i==2)Write(9,*) "echo", "    - U for small glaciers - global" 
 if(i==3)Write(9,*) "echo", "    - N for small glaciers - global" 
!
 Write(9,*) " "
 Write(9,*) "makecpt -Cno_green ", trim(adjustl(T_OPTION)), " > pale.cpt" 
 Write(9,*) "psbasemap -X3 -Y5 -Ba180/a40WSEn -Jm0.018i ", & 
             trim(adjustl(R_OPTION)), " -K > ",            & 
	     trim(adjustl(nameout))  
 Write(9,*) "pscontour -I -Jm -O -K ", & 
             trim(adjustl(R_OPTION)), " ", & 
             trim(adjustl(namein)), " -Cpale.cpt  >> ", &   
	     trim(adjustl(nameout)) 
 Write(9,*) "pscoast -Jm -Dc -B -W2/0  -A1000 -O -K ",  &
             trim(adjustl(R_OPTION)), " >> ",           & 
             trim(adjustl(nameout))  
!
! A file for "titles" 
 open (4,file=NAMEF,status='unknown')  
!
! Alma 
 If(CODE=='-1')  THEN 
		Write(4,*) "180 +84.2 26 0 3 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	 	  "  -Elasticity: ALMA ", trim(adjustl(SHORT_VISCO))  
		Write(4,*) "180 -87.5   13  0 0 BC ",          & 
 	   	 	   " -LMAX=", trim(adjustl(DEGREE)),   & 
 	    		   " -RES=", trim(adjustl(RESOLUTION)),& 
	    	 	   " -NV=", trim(adjustl(NV)),         & 
	    		   " -MODE=", trim(adjustl(MODE)),     &   
	    	 	   " -ITER=", trim(adjustl(ITER))
		close(4) 
!
! Taboo 
 ELSEIF(CODE/='-1'.and.CODE/='-2')THEN 
		Write(4,*) "180 +84.2 26 0 3 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	   	"  -Elasticity: TABOO "  
		Write(4,*) "180 -87.5   13  0 0 BC ",          & 
 	    		   " -LMAX=", trim(adjustl(DEGREE)),   & 
 	    		   " -RES=", trim(adjustl(RESOLUTION)),& 
	    		   " -NV=", trim(adjustl(NV)),         & 
	    		   " -CODE=", trim(adjustl(CODE)),     & 
	    		   " -MODE=", trim(adjustl(MODE)),     &   
	    		   " -ITER=", trim(adjustl(ITER))
		close(4) 
! External 
 ELSEIF(CODE=='-2')THEN 
		Write(4,*) "180 +84.2 26 0 3 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	   	"  -Elasticity: External "  
		Write(4,*) "180 -87.5   13  0 0 BC ",          & 
 	    		   " -LMAX=", trim(adjustl(DEGREE)),   & 
 	    		   " -RES=", trim(adjustl(RESOLUTION)),& 
	    		   " -CODE=", trim(adjustl(CODE)),     & 
	    		   " -MODE=", trim(adjustl(MODE)),     &   
	    		   " -ITER=", trim(adjustl(ITER))
		close(4) 
 ENDIF
		Write(9,*) "pstext -N -R -Jm -B ", trim(adjustl(NAMEF)),    & 
		           " -G0 -O -K >> ", trim(adjustl(nameout)) 
		Write(9,*) "psscale -U/0.5/0.5/'SELEN 3.2' -E -Cpale.cpt ", & 
		            trim(adjustl(B_OPTION)), " -D8.25/-1/10/1h -O -K >> ",   & 
	                    trim(adjustl(nameout)) 
 		Write(9,*) "psbasemap -J -R -B::/::EWSN -O >> ", trim(adjustl(nameout)) 
 10 CONTINUE
!
!
!
! \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/
!   A script for the Normalized Sea Level change 
! \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/
! \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/
!   A script for the Normalized Sea Level change 
! \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/
! \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/
!   A script for the Normalized Sea Level change 
! \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/

!
!
 Write(9,*) "# "
 Write(9,*) "# "
 Write(9,*) "# "
 Write(9,*) "# GLOBAL NORMALIZED SEA LEVEL CHANGE MAP - GLOBAL NORMALIZED SEA LEVEL CHANGE MAP "
 Write(9,*) "# GLOBAL NORMALIZED SEA LEVEL CHANGE MAP - GLOBAL NORMALIZED SEA LEVEL CHANGE MAP "
 Write(9,*) "# GLOBAL NORMALIZED SEA LEVEL CHANGE MAP - GLOBAL NORMALIZED SEA LEVEL CHANGE MAP "
 Write(9,*) "# "
!
 namein  = "smap_glac_glob_norm.dat" 
 nameout = "smap_glac_glob_norm.ps" 
 namef   = "glc_ngtmpfs.dat"  
 TITRE   = "Sea Level Change / Eustatic Change"
!
 Write(9,*)" "
 Write(9,*)"gmtset PAPER_MEDIA A4+"
 Write(9,*)"gmtset LABEL_FONT_SIZE 24p"
 Write(9,*)"gmtset ANOT_FONT_SIZE 16p"
 Write(9,*)"gmtset FRAME_WIDTH 0.1c"
 Write(9,*)" "
!
 T_OPTION="-T-0.1/1.3/0.1"         ! Range of the palette 
 R_OPTION="-R0/360/-80/85"         ! Range of plot  
!
 Write(9,*) "echo S normalized for small glaciers - global" 
!
 Write(9,*) " "
 Write(9,*) "makecpt -Cno_green ", trim(adjustl(T_OPTION)), " > pale_norm.cpt" 
 Write(9,*) "psbasemap -X3 -Y5 -JQ0/18 -Ba90/a30WSEn  ", & 
             trim(adjustl(R_OPTION)), " -K > ",            & 
	     trim(adjustl(nameout))  
 Write(9,*) "pscontour -I -O -K  -JQ ",         & 
             trim(adjustl(R_OPTION)), " ",      & 
             trim(adjustl(namein)), " -Cpale_norm.cpt  >> ", &   
	     trim(adjustl(nameout)) 
 Write(9,*) "pscoast -Dc -B -W2/0 -JQ -A10000 -G200 -O -K ", & 
             trim(adjustl(R_OPTION)), " >> ",      & 
             trim(adjustl(nameout))  
!
!
! A file for "titles" 
 open (4,file=NAMEF,status='unknown')  
!
! Alma 
 If(CODE=='-1')  THEN 
		Write(4,*) "0 +100     22 0 0 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "0 -175     14 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	 	  "  -Elasticity: ALMA ", trim(adjustl(SHORT_VISCO))  
		Write(4,*) "0 -160.5   13  0 0 BC  ",          & 
 	   	 	   " -LMAX=", trim(adjustl(DEGREE)),   & 
 	    		   " -RES=", trim(adjustl(RESOLUTION)),& 
	    	 	   " -NV=", trim(adjustl(NV)),         & 
	    		   " -MODE=", trim(adjustl(MODE)),     &   
	    	 	   " -ITER=", trim(adjustl(ITER))
		close(4) 
!
! Taboo 
 ELSEIF(CODE/='-1'.and.CODE/='-2')THEN 
		Write(4,*) "0 +100     22 0 0 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "0 -175     14 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	   	"  -Elasticity: TABOO "  
		Write(4,*) "0 -160.5   13  0 0 BC ",          & 
 	    		   " -LMAX=", trim(adjustl(DEGREE)),   & 
 	    		   " -RES=", trim(adjustl(RESOLUTION)),& 
	    		   " -NV=", trim(adjustl(NV)),         & 
	    		   " -CODE=", trim(adjustl(CODE)),     & 
	    		   " -MODE=", trim(adjustl(MODE)),     &   
	    		   " -ITER=", trim(adjustl(ITER))
		close(4) 
! External 
 ELSEIF(CODE=='-2')THEN 
		Write(4,*) "0 +100     22 0 0 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "0 -175     14 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	   	"  -Elasticity: External "  
		Write(4,*) "0 -160.5   13  0 0 BC ",          & 
 	    		   " -LMAX=", trim(adjustl(DEGREE)),   & 
 	    		   " -RES=", trim(adjustl(RESOLUTION)),& 
	    		   " -CODE=", trim(adjustl(CODE)),     & 
	    		   " -MODE=", trim(adjustl(MODE)),     &   
	    		   " -ITER=", trim(adjustl(ITER))
		close(4) 
 ENDIF
!
!
!
 Write(9,*) "pstext -N -R -JQ -B ", trim(adjustl(NAMEF)),    & 
            " -G0 -O -K >> ", trim(adjustl(nameout)) 
 Write(9,*) "psscale -U/0.5/0.5/'SELEN 3.2' -E -Cpale_norm.cpt ", & 
		           "-B1f0.1a0.1/:S/S@-E@-: -D9/-1.5/-16/1h -O -K >> ",   & 
	                    trim(adjustl(nameout)) 
 Write(9,*) "psbasemap -JQ -R -B::/::EWSN -O >> ", trim(adjustl(nameout)) 


! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ---------- Global Glaciers!
 ENDIF    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ---------- Global Glaciers!
!
 END SUBROUTINE MAKE_ELAREB_MAPS_GLAC
!
!
!
!
!
!
 SUBROUTINE MAKE_ELAREB_MAPS_GREENLAND (option_reb_gg, & 
	  			        option_reb_gr, & 
 	  			        resolution,    & 
	   				nv,            & 
	  				code,          & 
	  				iter,          & 
	  				mode,          &
	  				degree,        & 
	  				vstring,       & 
					TITLICE,       & 
					FILE_GMT) 
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! Prepares a GMT script for *** GLOBAL and REGIONAL maps *** showing variations of 
! of S, U & N in response to melting of the Greenland ice sheet *** GS July 2010 ** 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! Based on "SUBROUTINE MAKE_GMAPS"
! Last change GS March 22 2010 
!
!


 IMPLICIT NONE
 CHARACTER*20  FILE_GMT
 CHARACTER*1   OPTION_REB_GG, OPTION_REB_GR
 CHARACTER*30  NAMEIN, NAMEOUT, NAMEF   
 CHARACTER*80  TITRE
 INTEGER I 

 CHARACTER*100 SHORT_VISCO 
 CHARACTER*10  TITLICE, RESOLUTION, DEGREE
 CHARACTER*30  VSTRING
 CHARACTER*20  R_OPTION, T_OPTION, J_OPTION
 CHARACTER*3   NV, CODE
 CHARACTER*1   OPTION_ROF
 CHARACTER*1   ITER, MODE 

 INTEGER K
!
!
!
!
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
! A unique file for global and regional analyses 
!
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
 OPEN (9,FILE=FILE_GMT,STATUS='UNKNOWN')
! 
 Write(9,*)"gmtset PAPER_MEDIA A4+"
 Write(9,*)"gmtset LABEL_FONT_SIZE 24p"
 Write(9,*)"gmtset ANOT_FONT_SIZE 16p"
 Write(9,*)"gmtset FRAME_WIDTH 0.1c"
 Write(9,*)" "
!
!
!




!
! ///////////////////////////////////////// ---------- Global Greenland!
!        Global map for Greenland 
! ///////////////////////////////////////// ---------- Global Greenland!
! 
 IF(OPTION_REB_GG=='y') THEN 
! 
  T_OPTION="-T-0.4/0.4/0.05"         ! Range of the palette 
  R_OPTION="-R0/360/-80/80"     ! Range of plot  
! 
 Write(9,*) "# "
 Write(9,*) "# "
 Write(9,*) "# "
 Write(9,*) "# GLOBAL MAP - GLOBAL MAP - GLOBAL MAP - GLOBAL MAP - GLOBAL MAP "
 Write(9,*) "# GLOBAL MAP - GLOBAL MAP - GLOBAL MAP - GLOBAL MAP - GLOBAL MAP "
 Write(9,*) "# GLOBAL MAP - GLOBAL MAP - GLOBAL MAP - GLOBAL MAP - GLOBAL MAP "
 Write(9,*) "# "
!
 DO 10 I=1, 3 
!
 if(i==1)namein ="smap_gree_glob.dat" 
 if(i==2)namein ="umap_gree_glob.dat" 
 if(i==3)namein ="nmap_gree_glob.dat" 
 if(i==1)nameout="smap_gree_glob.ps" 
 if(i==2)nameout="umap_gree_glob.ps" 
 if(i==3)nameout="nmap_gree_glob.ps" 
 if(i==1)namef="gtmpfs.dat"
 if(i==2)namef="gtmpfu.dat" 
 if(i==3)namef="gtmpfn.dat"  
 if(i==1)TITRE="Rate of relative SL change (\dot S)"
 if(i==2)TITRE="Rate of vertical uplift (\dot U)"
 if(i==3)TITRE="Rate of absolute SL change (\dot N)"
!
 Write(9,*) " "
 if(i==1)Write(9,*)"# ---- Global map of S at present time for ELASTIC REBOUND in GREENLAND ----" 		      	   
 if(i==2)Write(9,*)"# ---- Global map of U at present time for ELASTIC REBOUND in GREENLAND ----" 		      	   
 if(i==3)Write(9,*)"# ---- Global map of N at present time for ELASTIC REBOUND in GREENLAND ----" 		      	   
!
 if(i==1)Write(9,*) "echo", "    - S for Greenland - global" 
 if(i==2)Write(9,*) "echo", "    - U for Greenland - global" 
 if(i==3)Write(9,*) "echo", "    - N for Greenland - global" 
!
 Write(9,*) " "
 Write(9,*) "makecpt -Cno_green ", trim(adjustl(T_OPTION)), " > pale.cpt" 
 Write(9,*) "psbasemap -X3 -Y5 -Ba180/a40WSEn -Jm0.018i ", & 
             trim(adjustl(R_OPTION)), " -K > ",            & 
	     trim(adjustl(nameout))  
!
 Write(9,*) "pscontour -I -Jm -O -K ", & 
             trim(adjustl(R_OPTION)), " ", & 
             trim(adjustl(namein)), " -Cpale.cpt  >> ", &   
	     trim(adjustl(nameout)) 
!
 Write(9,*) "pscontour -O -K -G8 -W2/0/100/0 -Jm -A+g255 ",         & 
             trim(adjustl(R_OPTION)), " ", & 
             trim(adjustl(namein)), " -Cpale.cpt  >> ", & 
	     trim(adjustl(nameout)) 
!
 Write(9,*) "pscoast -Jm -Dc -B -W2/0  -A1000 -O -K ",  &
             trim(adjustl(R_OPTION)), " >> ",           & 
             trim(adjustl(nameout))  
!
! A file for "titles" 
 open (4,file=NAMEF,status='unknown')  
!
! Alma 
 If(CODE=='-1')  THEN 
		Write(4,*) "180 +84.2 26 0 3 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	 	  "  -Elasticity: ALMA ", trim(adjustl(SHORT_VISCO))  
		Write(4,*) "180 -87.5   13  0 0 BC ",          & 
 	   	 	   " -LMAX=", trim(adjustl(DEGREE)),   & 
 	    		   " -RES=", trim(adjustl(RESOLUTION)),& 
	    	 	   " -NV=", trim(adjustl(NV)),         & 
	    		   " -MODE=", trim(adjustl(MODE)),     &   
	    	 	   " -ITER=", trim(adjustl(ITER))
		close(4) 
!
! Taboo 
 ELSEIF(CODE/='-1'.and.CODE/='-2')THEN 
		Write(4,*) "180 +84.2 26 0 3 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	   	"  -Elasticity: TABOO "  
		Write(4,*) "180 -87.5   13  0 0 BC ",          & 
 	    		   " -LMAX=", trim(adjustl(DEGREE)),   & 
 	    		   " -RES=", trim(adjustl(RESOLUTION)),& 
	    		   " -NV=", trim(adjustl(NV)),         & 
	    		   " -CODE=", trim(adjustl(CODE)),     & 
	    		   " -MODE=", trim(adjustl(MODE)),     &   
	    		   " -ITER=", trim(adjustl(ITER))
		close(4) 
! External 
 ELSEIF(CODE=='-2')THEN 
		Write(4,*) "180 +84.2 26 0 3 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	   	"  -Elasticity: External "  
		Write(4,*) "180 -87.5   13  0 0 BC ",          & 
 	    		   " -LMAX=", trim(adjustl(DEGREE)),   & 
 	    		   " -RES=", trim(adjustl(RESOLUTION)),& 
	    		   " -CODE=", trim(adjustl(CODE)),     & 
	    		   " -MODE=", trim(adjustl(MODE)),     &   
	    		   " -ITER=", trim(adjustl(ITER))
		close(4) 
 ENDIF
		Write(9,*) "pstext -N -R -Jm -B ", trim(adjustl(NAMEF)),    & 
		           " -G0 -O -K >> ", trim(adjustl(nameout)) 
		Write(9,*) "psscale -U/0.5/0.5/'SELEN 3.2' -E -Cpale.cpt ", & 
		           "-B1f0.1a0.2/:mm/yr: -D8.25/-1/10/1h -O -K >> ",   & 
	                    trim(adjustl(nameout)) 
 		Write(9,*) "psbasemap -J -R -B::/::EWSN -O >> ", trim(adjustl(nameout)) 
 10 CONTINUE
!
!
!
! \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/
!   A script for the Normalized Sea Level change 
! \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/
! \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/
!   A script for the Normalized Sea Level change 
! \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/
! \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/
!   A script for the Normalized Sea Level change 
! \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/

!
!
 Write(9,*) "# "
 Write(9,*) "# "
 Write(9,*) "# "
 Write(9,*) "# GLOBAL NORMALIZED SEA LEVEL CHANGE MAP - GLOBAL NORMALIZED SEA LEVEL CHANGE MAP "
 Write(9,*) "# GLOBAL NORMALIZED SEA LEVEL CHANGE MAP - GLOBAL NORMALIZED SEA LEVEL CHANGE MAP "
 Write(9,*) "# GLOBAL NORMALIZED SEA LEVEL CHANGE MAP - GLOBAL NORMALIZED SEA LEVEL CHANGE MAP "
 Write(9,*) "# "
!
 namein  = "smap_gree_glob_norm.dat" 
 nameout = "smap_gree_glob_norm.ps" 
 namef   = "ngtmpfs.dat"  
 TITRE   = "Sea Level Change / Eustatic Change"
!
 Write(9,*)" "
 Write(9,*)"gmtset PAPER_MEDIA A4+"
 Write(9,*)"gmtset LABEL_FONT_SIZE 24p"
 Write(9,*)"gmtset ANOT_FONT_SIZE 16p"
 Write(9,*)"gmtset FRAME_WIDTH 0.1c"
 Write(9,*)" "
!
 T_OPTION="-T-0.1/1.3/0.1"         ! Range of the palette 
 R_OPTION="-R0/360/-80/85"         ! Range of plot  
!
 Write(9,*) "echo", "    - S normalized for  Greenland - global" 
!
 Write(9,*) " "
 Write(9,*) "makecpt -Cno_green ", trim(adjustl(T_OPTION)), " > pale_norm.cpt" 
 Write(9,*) "psbasemap -X3 -Y5 -JQ0/18 -Ba90/a30WSEn  ", & 
             trim(adjustl(R_OPTION)), " -K > ",            & 
	     trim(adjustl(nameout))  
 Write(9,*) "pscontour -I -O -K  -JQ ",         & 
             trim(adjustl(R_OPTION)), " ",      & 
             trim(adjustl(namein)), " -Cpale_norm.cpt  >> ", &   
	     trim(adjustl(nameout)) 
 Write(9,*) "pscoast -Dc -B -W2/0 -JQ -A10000 -G200 -O -K ", & 
             trim(adjustl(R_OPTION)), " >> ",      & 
             trim(adjustl(nameout))  
!
!
! A file for "titles" 
 open (4,file=NAMEF,status='unknown')  
!
! Alma 
 If(CODE=='-1')  THEN 
		Write(4,*) "0 +100     22 0 0 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "0 -175     14 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	 	  "  -Elasticity: ALMA ", trim(adjustl(SHORT_VISCO))  
		Write(4,*) "0 -160.5   13  0 0 BC  ",          & 
 	   	 	   " -LMAX=", trim(adjustl(DEGREE)),   & 
 	    		   " -RES=", trim(adjustl(RESOLUTION)),& 
	    	 	   " -NV=", trim(adjustl(NV)),         & 
	    		   " -MODE=", trim(adjustl(MODE)),     &   
	    	 	   " -ITER=", trim(adjustl(ITER))
		close(4) 
!
! Taboo 
 ELSEIF(CODE/='-1'.and.CODE/='-2')THEN 
		Write(4,*) "0 +100     22 0 0 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "0 -175     14 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	   	"  -Elasticity: TABOO "  
		Write(4,*) "0 -160.5   13  0 0 BC ",          & 
 	    		   " -LMAX=", trim(adjustl(DEGREE)),   & 
 	    		   " -RES=", trim(adjustl(RESOLUTION)),& 
	    		   " -NV=", trim(adjustl(NV)),         & 
	    		   " -CODE=", trim(adjustl(CODE)),     & 
	    		   " -MODE=", trim(adjustl(MODE)),     &   
	    		   " -ITER=", trim(adjustl(ITER))
		close(4) 
! External 
 ELSEIF(CODE=='-2')THEN 
		Write(4,*) "0 +100     22 0 0 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "0 -175     14 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	   	"  -Elasticity: External "  
		Write(4,*) "0 -160.5   13  0 0 BC ",          & 
 	    		   " -LMAX=", trim(adjustl(DEGREE)),   & 
 	    		   " -RES=", trim(adjustl(RESOLUTION)),& 
	    		   " -CODE=", trim(adjustl(CODE)),     & 
	    		   " -MODE=", trim(adjustl(MODE)),     &   
	    		   " -ITER=", trim(adjustl(ITER))
		close(4) 
 ENDIF
!
!
!
 Write(9,*) "pstext -N -R -JQ -B ", trim(adjustl(NAMEF)),    & 
            " -G0 -O -K >> ", trim(adjustl(nameout)) 
 Write(9,*) "psscale -U/0.5/0.5/'SELEN 3.2' -E -Cpale_norm.cpt ", & 
		           "-B1f0.1a0.1/:S/S@-E@-: -D9/-1.5/-16/1h -O -K >> ",   & 
	                    trim(adjustl(nameout)) 
 Write(9,*) "psbasemap -JQ -R -B::/::EWSN -O >> ", trim(adjustl(nameout)) 


! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ---------- Global Greenland!
 ENDIF    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ---------- Global Greenland!
!
!
!
 
 
!
! ///////////////////////////////////////// ---------- Regional Greenland!
!        Regional map for Greenland 
! ///////////////////////////////////////// ---------- Regional Greenland!
! 
 IF(OPTION_REB_GR=='y') THEN 
!
! Here I use a GMT projection suggested by Louise Sorensen  
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! July 2010
!
!
 Write(9,*) "# "
 Write(9,*) "# "
 Write(9,*) "# "
 Write(9,*) "# REGIONAL MAP - REGIONAL MAP - REGIONAL MAP - REGIONAL MAP - REGIONAL MAP# "
 Write(9,*) "# REGIONAL MAP - REGIONAL MAP - REGIONAL MAP - REGIONAL MAP - REGIONAL MAP# "
 Write(9,*) "# REGIONAL MAP - REGIONAL MAP - REGIONAL MAP - REGIONAL MAP - REGIONAL MAP# "
 Write(9,*) "# "
!
 Write(9,*)"gmtset PAPER_MEDIA A4+"
 Write(9,*)"gmtset LABEL_FONT_SIZE 24p"
 Write(9,*)"gmtset ANOT_FONT_SIZE 16p"
 Write(9,*)"gmtset FRAME_WIDTH 0.1c"
 Write(9,*)" "
!
  T_OPTION="-T-10/10/1"         ! Range of the palette 
  R_OPTION="-R-55/57/18/82r"    ! Range of plot (Louise)  
  J_OPTION="-JT-40/8c"
!
 DO 20 I=1, 3 
!
 if(i==1)namein ="smap_gree_reg.dat" 
 if(i==2)namein ="umap_gree_reg.dat" 
 if(i==3)namein ="nmap_gree_reg.dat" 
 if(i==1)nameout="smap_gree_reg.ps" 
 if(i==2)nameout="umap_gree_reg.ps" 
 if(i==3)nameout="nmap_gree_reg.ps" 
 if(i==1)namef="rrtmpfs.dat"
 if(i==2)namef="rrtmpfu.dat" 
 if(i==3)namef="rrtmpfn.dat"  
 if(i==1)TITRE="Rate of sea level change"
 if(i==2)TITRE="Vertical velocity"
 if(i==3)TITRE="Rate of geoid change"
!
 Write(9,*) " "
 if(i==1)Write(9,*)"# ---- Regional map of S at present time for ELASTIC REBOUND in GREENLAND ----" 		      	   
 if(i==2)Write(9,*)"# ---- Regional map of U at present time for ELASTIC REBOUND in GREENLAND ----" 		      	   
 if(i==3)Write(9,*)"# ---- Regional map of N at present time for ELASTIC REBOUND in GREENLAND ----" 		      	   
!
!
!
 if(i==1)Write(9,*) "echo", "    - S for Greenland - regional" 
 if(i==2)Write(9,*) "echo", "    - U for Greenland - regional" 
 if(i==3)Write(9,*) "echo", "    - N for Greenland - regional" 
!
!
 Write(9,*) " "
 Write(9,*) "makecpt -Cno_green ", trim(adjustl(T_OPTION)), " > pale.cpt" 
 Write(9,*) "psbasemap -X4 -Y4 -Ba10g10/a10g5WSen ",     trim(adjustl(R_OPTION)), " ", & 
             trim(adjustl(J_OPTION)), " -K > ", " ", trim(adjustl(nameout))  
 Write(9,*) "pscontour -I -O -K ",         & 
             trim(adjustl(R_OPTION)), " ", & 
   	     trim(adjustl(J_OPTION)), " ", &
             trim(adjustl(namein)), " -Cpale.cpt  >> ", & 
	     trim(adjustl(nameout)) 
!
 Write(9,*) "pscontour -O -K -G8 -W2/0/100/0 -A+g255 ",         & 
             trim(adjustl(R_OPTION)), " ", & 
   	     trim(adjustl(J_OPTION)), " ", &
             trim(adjustl(namein)), " -Cpale.cpt  >> ", & 
	     trim(adjustl(nameout)) 
!	     
 Write(9,*) "pscoast -Dh -B -W2/0 -A1000 -O -K ", & 
             trim(adjustl(R_OPTION)), " ",         & 
             trim(adjustl(J_OPTION)), " ", " >> ", & 
	     trim(adjustl(nameout))  
!
! A file for "titles" 
 open (4,file=NAMEF,status='unknown')  
!
! Alma 
 If(CODE=='-1')  THEN 
		Write(4,*) "307 +87     24  0 2 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "319 +54     16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	 	  "  -Elasticity: ALMA ", trim(adjustl(SHORT_VISCO))  
		Write(4,*) "319 +52.5   13  0 0 BC ",& 
 	   	 	   " -LMAX=", trim(adjustl(DEGREE)),& 
 	    		   " -RES=", trim(adjustl(RESOLUTION)),& 
	    		   " -NV=", trim(adjustl(NV)),& 
	    		   " -MODE=", trim(adjustl(MODE)),&   
	    		   " -ITER=", trim(adjustl(ITER))
		close(4) 
!
! Taboo 
 ELSEIF(CODE/='-1'.and.CODE/='-2')THEN 
 		Write(4,*) "307 +87     24  0 2 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "319 +54     16 0 0 BC   -Ice model: ", trim(adjustl(TITLICE)),& 
 	   	"  -Elasticity: TABOO "  
		Write(4,*) "319 +52.5   13  0 0 BC ",& 
 	    		   " -LMAX=", trim(adjustl(DEGREE)),& 
 	    		   " -RES=", trim(adjustl(RESOLUTION)),& 
	    	 	   " -NV=", trim(adjustl(NV)),& 
	    		   " -CODE=", trim(adjustl(CODE)),& 
	    		   " -MODE=", trim(adjustl(MODE)),&   
	    		   " -ITER=", trim(adjustl(ITER))
		close(4) 
!
 ELSEIF(CODE=='-2')THEN 
 		Write(4,*) "307 +87     24  0 2 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "319 +54     16 0 0 BC   -Ice model: ", trim(adjustl(TITLICE)),& 
 	   	"  -Elasticity: External "  
		Write(4,*) "319 +52.5   13  0 0 BC ",& 
 	    		   " -LMAX=", trim(adjustl(DEGREE)),& 
 	    		   " -RES=", trim(adjustl(RESOLUTION)),& 
	    		   " -CODE=", trim(adjustl(CODE)),& 
	    		   " -MODE=", trim(adjustl(MODE)),&   
	    		   " -ITER=", trim(adjustl(ITER))
		close(4) 
 ENDIF
!
		Write(9,*) "pstext -N -R -B ", trim(adjustl(NAMEF)), " ", &
		                               trim(adjustl(J_OPTION)), " -G0 -O -K >> ", & 
					       trim(adjustl(nameout)) 
		Write(9,*) "psscale -U/0.5/0.5/'SELEN 3.2' -E -Cpale.cpt -B1f2a2/:mm/yr: -D9/7.5/6/1 -O -K >> ", & 
	                    trim(adjustl(nameout)) 
 		Write(9,*) "psbasemap ", trim(adjustl(J_OPTION)), " ", & 
		           " -R -B::/::EWSN -O >> ", trim(adjustl(nameout)) 
!
!
20 CONTINUE 
!
!

!
!
 IF(OPTION_REB_GG=='y'.or.OPTION_REB_GR=='y') CLOSE(9)
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ---------- On Regional Greenland!
 ENDIF    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ---------- On Regional Greenland!
!
 END SUBROUTINE MAKE_ELAREB_MAPS_GREENLAND
!
!
!
!
!
!
   SUBROUTINE MAKE_TG_PLOTS_VAR  (TGFILE,         & 
                                  TG_FILE_FORMAT, &
				  FILE_GMT)
   IMPLICIT NONE
!
! ----- GS February 22 2012 
!       GS April 2014 
!
   CHARACTER*20  FILE_GMT
   CHARACTER*40  NAMEOUT
   CHARACTER*30  TGFILE 
   CHARACTER*1   TG_FILE_FORMAT
   CHARACTER*100 RIGA		    ! A row 
   INTEGER IJ, STAT_ID, TMIN, TMAX, NVAL
   REAL*8 STAT_LON, STAT_LAT, TREND, DTREND
   CHARACTER*30 STAT_NAME
   CHARACTER*15 INDIVIDUAL_FILENAME
   CHARACTER*18 INDIVIDUAL_FILENAME_PS    
   CHARACTER*15 FAKE_NAME
   CHARACTER*18 FAKE_NAME_PS  
   CHARACTER*1 CA(6) 
   INTEGER, PARAMETER :: MAXN=10000  ! A large number   
   INTEGER I, J, K, N, NH

   CHARACTER*80, PARAMETER :: B_OPTION="-Ba10f5:'time (yr)':/a0.2f0.1:'Sea level change (m)':SWen"
   CHARACTER*80, PARAMETER :: U_OPTION="-U'GS for ice2sea'"
   CHARACTER*80, PARAMETER :: R_OPTION="-R-5/205/-0.8/0.8"
   CHARACTER*80, PARAMETER :: J_OPTION="-JX17.33/7.4164"
   CHARACTER*80, PARAMETER :: XY_OPTION1="-X4 -Y4"
   CHARACTER*80, PARAMETER :: XY_OPTION2="-X0 -Y0"
   CHARACTER*80, PARAMETER :: W_OPTION="-W4"
   CHARACTER*80, PARAMETER :: W_OPTION2="-W4/-"
!
!
!
!
   OPEN( 9,file=file_gmt,status='unknown')
!
   Write(9,*)" "
   Write(9,*)"gmtset PAPER_MEDIA A4+"
   Write(9,*)"gmtset LABEL_FONT_SIZE 16p"
   Write(9,*)"gmtset ANOT_FONT_SIZE 12p"
   Write(9,*)"gmtset FRAME_WIDTH 0.01c"
   Write(9,*)"gmtset TICK_LENGTH -0.2c"
   Write(9,*)"gmtset FRAME_PEN = 0.5p"
   Write(9,*)" "
!
!
!
!
!
!
! --- Analysing the T/G database. 
!     Counting the header lines beginning with '#' and the data lines
!
!       write(*,*) tgfile
        open(10,file=TGFILE,status='unknown')
        nh=0 
        n=0
        do i=1, maxn 
       		read(10,'(a100)',end=1) riga
       		if(riga(1:1)=='#') nh=nh+1
       		if(riga(1:1)/='#') n=n+1
        enddo 
 1      close(10)  
!
! --- Opening the PSMSL data file... 
        open(10,file=TGFILE,status='unknown')
!
! --- Reading again the header lines  
        do i=1, nh 
	    read(10,'(a100)',end=1) riga 
      enddo  
!
!write(*,*) tgfile, N, TG_FILE_FORMAT
!
!
!
! ***********************
         DO 6 I=1, N       ! Loop on the tide-gauges 
! ***********************
!
! --- Reading one line      		     	    
         Read(10,119)  IJ, STAT_ID, STAT_LON, STAT_LAT, & 
	  	           TMIN,    TMAX,  NVAL,        & 
		           TREND,    DTREND,	        & 
		           STAT_NAME
!
! --- Name of individual files 
!
         DO K=1, 6 
	        IF(STAT_NAME(K:K)==' '.OR.STAT_NAME(K:K)=="/") then 
		     CA(K)="_"
		else
		     CA(K)=STAT_NAME(K:K) 
		endif	 
	 ENDDO
	 FAKE_NAME=   CA(1)//CA(2)//CA(3)//CA(4)//CA(5)//CA(6)//".tga"
	 FAKE_NAME_PS=CA(1)//CA(2)//CA(3)//CA(4)//CA(5)//CA(6)//".tga.ps"
!
         open(3,file='junk.dat',status='unknown') 
             if(              i.le.9)     write(3,'(A3,i1,A1,A10)') "000",i,"-",FAKE_NAME 
             if(10  .le.i.and.i.le.99)    write(3,'(A2,i2,A1,A10)')  "00",i,"-",FAKE_NAME 
             if(100 .le.i.and.i.le.999)   write(3,'(A1,i3,A1,A10)')   "0",i,"-",FAKE_NAME  
             if(1000.le.i.and.i.le.9999)  write(3,'(   i4,A1,A10)')       i,"-",FAKE_NAME
             if(              i.le.9)     write(3,'(A3,i1,A1,A13)') "000",i,"-",FAKE_NAME_PS
             if(10  .le.i.and.i.le.99)    write(3,'(A2,i2,A1,A13)')  "00",i,"-",FAKE_NAME_PS 
             if(100 .le.i.and.i.le.999)   write(3,'(A1,i3,A1,A13)')   "0",i,"-",FAKE_NAME_PS  
             if(1000.le.i.and.i.le.9999)  write(3,'(   i4,A1,A13)')       i,"-",FAKE_NAME_PS
         close(3)
	     open(3,file='junk.dat',status='unknown')
	     read(3,'(a15)') individual_filename
	     read(3,'(a18)') individual_filename_ps	    
	 close(3)    
!
!--- --- --- --- --- write(*,*) individual_filename, individual_filename_ps
!
         Write(9,*) " "
         IF(i==1.or.i==N.or.mod(i,100).eq.0) Write(9,*) "echo - Tide gauge plot [variations]", i, "/", n
         Write(9,*) "FILE_IN=", '"'//trim(adjustl(individual_filename))//'"'          
         Write(9,*) "FILE_OUT=", '"'//trim(adjustl(individual_filename_ps))//'"'       
         Write(9,*) " "
         Write(9,*) "sed '/#/ d' ", '$FILE_IN > d.tmp' 
         Write(9,*) "awk '{print $2, $3}' d.tmp > nd.tmp"
         Write(9,*) "psxy ", trim(adjustl(XY_OPTION1)), " ", trim(adjustl(R_OPTION)), " nd.tmp ", & 
	                     trim(adjustl(W_OPTION)),   " ", trim(adjustl(J_OPTION)), " "       , & 
			     trim(adjustl(B_OPTION)),   " -K >  $FILE_OUT" 
         Write(9,*) "awk '{print $2, $4}' d.tmp > nd.tmp"
         Write(9,*) "psxy ", trim(adjustl(XY_OPTION2)), " ", trim(adjustl(R_OPTION)), " nd.tmp ", & 
	                     trim(adjustl(W_OPTION2)),   " ", trim(adjustl(J_OPTION)), " "       , & 
			     trim(adjustl(B_OPTION)),   " -O -K >> $FILE_OUT" 
!
! --- T/G name label 
         Write(9,*) "pstext -N -G0 -W255 ", trim(adjustl(R_OPTION)), " ", trim(adjustl(J_OPTION)), " ", " -O -K <<END >> $FILE_OUT"
         Write(9,*) "0 0.4 12 0 2 ML ", trim(adjustl(individual_filename))
	 Write(9,'(A3)') "END"
!
! --- Label for "S" 
        Write(9,*) "psxy -B -R -JX -O ", trim(adjustl(W_OPTION)), " -K   <<END >>  $FILE_OUT"
        Write(9,*) "72  -0.3"
        Write(9,*) "85  -0.3"  
	Write(9,'(A3)') "END"
        Write(9,*) "echo '66 -0.3 11 0 0 LM S' | pstext  -N -R -JX -O -K >> $FILE_OUT" 
!
! --- Label for "N" 
        Write(9,*) "psxy -B -R -JX -O ", trim(adjustl(W_OPTION2)), " -K   <<END >>  $FILE_OUT"
        Write(9,*) "72  -0.4"
        Write(9,*) "85  -0.4"  
	Write(9,'(A3)') "END"
        Write(9,*) "echo '66 -0.4 11 0 0 LM N' | pstext  -N -R -JX -O -K >> $FILE_OUT" 
!
! --- Zero line 
        Write(9,*) "psxy -N -B -R -JX -O -W1/. <<END >>  $FILE_OUT" 
        Write(9,*) "0  0"
        Write(9,*) "100 0"  
	Write(9,'(A3)') "END"
!
!
! ***********************
 6  CONTINUE 
! ***********************
!
    CLOSE(9) 
!    
!
! --- Reading format for type '1'
119  FORMAT (2(i5, 1x), 2(f10.4, 1x), 3(i5, 1x), F7.2, 1x, F7.2, 4x, A30)
!
!
   END SUBROUTINE MAKE_TG_PLOTS_VAR 
!
!
!
!
!
!
 SUBROUTINE MAKE_GMAPS_ELA_MS     (DEPOT,       & 
                                   TITLICE,     & 
                                   RESOLUTION,  & 
				   DEGREE,      & 
				   NINC,        & 
				   CODE,        & 
				   MODE,        & 
				   ITER,        & 
				   file1_gmt,   & 		
				   file2_gmt,   & 
				   option_gia,  &		
                                   option_s,    & 
				   option_n,    & 
				   option_u,    & 			 
				   option_sdot, &
				   option_ndot, & 
				   option_udot)  
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Prepares a GMT script for *** GLOBAL MAPS *** maps of S, U & N for elastic multi-step ice models
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! GS February 22 2012 
! Depot in the title on April 12
! Depot in the title on April 12
! Option GIA on April 19 (Urbino Hospital...) 
!
 IMPLICIT NONE
 CHARACTER*12  DEPOT 
 CHARACTER*80  TITRE
 CHARACTER*30  NAMEIN, NAMEOUT, NAMEF   
 CHARACTER*10  TITLICE, RESOLUTION, DEGREE
 CHARACTER*20  FILE1_GMT, FILE2_GMT
 CHARACTER*20  R_OPTION, T_OPTION
 CHARACTER*3   NV, NINC, CODE
 CHARACTER*3   LABCHAR 
 CHARACTER*1   ITER, MODE 
 CHARACTER*1   OPTION_ROF
 CHARACTER*1   OPTION_S,    OPTION_N,     OPTION_U 
 CHARACTER*1   OPTION_SDOT, OPTION_NDOT,  OPTION_UDOT 
 CHARACTER*1  OPTION_GIA 
 INTEGER K, NN, KT 
!
!
 open( 9,file=file1_gmt,status='unknown')
 open(19,file=file2_gmt,status='unknown')
!
 Write(9,*)"gmtset PAPER_MEDIA A4+"
 Write(9,*)"gmtset LABEL_FONT_SIZE 24p"
 Write(9,*)"gmtset ANOT_FONT_SIZE 16p"
 Write(9,*)"gmtset FRAME_WIDTH 0.1c"

 Write(19,*)"gmtset PAPER_MEDIA A4+"
 Write(19,*)"gmtset LABEL_FONT_SIZE 24p"
 Write(19,*)"gmtset ANOT_FONT_SIZE 16p"
 Write(19,*)"gmtset FRAME_WIDTH 0.1c"

!
! ********************************************************************** SSSSSSSSSSSS
! ********************************************************************** SSSSSSSSSSSS
!
 IF (OPTION_S=='y') THEN
!
 Write(9,*)" "
 Write(9,*)"# ---- Global maps of S (VARIATIONS, NOT RATES) at all times  ----" 
 Write(9,*) "echo - Sea level variation S"		      		      
!
 R_OPTION="-R0/360/-80/80"      ! Range of the *Mercator* map  
 T_OPTION="-T-1.0/1.0/0.1"      ! Range of color table  

!
 open(3,file='junk.dat',status='unknown'); write(3,'(a3)') ninc ; close(3)
 open(3,file='junk.dat',status='unknown'); read (3,*)      nn   ; close(3)         
!
!
 DO 1 KT=0,NN+1  
open(3,file='junk.dat',status='unknown')
if(kt<=9)              write(3,'(a2,i1)') '00',kt  
if(kt>=10.and.kt <=99) write(3,'(a1,i2)')  '0',kt 
if(kt>=100)            write(3,   '(i3)')      kt ; close(3)
open(3,file='junk.dat',status='unknown'); read(3,'(a3)') labchar; close(3) 
!
!open(3,file='junk.dat',status='unknown') 
!if(kt<=9) write(3,'(a1,i1)') '0',kt  
!if(kt> 9) write(3,'(i2)')        kt ; close(3)
!open(3,file='junk.dat',status='unknown'); read(3,'(a2)') labchar; close(3)         
!
		NAMEIN ='svarmap-'//labchar//'.dat'
		NAMEOUT='svarmap-'//labchar//'.ps'
		namef="tmpfs-s"//labchar//".dat"
		open(101,file=NAMEIN,status='unknown')
	        TITRE="RSL variation (S), "//" time = "//trim(adjustl(labchar))
		Write(9,*) " "     
		Write(9,*) "makecpt -Cno_green ", trim(adjustl(T_OPTION)), " > pale.cpt" 
		Write(9,*) "psbasemap -X3 -Y5 -Ba180/a40WSEn -Jm0.018i ", & 
           	 trim(adjustl(R_OPTION)), " -K > ", " ", trim(adjustl(nameout))  
		Write(9,*) "pscontour -I -Jm -O -K ", trim(adjustl(R_OPTION)), " ", & 
            	            trim(adjustl(namein)), " -Cpale.cpt  >> ", trim(adjustl(nameout)) 
 	 	Write(9,*) "pscoast -Jm -Dc -B -W2/0 -G220  -A1000 -O -K ", trim(adjustl(R_OPTION)), " >> ", trim(adjustl(nameout))  
		open (4,file=NAMEF,status='unknown')  
		Write(4,*) "180 +84.2 26 0 3 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)), " ", trim(adjustl(DEPOT)) 
		Write(4,*) "180 -87.5   13  0 0 BC ",& 
 	    	" -LMAX=", trim(adjustl(DEGREE)),& 
 	    	" -RES=", trim(adjustl(RESOLUTION)),& 
	    	" -CODE=", trim(adjustl(CODE)),& 
	    	" -MODE=", trim(adjustl(MODE)),&   
	    	" -ITER=", trim(adjustl(ITER))
		close(4) 
		Write(9,*) "pstext -N -R -Jm -B ", trim(adjustl(NAMEF)), " -G0 -O -K >> ", trim(adjustl(nameout)) 		
        	Write(9,*) "psscale -U/0.5/0.5/'SELEN 3.2' -E -Cpale.cpt -B1f0.25a0.5/:m: -D8.25/-1/10/1h -O -K >> ", & 
                trim(adjustl(nameout))
 		Write(9,*) "psbasemap -J -R -B::/::EWSN -O >> ", trim(adjustl(nameout)) 
! 
1 CONTINUE 
ENDIF    
!
!
!
!
! ********************************************************************** NNNNNNNNNNNN
! ********************************************************************** NNNNNNNNNNNN
!
 IF (OPTION_N=='y') THEN
!
 Write(9,*)" "
 Write(9,*)"# ---- Global maps of N (VARIATIONS, NOT RATES) at all times  ----" 
 Write(9,*) "echo - Sea surface variation N"		      
!
 R_OPTION="-R0/360/-80/80"      ! Range of the *Mercator* map   
 T_OPTION="-T-1.0/1.0/0.1"      ! Range of color table  
!
 open(3,file='junk.dat',status='unknown'); write(3,'(a3)') ninc ; close(3)
 open(3,file='junk.dat',status='unknown'); read (3,*)      nn   ; close(3)         
!
!
 DO 2 KT=0,NN+1  
open(3,file='junk.dat',status='unknown')
if(kt<=9)              write(3,'(a2,i1)') '00',kt  
if(kt>=10.and.kt <=99) write(3,'(a1,i2)')  '0',kt 
if(kt>=100)            write(3,   '(i3)')      kt ; close(3)
open(3,file='junk.dat',status='unknown'); read(3,'(a3)') labchar; close(3) 
		NAMEIN ='nvarmap-'//labchar//'.dat'
		NAMEOUT='nvarmap-'//labchar//'.ps'
		namef="tmpfs-n"//labchar//".dat"
		open(101,file=NAMEIN,status='unknown')
	        TITRE="ASL variation (N), "//" time = "//trim(adjustl(labchar))
		Write(9,*) " "     
		Write(9,*) "makecpt -Cno_green ", trim(adjustl(T_OPTION)), " > pale.cpt" 
		Write(9,*) "psbasemap -X3 -Y5 -Ba180/a40WSEn -Jm0.018i ", & 
           	 trim(adjustl(R_OPTION)), " -K > ", " ", trim(adjustl(nameout))  
		Write(9,*) "pscontour -I -Jm -O -K ", trim(adjustl(R_OPTION)), " ", & 
            	            trim(adjustl(namein)), " -Cpale.cpt  >> ", trim(adjustl(nameout)) 
 	 	Write(9,*) "pscoast -Jm -Dc -B -W2/0 -G220  -A1000 -O -K ", trim(adjustl(R_OPTION)), " >> ", trim(adjustl(nameout))  
		open (4,file=NAMEF,status='unknown')  
		Write(4,*) "180 +84.2 26 0 3 BC ", trim(adjustl(TITRE)) 
!Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE))
		Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)), " ", trim(adjustl(DEPOT)) 
		Write(4,*) "180 -87.5   13  0 0 BC ",& 
 	    	" -LMAX=", trim(adjustl(DEGREE)),& 
 	    	" -RES=", trim(adjustl(RESOLUTION)),& 
	    	" -CODE=", trim(adjustl(CODE)),& 
	    	" -MODE=", trim(adjustl(MODE)),&   
	    	" -ITER=", trim(adjustl(ITER))
		close(4) 
		Write(9,*) "pstext -N -R -Jm -B ", trim(adjustl(NAMEF)), " -G0 -O -K >> ", trim(adjustl(nameout)) 		
        	Write(9,*) "psscale -U/0.5/0.5/'SELEN 3.2' -E -Cpale.cpt -B1f0.25a0.5/:m: -D8.25/-1/10/1h -O -K >> ", & 
                trim(adjustl(nameout))
 		Write(9,*) "psbasemap -J -R -B::/::EWSN -O >> ", trim(adjustl(nameout)) 
! 
2 CONTINUE 
ENDIF    
!
!
!
! ********************************************************************** UUUUUUUUUUUU
! ********************************************************************** UUUUUUUUUUUU
!
 IF (OPTION_U=='y') THEN
!
 Write(9,*)" "
 Write(9,*)"# ---- Global maps of U (VARIATIONS, NOT RATES) at all times  ----" 
 Write(9,*) "echo - Vertical displacement U"		      
!
 R_OPTION="-R0/360/-80/80"      ! Range of the *Mercator* map   
 T_OPTION="-T-1.0/1.0/0.1"      ! Range of color table  
!
 open(3,file='junk.dat',status='unknown'); write(3,'(a3)') ninc ; close(3)
 open(3,file='junk.dat',status='unknown'); read (3,*)      nn   ; close(3)         
!
!
 DO 26 KT=0,NN+1  
open(3,file='junk.dat',status='unknown')
if(kt<=9)              write(3,'(a2,i1)') '00',kt  
if(kt>=10.and.kt <=99) write(3,'(a1,i2)')  '0',kt 
if(kt>=100)            write(3,   '(i3)')      kt ; close(3)
open(3,file='junk.dat',status='unknown'); read(3,'(a3)') labchar; close(3) 
		NAMEIN ='uvarmap-'//labchar//'.dat'
		NAMEOUT='uvarmap-'//labchar//'.ps'
		namef="tmpfs-u"//labchar//".dat"
		open(101,file=NAMEIN,status='unknown')
	        TITRE="Vertical displacement (U), "//" time = "//trim(adjustl(labchar))
		Write(9,*) " "     
		Write(9,*) "makecpt -Cno_green ", trim(adjustl(T_OPTION)), " > pale.cpt" 
		Write(9,*) "psbasemap -X3 -Y5 -Ba180/a40WSEn -Jm0.018i ", & 
           	 trim(adjustl(R_OPTION)), " -K > ", " ", trim(adjustl(nameout))  
		Write(9,*) "pscontour -I -Jm -O -K ", trim(adjustl(R_OPTION)), " ", & 
            	            trim(adjustl(namein)), " -Cpale.cpt  >> ", trim(adjustl(nameout)) 
 	 	Write(9,*) "pscoast -Jm -Dc -B -W2/0  -G220 -A1000 -O -K ", trim(adjustl(R_OPTION)), " >> ", trim(adjustl(nameout))  
		open (4,file=NAMEF,status='unknown')  
		Write(4,*) "180 +84.2 26 0 3 BC ", trim(adjustl(TITRE)) 
!Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE))
		Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)), " ", trim(adjustl(DEPOT)) 
		Write(4,*) "180 -87.5   13  0 0 BC ",& 
 	    	" -LMAX=", trim(adjustl(DEGREE)),& 
 	    	" -RES=", trim(adjustl(RESOLUTION)),& 
	    	" -CODE=", trim(adjustl(CODE)),& 
	    	" -MODE=", trim(adjustl(MODE)),&   
	    	" -ITER=", trim(adjustl(ITER))
		close(4) 
		Write(9,*) "pstext -N -R -Jm -B ", trim(adjustl(NAMEF)), " -G0 -O -K >> ", trim(adjustl(nameout)) 		
        	Write(9,*) "psscale -U/0.5/0.5/'SELEN 3.2' -E -Cpale.cpt -B1f0.25a0.5/:m: -D8.25/-1/10/1h -O -K >> ", & 
                trim(adjustl(nameout))
 		Write(9,*) "psbasemap -J -R -B::/::EWSN -O >> ", trim(adjustl(nameout)) 
! 
26 CONTINUE 
ENDIF    
!
!
!
! ********************************************************************** SSSSSSSSSSSS DOT-DOT-DOT-DOT-DOT
! ********************************************************************** SSSSSSSSSSSS DOT-DOT-DOT-DOT-DOT
!
 IF (OPTION_SDOT=='y') THEN
!
 Write(19,*)" "
 Write(19,*)"# ---- Global maps of S-dot (RATE, not VARIATION) at all times  ----" 
 Write(19,*) "echo - Rate of Sea level variation dot-S"		      		      
!
 R_OPTION="-R0/360/-80/80"      ! Range of the *Mercator* map common to all maps 
 T_OPTION="-T-5.0/5.0/0.5"      ! Range of color table for ALL plots 
!
 open(3,file='junk.dat',status='unknown'); write(3,'(a3)') ninc ; close(3)
 open(3,file='junk.dat',status='unknown'); read (3,*)      nn   ; close(3)         
!
!
 DO 3 KT=1,NN-1 ! Note: from "1" to "N-1"  
open(3,file='junk.dat',status='unknown')
if(kt<=9)              write(3,'(a2,i1)') '00',kt  
if(kt>=10.and.kt <=99) write(3,'(a1,i2)')  '0',kt 
if(kt>=100)            write(3,   '(i3)')      kt ; close(3)
open(3,file='junk.dat',status='unknown'); read(3,'(a3)') labchar; close(3) 
		NAMEIN ='sdotmap-'//labchar//'.dat'
		NAMEOUT='sdotmap-'//labchar//'.ps'
		namef="tmpfs-sd"//labchar//".dat"
		open(101,file=NAMEIN,status='unknown')
	        TITRE="Rate of RSL change"//" time = "//trim(adjustl(labchar))
		Write(19,*) " "     
		Write(19,*) "makecpt -Cno_green ", trim(adjustl(T_OPTION)), " > pale.cpt" 
		Write(19,*) "psbasemap -X3 -Y5 -Ba180/a40WSEn -Jm0.018i ", & 
           	 trim(adjustl(R_OPTION)), " -K > ", " ", trim(adjustl(nameout))  
		Write(19,*) "pscontour -I -Jm -O -K ", trim(adjustl(R_OPTION)), " ", & 
            	            trim(adjustl(namein)), " -Cpale.cpt  >> ", trim(adjustl(nameout)) 
 	 	Write(19,*) "pscoast -Jm -Dc -B -W2/0 -G220  -A1000 -O -K ", trim(adjustl(R_OPTION)), " >> ", trim(adjustl(nameout))  
		open (4,file=NAMEF,status='unknown')  
		Write(4,*) "180 +84.2 26 0 3 BC ", trim(adjustl(TITRE)) 
!Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE))
		Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)), " ", trim(adjustl(DEPOT)) 
		Write(4,*) "180 -87.5   13  0 0 BC ",& 
 	    	" -LMAX=", trim(adjustl(DEGREE)),& 
 	    	" -RES=", trim(adjustl(RESOLUTION)),& 
	    	" -CODE=", trim(adjustl(CODE)),& 
	    	" -MODE=", trim(adjustl(MODE)),&   
	    	" -ITER=", trim(adjustl(ITER))
		close(4) 
		Write(19,*) "pstext -N -R -Jm -B ", trim(adjustl(NAMEF)), " -G0 -O -K >> ", trim(adjustl(nameout)) 		
        	Write(19,*) "psscale -U/0.5/0.5/'SELEN 3.2' -E -Cpale.cpt -B1f0.5a1.0/:mm/yr: -D8.25/-1/10/1h -O -K >> ", & 
                trim(adjustl(nameout))
 		Write(19,*) "psbasemap -J -R -B::/::EWSN -O >> ", trim(adjustl(nameout)) 
! 
3 CONTINUE 
ENDIF    
!
!
!
!
!
! ********************************************************************** NNNNNNNNNNNN DOT-DOT-DOT-DOT-DOT
! ********************************************************************** NNNNNNNNNNNN DOT-DOT-DOT-DOT-DOT
!
 IF (OPTION_NDOT=='y') THEN
!
 Write(19,*)" "
 Write(19,*)"# ---- Global maps of S-dot (RATE, not VARIATION) at all times  ----" 
 Write(19,*) "echo - Rate of Sea surface variation dot-N"		      		      
!
 R_OPTION="-R0/360/-80/80"      ! Range of the *Mercator* map common to all maps 
 T_OPTION="-T-5.0/5.0/0.5"      ! Range of color table for ALL plots 
!
 open(3,file='junk.dat',status='unknown'); write(3,'(a3)') ninc ; close(3)
 open(3,file='junk.dat',status='unknown'); read (3,*)      nn   ; close(3)         
!
!
 DO 4 KT=1,NN-1  ! Note: from "1" to "N-1"  
open(3,file='junk.dat',status='unknown')
if(kt<=9)              write(3,'(a2,i1)') '00',kt  
if(kt>=10.and.kt <=99) write(3,'(a1,i2)')  '0',kt 
if(kt>=100)            write(3,   '(i3)')      kt ; close(3)
open(3,file='junk.dat',status='unknown'); read(3,'(a3)') labchar; close(3) 
		NAMEIN ='ndotmap-'//labchar//'.dat'
		NAMEOUT='ndotmap-'//labchar//'.ps'
		namef="tmpfs-nd"//labchar//".dat"
		open(101,file=NAMEIN,status='unknown')
	        TITRE="Rate of ASL variation"//" time = "//trim(adjustl(labchar))
		Write(19,*) " "     
		Write(19,*) "makecpt -Cno_green ", trim(adjustl(T_OPTION)), " > pale.cpt" 
		Write(19,*) "psbasemap -X3 -Y5 -Ba180/a40WSEn -Jm0.018i ", & 
           	 trim(adjustl(R_OPTION)), " -K > ", " ", trim(adjustl(nameout))  
		Write(19,*) "pscontour -I -Jm -O -K ", trim(adjustl(R_OPTION)), " ", & 
            	            trim(adjustl(namein)), " -Cpale.cpt  >> ", trim(adjustl(nameout)) 
 	 	Write(19,*) "pscoast -Jm -Dc -B -W2/0 -G220  -A1000 -O -K ", trim(adjustl(R_OPTION)), " >> ", trim(adjustl(nameout))  
		open (4,file=NAMEF,status='unknown')  
		Write(4,*) "180 +84.2 26 0 3 BC ", trim(adjustl(TITRE)) 
!Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE))
		Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)), " ", trim(adjustl(DEPOT)) 
		Write(4,*) "180 -87.5   13  0 0 BC ",& 
 	    	" -LMAX=", trim(adjustl(DEGREE)),& 
 	    	" -RES=", trim(adjustl(RESOLUTION)),& 
	    	" -CODE=", trim(adjustl(CODE)),& 
	    	" -MODE=", trim(adjustl(MODE)),&   
	    	" -ITER=", trim(adjustl(ITER))
		close(4) 
		Write(19,*) "pstext -N -R -Jm -B ", trim(adjustl(NAMEF)), " -G0 -O -K >> ", trim(adjustl(nameout)) 		
        	Write(19,*) "psscale -U/0.5/0.5/'SELEN 3.2' -E -Cpale.cpt -B1f0.5a1.0/:mm/yr: -D8.25/-1/10/1h -O -K >> ", & 
                trim(adjustl(nameout))
 		Write(19,*) "psbasemap -J -R -B::/::EWSN -O >> ", trim(adjustl(nameout)) 
! 
4 CONTINUE 
ENDIF    
!
!
!
! ********************************************************************** UUUUUUUUUUUU DOT-DOT-DOT-DOT-DOT
! ********************************************************************** UUUUUUUUUUUU DOT-DOT-DOT-DOT-DOT
!
 IF (OPTION_UDOT=='y') THEN
!
 Write(19,*)" "
 Write(19,*)"# ---- Global maps of U-dot (RATE, not VARIATION) at all times  ----" 
 Write(19,*) "echo - Vertical velocity dot-U"		      		      
!
 R_OPTION="-R0/360/-80/80"      ! Range of the *Mercator* map common to all maps 
 T_OPTION="-T-5.0/5.0/0.5"      ! Range of color table for ALL plots 
!
 open(3,file='junk.dat',status='unknown'); write(3,'(a3)') ninc ; close(3)
 open(3,file='junk.dat',status='unknown'); read (3,*)      nn   ; close(3)         
!
!
 DO 46 KT=1,NN-1  ! Note: from "1" to "N-1"  
open(3,file='junk.dat',status='unknown')
if(kt<=9)              write(3,'(a2,i1)') '00',kt  
if(kt>=10.and.kt <=99) write(3,'(a1,i2)')  '0',kt 
if(kt>=100)            write(3,   '(i3)')      kt ; close(3)
open(3,file='junk.dat',status='unknown'); read(3,'(a3)') labchar; close(3) 
		NAMEIN ='udotmap-'//labchar//'.dat'
		NAMEOUT='udotmap-'//labchar//'.ps'
		namef="tmpfs-ud"//labchar//".dat"
		open(101,file=NAMEIN,status='unknown')
	        TITRE="Vertical velocity"//" time = "//trim(adjustl(labchar))
		Write(19,*) " "     
		Write(19,*) "makecpt -Cno_green ", trim(adjustl(T_OPTION)), " > pale.cpt" 
		Write(19,*) "psbasemap -X3 -Y5 -Ba180/a40WSEn -Jm0.018i ", & 
           	 trim(adjustl(R_OPTION)), " -K > ", " ", trim(adjustl(nameout))  
		Write(19,*) "pscontour -I -Jm -O -K ", trim(adjustl(R_OPTION)), " ", & 
            	            trim(adjustl(namein)), " -Cpale.cpt  >> ", trim(adjustl(nameout)) 
 	 	Write(19,*) "pscoast -Jm -Dc -B -W2/0  -G220 -A1000 -O -K ", trim(adjustl(R_OPTION)), " >> ", trim(adjustl(nameout))  
		open (4,file=NAMEF,status='unknown')  
		Write(4,*) "180 +84.2 26 0 3 BC ", trim(adjustl(TITRE)) 
!Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE))
		Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)), " ", trim(adjustl(DEPOT)) 
		Write(4,*) "180 -87.5   13  0 0 BC ",& 
 	    	" -LMAX=", trim(adjustl(DEGREE)),& 
 	    	" -RES=", trim(adjustl(RESOLUTION)),& 
	    	" -CODE=", trim(adjustl(CODE)),& 
	    	" -MODE=", trim(adjustl(MODE)),&   
	    	" -ITER=", trim(adjustl(ITER))
		close(4) 
		Write(19,*) "pstext -N -R -Jm -B ", trim(adjustl(NAMEF)), " -G0 -O -K >> ", trim(adjustl(nameout)) 		
        	Write(19,*) "psscale -U/0.5/0.5/'SELEN 3.2' -E -Cpale.cpt -B1f0.5a1.0/:mm/yr: -D8.25/-1/10/1h -O -K >> ", & 
                trim(adjustl(nameout))
 		Write(19,*) "psbasemap -J -R -B::/::EWSN -O >> ", trim(adjustl(nameout)) 
! 
46 CONTINUE 
ENDIF    
!
!
!
 CLOSE(9)
 CLOSE(19)

!
!
 END SUBROUTINE MAKE_GMAPS_ELA_MS




!
!
!
 SUBROUTINE MAKE_GMAPS (TITLICE,    & 
 			RESOLUTION, &
			NV,   & 
			CODE, & 
			ITER, & 
 	                MODE, & 
			DEGREE,     & 
			VSTRING,    & 
			OPTION_ROF, & 
			FILECAP,    & 
			FILE1_GMT,  &
			SHORT_VISCO,   & 
                        OPTION_DOTS,   & 
			OPTION_DOTU,   & 
			OPTION_DOTN,   & 
			OPTION_DOTG,   & 
			OPTION_DOTW,   & 
			OPTION_DOTFA,  & 
			OPTION_DOTSS,  &			
			OPTION_DOTLOI, & 
			OPTION_DOTLOO, & 
			OPTION_DOTLOT) 
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Prepares a GMT script for *** GLOBAL MAPS *** maps of dot-S, U & N at present time 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Last change GS March 24 08 
! Modified on June 24, 2008 for SELEN 2.6 GS 
! Modified on July 26, 2008 for SELEN 2.6 GS 
! Upgraded April 2010 by GS for ALMA on g95 
! Upgraded June 2010 by GS for Free air gravity   
! Upgraded Marc 2011 by GS for \Phi \over \gamma 
! Modified GS May 28, 2011. A title (for the geoid) was wrong... 
! Modified GS April 15 2012 Munich airport for the ice2sea work
! Aug 5, 2012; Equivalent water height (understanding the Chambers et al. 2010 results)
!
 IMPLICIT NONE
 CHARACTER*100 SHORT_VISCO 
 CHARACTER*80  TITRE
 CHARACTER*30  NAMEIN, NAMEOUT, NAMEF   
 CHARACTER*10  TITLICE, RESOLUTION, DEGREE
 CHARACTER*30  VSTRING
 CHARACTER*20  FILECAP, FILE1_GMT
 CHARACTER*20  R_OPTION, T_OPTION
 CHARACTER*3   NV, CODE
 CHARACTER*1   OPTION_ROF
 CHARACTER*1   ITER, MODE 
 CHARACTER*1   OPTION_DOTS,     & 
               OPTION_DOTU,     & 
	       OPTION_DOTN,     & 
	       OPTION_DOTG,     & 
	       OPTION_DOTW,     &     ! News as August 2012 XXX
	       OPTION_DOTFA,    & 
	       OPTION_DOTSS,    & 
	       OPTION_DOTLOI,   & 
	       OPTION_DOTLOO,   & 
	       OPTION_DOTLOT
 INTEGER K
!
!
!
 open(9,file=file1_gmt,status='unknown')
 Write(9,*)"gmtset PAPER_MEDIA A4+"
 Write(9,*)"gmtset LABEL_FONT_SIZE 24p"
 Write(9,*)"gmtset ANOT_FONT_SIZE 16p"
 Write(9,*)"gmtset FRAME_WIDTH 0.1c"

 Write(9,*)" "
 Write(9,*)"# ---- Global maps of dot S, U, and N at present time  ----" 		      
!
 R_OPTION="-R0/360/-80/80"      ! Range of the *Mercator* map common to all maps 
!
!
!
 	if    (OPTION_DOTS=='y') THEN 
	        Write(9,*) ""
	        namein ="sdotmap.dat"; nameout="sdotmap.ps"; namef="tmpfs.dat"
	        TITRE="Rate of RSL change GIA"
		T_OPTION="-T-1.0/1.0/0.1"         ! Range of the palette 
		Write(9,*) " "
		Write(9,*) "makecpt -CGMT_panoply ", trim(adjustl(T_OPTION)), " > pale.cpt" 
		Write(9,*) "psbasemap -X3 -Y5 -Ba180/a40WSEn -Jm0.018i ", & 
           	 trim(adjustl(R_OPTION)), " -K > ", " ", trim(adjustl(nameout))  
		Write(9,*) "pscontour -I -Jm -O -K ", trim(adjustl(R_OPTION)), " ", & 
            	trim(adjustl(namein)), " -Cpale.cpt  >> ", trim(adjustl(nameout)) 
		If(option_rof=='r')then 
 	 	Write(9,*) "pscoast -Jm -Dc -B -W2/0  -G220 -A1000 -O -K ", trim(adjustl(R_OPTION)), " >> ", trim(adjustl(nameout))  
  		else
  		Write(9,*) "pscoast -Jm -Dc -B -W1/240 -G220 -A1000 -O -K ", trim(adjustl(R_OPTION)), " >> ", trim(adjustl(nameout))  
		Endif
		If(option_rof=='z')& 
		Write(9,*) "psxy -R -Jm ", trim(adjustl(filecap)), " -M -W4/0 -B -A -O -K >> ", trim(adjustl(nameout))  
		open (4,file=NAMEF,status='unknown')  
		If(CODE=='-1')    THEN 
		Write(4,*) "180 +84.2 26 0 3 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	 	  "  -ALMA rheology: ", trim(adjustl(SHORT_VISCO))  
		Write(4,*) "180 -87.5   13  0 0 BC ",& 
 	   	 " -LMAX=", trim(adjustl(DEGREE)),& 
 	    	" -RES=", trim(adjustl(RESOLUTION)),& 
	    	" -NV=", trim(adjustl(NV)),& 
	    	" -MODE=", trim(adjustl(MODE)),&   
	    	" -ITER=", trim(adjustl(ITER));  close(4) 
		Write(9,*) "pstext -N -R -Jm -B ", trim(adjustl(NAMEF)), " -G0 -O -K >> ", trim(adjustl(nameout)) 
		Write(9,*) "psscale -U/0.5/0.5/'SELEN 3.2' -E -Cpale.cpt -Bf0.1a0.5/:mm/yr: -D8.25/-1/10/1h -O >> ", & 
	            trim(adjustl(nameout)) 
		ELSEIF(CODE/='-1')THEN 
		Write(4,*) "180 +84.2 26 0 3 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	   	"  -Viscosity profile: ", trim(adjustl(VSTRING))  
		Write(4,*) "180 -87.5   13  0 0 BC ",& 
 	    	" -LMAX=", trim(adjustl(DEGREE)),& 
 	    	" -RES=", trim(adjustl(RESOLUTION)),& 
	    	" -NV=", trim(adjustl(NV)),& 
	    	" -CODE=", trim(adjustl(CODE)),& 
	    	" -MODE=", trim(adjustl(MODE)),&   
	    	" -ITER=", trim(adjustl(ITER));  close(4) 
		Write(9,*) "pstext -N -R -Jm -B ", trim(adjustl(NAMEF)), " -G0 -O -K >> ", trim(adjustl(nameout)) 		
        	Write(9,*) "psscale -U/0.5/0.5/'SELEN 3.2' -E -Cpale.cpt -Bf0.1a0.5/:mm/yr: -D8.25/-1/10/1h -O -K >> ", & 
                trim(adjustl(nameout))
 		Write(9,*) "psbasemap -J -R -B::/::EWSN -O >> ", trim(adjustl(nameout)) 
		ENDIF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	ENDIF    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 	if    (OPTION_DOTU=='y') THEN 
	        Write(9,*) ""
	        namein ="udotmap.dat"; nameout="udotmap.ps"; namef="tmpfu.dat"
	        TITRE="Rate of vertical uplift GIA"
		T_OPTION="-T-1.0/1.0/0.1"         ! Range of the palette 
		Write(9,*) " "
		Write(9,*) "makecpt -CGMT_panoply ", trim(adjustl(T_OPTION)), " > pale.cpt" 
		Write(9,*) "psbasemap -X3 -Y5 -Ba180/a40WSEn -Jm0.018i ", & 
           	 trim(adjustl(R_OPTION)), " -K > ", " ", trim(adjustl(nameout))  
		Write(9,*) "pscontour -I -Jm -O -K ", trim(adjustl(R_OPTION)), " ", & 
            	trim(adjustl(namein)), " -Cpale.cpt  >> ", trim(adjustl(nameout)) 
		If(option_rof=='r')then 
 	 	Write(9,*) "pscoast -Jm -Dc -B -W2/0   -A1000 -O -K ", trim(adjustl(R_OPTION)), " >> ", trim(adjustl(nameout))  
  		else
  		Write(9,*) "pscoast -Jm -Dc -B -W1/240 -A1000 -O -K ", trim(adjustl(R_OPTION)), " >> ", trim(adjustl(nameout))  
		Endif
		If(option_rof=='z')& 
		Write(9,*) "psxy -R -Jm ", trim(adjustl(filecap)), " -M -W4/0 -B -A -O -K >> ", trim(adjustl(nameout))  
		open (4,file=NAMEF,status='unknown')  
		If(CODE=='-1')    THEN 
		Write(4,*) "180 +84.2 26 0 3 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	 	  "  -ALMA rheology: ", trim(adjustl(SHORT_VISCO))  
		Write(4,*) "180 -87.5   13  0 0 BC ",& 
 	   	 " -LMAX=", trim(adjustl(DEGREE)),& 
 	    	" -RES=", trim(adjustl(RESOLUTION)),& 
	    	" -NV=", trim(adjustl(NV)),& 
	    	" -MODE=", trim(adjustl(MODE)),&   
	    	" -ITER=", trim(adjustl(ITER));  close(4) 
		Write(9,*) "pstext -N -R -Jm -B ", trim(adjustl(NAMEF)), " -G0 -O -K >> ", trim(adjustl(nameout)) 
		Write(9,*) "psscale -U/0.5/0.5/'SELEN 3.2' -E -Cpale.cpt -Bf0.1a0.5/:mm/yr: -D8.25/-1/10/1h -O >> ", & 
	            trim(adjustl(nameout)) 
		ELSEIF(CODE/='-1')THEN 
		Write(4,*) "180 +84.2 26 0 3 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	   	"  -Viscosity profile: ", trim(adjustl(VSTRING))  
		Write(4,*) "180 -87.5   13  0 0 BC ",& 
 	    	" -LMAX=", trim(adjustl(DEGREE)),& 
 	    	" -RES=", trim(adjustl(RESOLUTION)),& 
	    	" -NV=", trim(adjustl(NV)),& 
	    	" -CODE=", trim(adjustl(CODE)),& 
	    	" -MODE=", trim(adjustl(MODE)),&   
	    	" -ITER=", trim(adjustl(ITER));  close(4) 
		Write(9,*) "pstext -N -R -Jm -B ", trim(adjustl(NAMEF)), " -G0 -O -K >> ", trim(adjustl(nameout)) 		
        	Write(9,*) "psscale -U/0.5/0.5/'SELEN 3.2' -E -Cpale.cpt -Bf0.1a0.5/:mm/yr: -D8.25/-1/10/1h -O -K >> ", & 
                trim(adjustl(nameout))
 		Write(9,*) "psbasemap -J -R -B::/::EWSN -O >> ", trim(adjustl(nameout)) 
		ENDIF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	ENDIF     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 	if    (OPTION_DOTN=='y') THEN 
	        Write(9,*) ""
	        namein ="ndotmap.dat"; nameout="ndotmap.ps"; namef="tmpfn.dat"
	        TITRE="Rate of ASL change GIA"
		T_OPTION="-T-1.0/1.0/0.1"         ! Range of the palette 
		Write(9,*) " "
		Write(9,*) "makecpt -CGMT_panoply ", trim(adjustl(T_OPTION)), " > pale.cpt" 
		Write(9,*) "psbasemap -X3 -Y5 -Ba180/a40WSEn -Jm0.018i ", & 
           	 trim(adjustl(R_OPTION)), " -K > ", " ", trim(adjustl(nameout))  
		Write(9,*) "pscontour -I -Jm -O -K ", trim(adjustl(R_OPTION)), " ", & 
            	trim(adjustl(namein)), " -Cpale.cpt  >> ", trim(adjustl(nameout)) 
		If(option_rof=='r')then 
 	 	Write(9,*) "pscoast -Jm -Dc -B -W2/0   -A1000 -O -K ", trim(adjustl(R_OPTION)), " >> ", trim(adjustl(nameout))  
  		else
  		Write(9,*) "pscoast -Jm -Dc -B -W1/240 -A1000 -O -K ", trim(adjustl(R_OPTION)), " >> ", trim(adjustl(nameout))  
		Endif
		If(option_rof=='z')& 
		Write(9,*) "psxy -R -Jm ", trim(adjustl(filecap)), " -M -W4/0 -B -A -O -K >> ", trim(adjustl(nameout))  
		open (4,file=NAMEF,status='unknown')  
		If(CODE=='-1')    THEN 
		Write(4,*) "180 +84.2 26 0 3 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	 	  "  -ALMA rheology: ", trim(adjustl(SHORT_VISCO))  
		Write(4,*) "180 -87.5   13  0 0 BC ",& 
 	   	 " -LMAX=", trim(adjustl(DEGREE)),& 
 	    	" -RES=", trim(adjustl(RESOLUTION)),& 
	    	" -NV=", trim(adjustl(NV)),& 
	    	" -MODE=", trim(adjustl(MODE)),&   
	    	" -ITER=", trim(adjustl(ITER));  close(4) 
		Write(9,*) "pstext -N -R -Jm -B ", trim(adjustl(NAMEF)), " -G0 -O -K >> ", trim(adjustl(nameout)) 
		Write(9,*) "psscale -U/0.5/0.5/'SELEN 3.2' -E -Cpale.cpt -Bf0.1a0.5/:mm/yr: -D8.25/-1/10/1h -O >> ", & 
	            trim(adjustl(nameout)) 
		ELSEIF(CODE/='-1')THEN 
		Write(4,*) "180 +84.2 26 0 3 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	   	"  -Viscosity profile: ", trim(adjustl(VSTRING))  
		Write(4,*) "180 -87.5   13  0 0 BC ",& 
 	    	" -LMAX=", trim(adjustl(DEGREE)),& 
 	    	" -RES=", trim(adjustl(RESOLUTION)),& 
	    	" -NV=", trim(adjustl(NV)),& 
	    	" -CODE=", trim(adjustl(CODE)),& 
	    	" -MODE=", trim(adjustl(MODE)),&   
	    	" -ITER=", trim(adjustl(ITER));  close(4) 
		Write(9,*) "pstext -N -R -Jm -B ", trim(adjustl(NAMEF)), " -G0 -O -K >> ", trim(adjustl(nameout)) 		
        	Write(9,*) "psscale -U/0.5/0.5/'SELEN 3.2' -E -Cpale.cpt -Bf0.1a0.5/:mm/yr: -D8.25/-1/10/1h -O -K >> ", & 
                trim(adjustl(nameout))
 		Write(9,*) "psbasemap -J -R -B::/::EWSN -O >> ", trim(adjustl(nameout)) 
		ENDIF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	ENDIF    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 	if    (OPTION_DOTG=='y') THEN 
	        Write(9,*) ""
	        namein ="gdotmap.dat"; nameout="gdotmap.ps"; namef="tmpfg.dat"
	        TITRE="Rate of Phi/gamma today"
		T_OPTION="-T-1.0/1.0/0.1"         ! Range of the palette 
		Write(9,*) " "
		Write(9,*) "makecpt -CGMT_panoply ", trim(adjustl(T_OPTION)), " > pale.cpt" 
		Write(9,*) "psbasemap -X3 -Y5 -Ba180/a40WSEn -Jm0.018i ", & 
           	 trim(adjustl(R_OPTION)), " -K > ", " ", trim(adjustl(nameout))  
		Write(9,*) "pscontour -I -Jm -O -K ", trim(adjustl(R_OPTION)), " ", & 
            	trim(adjustl(namein)), " -Cpale.cpt  >> ", trim(adjustl(nameout)) 
		If(option_rof=='r')then 
 	 	Write(9,*) "pscoast -Jm -Dc -B -W2/0   -A1000 -O -K ", trim(adjustl(R_OPTION)), " >> ", trim(adjustl(nameout))  
  		else
  		Write(9,*) "pscoast -Jm -Dc -B -W1/240 -A1000 -O -K ", trim(adjustl(R_OPTION)), " >> ", trim(adjustl(nameout))  
		Endif
		If(option_rof=='z')& 
		Write(9,*) "psxy -R -Jm ", trim(adjustl(filecap)), " -M -W4/0 -B -A -O -K >> ", trim(adjustl(nameout))  
		open (4,file=NAMEF,status='unknown')  
		If(CODE=='-1')    THEN 
		Write(4,*) "180 +84.2 26 0 3 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	 	  "  -ALMA rheology: ", trim(adjustl(SHORT_VISCO))  
		Write(4,*) "180 -87.5   13  0 0 BC ",& 
 	   	 " -LMAX=", trim(adjustl(DEGREE)),& 
 	    	" -RES=", trim(adjustl(RESOLUTION)),& 
	    	" -NV=", trim(adjustl(NV)),& 
	    	" -MODE=", trim(adjustl(MODE)),&   
	    	" -ITER=", trim(adjustl(ITER));  close(4) 
		Write(9,*) "pstext -N -R -Jm -B ", trim(adjustl(NAMEF)), " -G0 -O -K >> ", trim(adjustl(nameout)) 
		Write(9,*) "psscale -U/0.5/0.5/'SELEN 3.2' -E -Cpale.cpt -Bf0.1a0.5/:mm/yr: -D8.25/-1/10/1h -O >> ", & 
	            trim(adjustl(nameout)) 
		ELSEIF(CODE/='-1')THEN 
		Write(4,*) "180 +84.2 26 0 3 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	   	"  -Viscosity profile: ", trim(adjustl(VSTRING))  
		Write(4,*) "180 -87.5   13  0 0 BC ",& 
 	    	" -LMAX=", trim(adjustl(DEGREE)),& 
 	    	" -RES=", trim(adjustl(RESOLUTION)),& 
	    	" -NV=", trim(adjustl(NV)),& 
	    	" -CODE=", trim(adjustl(CODE)),& 
	    	" -MODE=", trim(adjustl(MODE)),&   
	    	" -ITER=", trim(adjustl(ITER));  close(4) 
		Write(9,*) "pstext -N -R -Jm -B ", trim(adjustl(NAMEF)), " -G0 -O -K >> ", trim(adjustl(nameout)) 		
        	Write(9,*) "psscale -U/0.5/0.5/'SELEN 3.2' -E -Cpale.cpt -Bf0.1a0.5/:mm/yr: -D8.25/-1/10/1h -O -K >> ", & 
                trim(adjustl(nameout))
 		Write(9,*) "psbasemap -J -R -B::/::EWSN -O >> ", trim(adjustl(nameout)) 
		ENDIF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	ENDIF    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!
!
! In progress: Equivalent water height 
!
!
!
 	if    (OPTION_DOTW=='y') THEN 

                R_OPTION="-R0/360/-90/90"      ! Range of the *Mollweide* map 
!
	        Write(9,*) ""
                Write(9,*) "gmtset HEADER_FONT_SIZE 32p"
!
	        namein ="wdotmap.dat"; nameout="wdotmap.ps"; namef="tmpfw.dat"
!
		T_OPTION="-T-4.0/4.0/0.1"         ! Range of the palette 
!
		Write(9,*) "makecpt -Cno_green ", trim(adjustl(T_OPTION)), " > pale.cpt" 
!
		Write(9,*) 'psbasemap -X3 -Y5 -Ba180/a45WSEn:."Equivalent water height (GIA)": -JW0/16 ', & 
		trim(adjustl(R_OPTION)), " -K > ", " ", trim(adjustl(nameout))  
!
		Write(9,*) "pscontour -I -JW -O -K ", trim(adjustl(R_OPTION)), " ", & 
            	trim(adjustl(namein)), " -Cpale.cpt  >> ", trim(adjustl(nameout)) 
!
 	 	Write(9,*) "pscoast -JW -Dc -B -W2/0   -A1000 -O -K ", trim(adjustl(R_OPTION)), " >> ", trim(adjustl(nameout))  
!
		open (4,file=NAMEF,status='unknown')  
!
	        IF(CODE/='-1')THEN 
		Write(4,*) "180 -170 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)), & 
 	   	"  -Viscosity profile: ", trim(adjustl(VSTRING))  
		Write(4,*) "180 -185 12  0 0 BC ",& 
 	    	" -LMAX=", trim(adjustl(DEGREE)),& 
 	    	" -RES=", trim(adjustl(RESOLUTION)),& 
	    	" -NV=", trim(adjustl(NV)),& 
	    	" -CODE=", trim(adjustl(CODE)),& 
	    	" -MODE=", trim(adjustl(MODE)),&   
	    	" -ITER=", trim(adjustl(ITER));  close(4) 
        	Write(9,*) "psscale -U/0/-0.8/'SELEN 3.2' -E -Cpale.cpt -Bf0.5a1.0/:mm/yr: -D8.25/-1/10/1h -O -K >> ", & 
                trim(adjustl(nameout))
		Write(9,*) "pstext -N -R -JX16/8 ", trim(adjustl(NAMEF)), " -G0 -O >> ", trim(adjustl(nameout)) 		
                Write(9,*) "gmtset HEADER_FONT_SIZE 36p"
!
!               Write(9,*) "psbasemap -J -R -B::/::EWSN -O >> ", trim(adjustl(nameout)) 
!
		ENDIF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	ENDIF    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!









 	if    (OPTION_DOTFA=='y') THEN 
	        Write(9,*) ""
	        namein ="fadotmap.dat"; nameout="fadotmap.ps"; namef="tmpffa.dat"
	        TITRE="Rate of FA gravity variation today"
		T_OPTION="-T-1.0/1.0/0.2"         ! Range of the palette 
		Write(9,*) " "
		Write(9,*) "makecpt -Cno_green ", trim(adjustl(T_OPTION)), " > pale.cpt" 
		Write(9,*) "psbasemap -X3 -Y5 -Ba180/a40WSEn -Jm0.018i ", & 
           	 trim(adjustl(R_OPTION)), " -K > ", " ", trim(adjustl(nameout))  
		Write(9,*) "pscontour -I -Jm -O -K ", trim(adjustl(R_OPTION)), " ", & 
            	trim(adjustl(namein)), " -Cpale.cpt  >> ", trim(adjustl(nameout)) 
		If(option_rof=='r')then 
 	 	Write(9,*) "pscoast -Jm -Dc -B -W2/0   -A1000 -O -K ", trim(adjustl(R_OPTION)), " >> ", trim(adjustl(nameout))  
  		else
  		Write(9,*) "pscoast -Jm -Dc -B -W1/240 -A1000 -O -K ", trim(adjustl(R_OPTION)), " >> ", trim(adjustl(nameout))  
		Endif
		If(option_rof=='z')& 
		Write(9,*) "psxy -R -Jm ", trim(adjustl(filecap)), " -M -W4/0 -B -A -O -K >> ", trim(adjustl(nameout))  
		open (4,file=NAMEF,status='unknown')  
		If(CODE=='-1')    THEN 
		Write(4,*) "180 +84.2 26 0 3 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	 	  "  -ALMA rheology: ", trim(adjustl(SHORT_VISCO))  
		Write(4,*) "180 -87.5   13  0 0 BC ",& 
 	   	 " -LMAX=", trim(adjustl(DEGREE)),& 
 	    	" -RES=", trim(adjustl(RESOLUTION)),& 
	    	" -NV=", trim(adjustl(NV)),& 
	    	" -MODE=", trim(adjustl(MODE)),&   
	    	" -ITER=", trim(adjustl(ITER));  close(4) 
		Write(9,*) "pstext -N -R -Jm -B ", trim(adjustl(NAMEF)), " -G0 -O -K >> ", trim(adjustl(nameout)) 
		Write(9,*) "psscale -U/0.5/0.5/'SELEN 3.2' -E -Cpale.cpt -B1f0.25a0.5/:microgal/yr: -D8.25/-1/10/1h -O >> ", & 
	            trim(adjustl(nameout)) 
		ELSEIF(CODE/='-1')THEN 
		Write(4,*) "180 +84.2 26 0 3 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	   	"  -Viscosity profile: ", trim(adjustl(VSTRING))  
		Write(4,*) "180 -87.5   13  0 0 BC ",& 
 	    	" -LMAX=", trim(adjustl(DEGREE)),& 
 	    	" -RES=", trim(adjustl(RESOLUTION)),& 
	    	" -NV=", trim(adjustl(NV)),& 
	    	" -CODE=", trim(adjustl(CODE)),& 
	    	" -MODE=", trim(adjustl(MODE)),&   
	    	" -ITER=", trim(adjustl(ITER));  close(4) 
		Write(9,*) "pstext -N -R -Jm -B ", trim(adjustl(NAMEF)), " -G0 -O -K >> ", trim(adjustl(nameout)) 		
		Write(9,*) & 
                "psscale -U/0.5/0.5/'SELEN 3.2' -E -Cpale.cpt -B1f0.25a0.5/:microgal/yr: -D8.25/-1/10/1h -O -K >> ", & 
         	trim(adjustl(nameout))
 		Write(9,*) "psbasemap -J -R -B::/::EWSN -O >> ", trim(adjustl(nameout)) 
		ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	ENDIF    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
 	if    (OPTION_DOTSS=='y') THEN 
	        Write(9,*) ""
	        namein ="ssdotmap.dat"; nameout="ssdotmap.ps"; namef="tmpfss.dat"
	        TITRE="Rate of SS gravity variation today"
		T_OPTION="-T-1.0/1.0/0.2"         ! Range of the palette 
		Write(9,*) " "
		Write(9,*) "makecpt -Cno_green ", trim(adjustl(T_OPTION)), " > pale.cpt" 
		Write(9,*) "psbasemap -X3 -Y5 -Ba180/a40WSEn -Jm0.018i ", & 
           	 trim(adjustl(R_OPTION)), " -K > ", " ", trim(adjustl(nameout))  
		Write(9,*) "pscontour -I -Jm -O -K ", trim(adjustl(R_OPTION)), " ", & 
            	trim(adjustl(namein)), " -Cpale.cpt  >> ", trim(adjustl(nameout)) 
		If(option_rof=='r')then 
 	 	Write(9,*) "pscoast -Jm -Dc -B -W2/0   -A1000 -O -K ", trim(adjustl(R_OPTION)), " >> ", trim(adjustl(nameout))  
  		else
  		Write(9,*) "pscoast -Jm -Dc -B -W1/240 -A1000 -O -K ", trim(adjustl(R_OPTION)), " >> ", trim(adjustl(nameout))  
		Endif
		If(option_rof=='z')& 
		Write(9,*) "psxy -R -Jm ", trim(adjustl(filecap)), " -M -W4/0 -B -A -O -K >> ", trim(adjustl(nameout))  
		open (4,file=NAMEF,status='unknown')  
		If(CODE=='-1')    THEN 
		Write(4,*) "180 +84.2 26 0 3 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	 	  "  -ALMA rheology: ", trim(adjustl(SHORT_VISCO))  
		Write(4,*) "180 -87.5   13  0 0 BC ",& 
 	   	 " -LMAX=", trim(adjustl(DEGREE)),& 
 	    	" -RES=", trim(adjustl(RESOLUTION)),& 
	    	" -NV=", trim(adjustl(NV)),& 
	    	" -MODE=", trim(adjustl(MODE)),&   
	    	" -ITER=", trim(adjustl(ITER));  close(4) 
		Write(9,*) "pstext -N -R -Jm -B ", trim(adjustl(NAMEF)), " -G0 -O -K >> ", trim(adjustl(nameout)) 
		Write(9,*) "psscale -U/0.5/0.5/'SELEN 3.2' -E -Cpale.cpt -B1f0.25a0.5/:microgal/yr: -D8.25/-1/10/1h -O >> ", & 
	            trim(adjustl(nameout)) 
		ELSEIF(CODE/='-1')THEN 
		Write(4,*) "180 +84.2 26 0 3 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	   	"  -Viscosity profile: ", trim(adjustl(VSTRING))  
		Write(4,*) "180 -87.5   13  0 0 BC ",& 
 	    	" -LMAX=", trim(adjustl(DEGREE)),& 
 	    	" -RES=", trim(adjustl(RESOLUTION)),& 
	    	" -NV=", trim(adjustl(NV)),& 
	    	" -CODE=", trim(adjustl(CODE)),& 
	    	" -MODE=", trim(adjustl(MODE)),&   
	    	" -ITER=", trim(adjustl(ITER));  close(4) 
		Write(9,*) "pstext -N -R -Jm -B ", trim(adjustl(NAMEF)), " -G0 -O -K >> ", trim(adjustl(nameout)) 		
		Write(9,*) & 
                "psscale -U/0.5/0.5/'SELEN 3.2' -E -Cpale.cpt -B1f0.25a0.5/:microgal/yr: -D8.25/-1/10/1h -O -K >> ", & 
         	trim(adjustl(nameout))
 		Write(9,*) "psbasemap -J -R -B::/::EWSN -O >> ", trim(adjustl(nameout)) 
		ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	ENDIF    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

1000    CONTINUE 
!
 close(9) 		      
!
 END SUBROUTINE MAKE_GMAPS		      
!

!
!
!
 SUBROUTINE MAKE_WNW (RESOLUTION, DEGREE, FILE1_GMT)
 IMPLICIT NONE
 CHARACTER*3 RESOLUTION, DEGREE
 CHARACTER*20 R_OPTION, FILE1_GMT 
 CHARACTER*6 NPC
 CHARACTER*5 JMAXC
 INTEGER LMAX, NP, RES, JMAX, JMAXP  
!
! ---- Prepares a GMT script for plotting the window function
!         Last modification - Giorgio Spada October 2007 
!
 open(9,file=file1_gmt,status='unknown')
 Write(9,*) "gmtset PAPER_MEDIA A4+"
 Write(9,*) "gmtset ANOT_FONT_SIZE 16p"
 Write(9,*) "gmtset LABEL_FONT_SIZE 26p"
 Write(9,*) "#"
 Write(9,*) "# ---- Window function ----" 
 open(3,file='tmpw0.dat',status='unknown'); Write(3,'(a3)')degree; close(3) 
 open(3,file='tmpw0.dat',status='unknown'); Read(3,*)lmax; close(3) 
 open(3,file='tmpw1.dat',status='unknown'); Write(3,'(a3)')resolution; close(3) 
 open(3,file='tmpw1.dat',status='unknown'); Read(3,*)res; close(3) 
 NP=2*RES*(RES-1)*20+12 
 JMAX=(LMAX+1)*(LMAX+2)/2
 JMAXP=JMAX+30
 open(3,file='tmpw2.dat',status='unknown'); Write(3,'(i5)') jmaxp; close(3) 
 open(3,file='tmpw2.dat',status='unknown'); Read (3,'(a5)') jmaxc; close(3) 
 open(3,file='tmpw3.dat',status='unknown'); Write(3,'(i5)') NP;    close(3) 
 open(3,file='tmpw3.dat',status='unknown'); Read (3,'(a5)') NPC;   close(3) 
!
 R_OPTION="-R3.5/"//trim(adjustl(JMAXC))//"/1e-4/1e2"
!
 Write(9,*) "psbasemap -X6 -Y4 -U/0.5/12.75/'SELEN 3.2' ", R_OPTION, & 
            "-Ba1f1p:'degree*(degree+1)/2+order+1':/f1a10:@~e@~,%::,%::.'Window function':WSne -JX-18l/14l -K > wnw.ps"
 open (4,file='tmpp4.dat',status='unknown') 
 Write(4,*) " 4 3e1 16 0 2 BR -RES=", trim(adjustl(RESOLUTION)), & 
            " -NP=", trim(adjustl(NPC)), " -LMAX= ", trim(adjustl(DEGREE)), "  " ; close(4)
 Write(9,*) "pstext tmpp4.dat -N -JX -R -G0 -O -K >> wnw.ps"	     
 Write(9,*) "psxy wnw.dat -B -R -JX -Sc0.1 -G0 -K -O >> wnw.ps"
 open (4,file='tmpp5.dat',status='unknown') 
 Write(4,*) JMAXP, " 1"
 Write(4,*)  " 1 1 " ; CLOSE(4) 
 Write(9,*) "psxy tmpp5.dat -B -R -JX -W4ta -G0 -K -O >> wnw.ps"
 open (4,file='tmpp6.dat',status='unknown') 
 Write(4,*) JMAXP, " 5"
 Write(4,*)  " 1 5 " ; CLOSE(4) 
 Write(9,*) "psxy tmpp6.dat -B -R -JX -W4to -G0 -K -O >> wnw.ps"
 open (4,file='tmpp7.dat',status='unknown') 
 Write(4,*) JMAXP, " 10"
 Write(4,*)  " 1 10 " ; CLOSE(4) 
 Write(9,*) "psxy tmpp7.dat -B -R -JX -W2to -G0 -O >> wnw.ps"
!
 close(9) 
!
 END SUBROUTINE MAKE_WNW
!
!
!
!
!
!
 SUBROUTINE MAKE_PLOT_LDC (NV, CODE, DEGREE, VSTRING, FILE1_GMT)
 IMPLICIT NONE
 CHARACTER*3 NV, CODE
 CHARACTER*3 DEGREE
 CHARACTER*20 FILE1_GMT
 CHARACTER*22 R_OPTION1, R_OPTION2, R_OPTION3, R_OPTION4 
 CHARACTER*30 VSTRING
 INTEGER LMAX
!
! ---- Prepares the GMT script for plotting the LDCs and the spectrum 
!             Last modification Giorgio Spada November 2007 
!	      Revised  June 15, 2007 (V 2.6 INTEL) 
!             Revised May 2010 g95
!
!
! ---- Target GMT file for plotting the spectral properties of the chosen Earth model 
!
open(9,file=file1_gmt,status='unknown')
!	
Write(9,*) "gmtset PAPER_MEDIA A4+"
Write(9,*) "gmtset ANOT_FONT_SIZE 16p"
Write(9,*) "gmtset LABEL_FONT_SIZE 18p"
Write(9,*) "#"
!
! ---- The x-range varies according to the maximum degree but the y-range is fixed ... 
open(3,file='tmpg0.dat',status='unknown'); Write(3,'(a3)')degree; close(3) 
open(3,file='tmpg0.dat',status='unknown'); Read(3,*)LMAX; close(3) 
! for the elastic LDC  	
if(0<=LMAX.and.LMAX<=10)  R_OPTION1="-R0.8/10/-0.40/0.05"
if(0<=LMAX.and.LMAX<=100) R_OPTION1="-R0.8/100/-0.40/0.05"
if(LMAX>100) 		  R_OPTION1="-R0.8/1000/-0.40/0.05"
! for the fluid LDC
if(0<=LMAX.and.LMAX<=10)  R_OPTION2="-R0.8/10/-1.25/0.25"
if(0<=LMAX.and.LMAX<=100) R_OPTION2="-R0.8/100/-1.25/0.25"
if(LMAX>100) 		  R_OPTION2="-R0.8/1000/-1.25/0.25"
! for the Spectrum 
if(0<=LMAX.and.LMAX<=10)  R_OPTION3="-R0.8/10/1e1/1e9"
if(0<=LMAX.and.LMAX<=100) R_OPTION3="-R0.8/100/1e1/1e9"
if(LMAX>100) 		  R_OPTION3="-R0.8/1000/1e1/1e9"
! for the Normalized residues  
if(0<=LMAX.and.LMAX<=10)  R_OPTION4="-R0.8/10/1e-7/1e1"
if(0<=LMAX.and.LMAX<=100) R_OPTION4="-R0.8/100/1e-7/1e1"
if(LMAX>100) 		  R_OPTION4="-R0.8/1000/1e-7/1e1"
!
Write(9,*) "#"
Write(9,*) "# ---------- Elastic LDCs vs degree (left frame of 'ela-flu.ps')"
Write(9,*) "psbasemap -X5 -Y5 -U/8/-3/'SELEN 3.2' ", R_OPTION1, & 
	   " -Ba1f1p:'Harmonic degree, n':/f0.05a0.05:'Elastic LDC':WSen -JX9l/12 -K > ela-flu.ps" 
open (4,file='tmpg1.dat',status='unknown') 
Write(4,*) "150 0.09 18 0 2 BC Earth model: NV=", trim(adjustl(NV)), " CODE=", trim(adjustl(CODE))
Write(4,*) "150 0.07 16 0 2 BC Viscosity profile: ", trim(adjustl(VSTRING)), "*1E21 Pa.s"; close(4) 
Write(9,*) "pstext tmpg1.dat -N -JX -R -G0 -CM -O -K >> ela-flu.ps"
!
Write(9,*) "awk '{print $1, $2/(2*$1+1)}' hh.dat > h.tmp" 
Write(9,*) "awk '{print $1, $2         }' ll.dat > l.tmp"  
Write(9,*) "awk '{print $1, $2         }' kk.dat > k.tmp" 
Write(9,*) "psxy -H2 h.tmp -B -R -JX -Ss0.4  -G0  -K -O >> ela-flu.ps" 
Write(9,*) "psxy -H2 l.tmp -B -R -JX -Sc0.35 -G0  -K -O >> ela-flu.ps" 
Write(9,*) "psxy -H2 k.tmp -B -R -JX -Si0.4  -G0  -K -O >> ela-flu.ps" 
!
open (4,file='tmpg2.dat',status='unknown') 
Write(4,*) "20 -0.20 18 0 0 BL h/(2*n+1)"
Write(4,*) "20 -0.23 18 0 0 BL l" 
Write(4,*) "20 -0.26 18 0 0 BL k"; close(4) 
Write(9,*) "pstext tmpg2.dat -N -JX -R -G0 -CM -O -K >> ela-flu.ps"
!
open (4,file='tmpg3.dat',status='unknown') 
Write(4,*) "15 -0.20"; close(4) 
Write(9,*) "psxy -JX -R -G0 -Ss0.4  -O -K tmpg3.dat >> ela-flu.ps"
open (4,file='tmpg4.dat',status='unknown') 
Write(4,*) "15 -0.23"; close(4) 
Write(9,*) "psxy -JX -R -G0 -Sc0.35 -O -K tmpg4.dat >> ela-flu.ps"
open (4,file='tmpg5.dat',status='unknown') 
Write(4,*) "15 -0.26"; close(4)
Write(9,*) "psxy -JX -R -G0 -Si0.4  -O -K tmpg5.dat >> ela-flu.ps"
! 
Write(9,*) "#"
Write(9,*) "# ---------- Fluid LDCs vs degree (right frame of 'ela-flu.ps')"
Write(9,*) "psbasemap -X11 ", R_OPTION2, & 
	    " -Ba1f1p:'Harmonic degree, n':/f0.25a0.25:'Fluid LDC':wSEn -JX -K -O >> ela-flu.ps" 
Write(9,*) "awk '{print $1, $3/(2*$1+1)}' hh.dat > h.tmp" 
Write(9,*) "awk '{print $1, $3         }' ll.dat > l.tmp"  
Write(9,*) "awk '{print $1, $3         }' kk.dat > k.tmp" 
Write(9,*) "psxy -H2 h.tmp -B -R -JX -Ss0.4  -G0  -K -O >> ela-flu.ps" 
Write(9,*) "psxy -H2 l.tmp -B -R -JX -Sc0.35 -G0  -K -O >> ela-flu.ps"
Write(9,*) "psxy -H2 k.tmp -B -R -JX -Si0.4  -G0     -O >> ela-flu.ps"
!
Write(9,*) "#"
Write(9,*) "# ----  Relaxation spectrum (i. e., relaxation times vs. degree), 'spectrum.ps'"
Write(9,*) "psbasemap -X5 -Y5 -U/0.5/11/'SELEN 3.2' ", R_OPTION3, & 
	   " -Ba1f1p:'Harmonic degree, n':/f1a1p:'Relaxation time (yrs)':WSen -JX9l/12l  -K >  spectrum.ps" 
open (4,file='tmpg6.dat',status='unknown') 
Write(4,*) "10 0.8e10 18 0 2 BC Earth model: NV=", trim(adjustl(NV)), " CODE=", trim(adjustl(CODE))
Write(4,*) "10 0.3e10 16 0 2 BC Viscosity profile: ", trim(adjustl(VSTRING)), "*1E21 Pa.s"; close(4) 
Write(9,*) "pstext -N -JX -R -G0 -CM -O -K tmpg6.dat >> spectrum.ps"
Write(9,*) "awk '{print $1, $5}' ss.dat > spe.tmp" 
Write(9,*) "psxy -M -H7 spe.tmp -B -R -JX -Sc0.15 -G0 -O >> spectrum.ps" 
!
Write(9,*) "#"
Write(9,*) "# ---------- Normalized residues vs degree, all three in figure 'n-residues.ps'"
Write(9,*) "# ---------- h"
Write(9,*) "psbasemap  -X3 -Y4 -U/0.25/0.25/'SELEN 3.2' ", R_OPTION4, & 
	   " -Ba1f1p/f1a1p:'Normalized residue':WSen -JX7l/11l  -K > n-residues.ps" 
Write(9,*) "psxy -M -H2 ihh.dat -B -R -JX -Ss0.25 -W4 -G0   -O -K >> n-residues.ps" 
open (4,file='tmpg7.dat',status='unknown') 
Write(4,*) "1.3 2e0 20 0 0 BL (h)"; close(4) 
Write(9,*) "pstext -N tmpg7.dat -JX -R -G0 -CM -O -K >> n-residues.ps"
Write(9,*) "#"
Write(9,*) "# ---------- l"
Write(9,*) "psbasemap  -X8 -R -Ba1f1p:'Harmonic degree, n':/f1a1p:'Normalized residue':wSen -JX -K -O >> n-residues.ps"
Write(9,*) "psxy -M -H2 ill.dat -B -R -JX -Sc0.2 -W4 -G0 -O -K >> n-residues.ps"
open (4,file='tmpg8.dat',status='unknown') 
Write(4,*) "10 0.8e2 18 0 2 BC Earth model: NV=", trim(adjustl(NV)), " CODE=", trim(adjustl(CODE))
Write(4,*) "10 0.3e2 16 0 2 BC Viscosity profile: ", trim(adjustl(VSTRING)), "*1E21 Pa.s"; close(4) 
Write(9,*) "pstext -N tmpg8.dat -JX -R -G0 -CM -O -K >> n-residues.ps"
open (4,file='tmpg9.dat',status='unknown') 
Write(4,*) "1.3 2e0 20 0 0 BL (l)"; close(4)
Write(9,*) "pstext -N tmpg9.dat -JX -R -G0 -CM -O -K >> n-residues.ps"
Write(9,*) "#"
Write(9,*) "# ---------- k"
Write(9,*) "psbasemap -X8 -R -Ba1f1p/f1a1p:'Normalized residue':wSen -JX -O -K >> n-residues.ps" 
Write(9,*) "psxy -M -H2 ikk.dat -B -R -JX -Si0.25 -W4 -G0   -O -K >> n-residues.ps" 
open (4,file='tmpg10.dat',status='unknown') 
Write(4,*) "1.3 2e0 20 0 0 BL (k)"; close(4) 
Write(9,*) "pstext -N tmpg10.dat -JX -R -G0 -CM -O >> n-residues.ps"
!
 close(9)
!
 END SUBROUTINE make_plot_ldc
!
!
!
!
!
!
 SUBROUTINE MAKE_RSLZONES (NV, CD, RUN, NINC, NRSL, TITLICE, RESOLUTION, ITER, & 
 		          MODE, DEGREE, VSTRING, FILE1_GMT, FILE2_GMT, SHORT_VISCO_FILENAME)
!
! --- A GMT script for plotting the RSL zones - GS November 12 2007 
!     Revised JULY 04 2008 for version 2.6 of SELEN
!     Revised by GS April 2010 version 3.2 ALMA & g95  
!
 IMPLICIT NONE
 INTEGER, PARAMETER :: NREGIONS=4 
 CHARACTER*20 R_OPTION(NREGIONS), G_OPTION(NREGIONS), W_OPTION(NREGIONS) 
 CHARACTER*20 FILE1_GMT, FILE2_GMT
 CHARACTER*30 VSTRING
 CHARACTER*10 DEGREE, RESOLUTION, TITLICE   
 CHARACTER*3 STRING, NV, CD
 CHARACTER*4 RUN    
 CHARACTER*2 LABEL(NREGIONS)
 CHARACTER*3 NINC
 CHARACTER*1 ITER, MODE
 CHARACTER*100 SHORT_VISCO_FILENAME 
 INTEGER I, NRSL 
!
!
! ---- A range for each  region 
 DATA LABEL  /'10','11','20','21'/
!  
! ---- A colour for each region 
 DATA G_OPTION /" -G0/0/255", " -G107/199/231", " -G226/0/122", " -G105/62/142"/ 
!
! ---- A colour for the "average RSL curve" 
 DATA W_OPTION /" -W8ta/0/0/255", " -W8ta/107/199/231 ", " -W8ta/226/0/122 ", " -W8ta/105/62/142 "/ 
! 
! ---- An R option for each region 
        R_OPTION(1) = " -R0/"//trim(adjustl(ninc))//"/-200/50"
	R_OPTION(2) = " -R0/"//trim(adjustl(ninc))//"/-200/50"
	R_OPTION(3) = " -R0/"//trim(adjustl(ninc))//"/-50/500"
	R_OPTION(4) = " -R0/"//trim(adjustl(ninc))//"/-50/150"
!
!
!
! ======== PART 1 - a single frame for each zone 
!
! --- Target GMT file for plotting <<RSL zones>> 
!
 	open(19,file=file1_gmt,status='unknown')

	Write(19,*) " gmtset PAPER_MEDIA A4+" 
	Write(19,*) " gmtset FRAME_WIDTH 0.1c"
	Write(19,*) " " 
!Write(19,*) "echo '     - ", trim(adjustl(file1_gmt))//":", " Creating ps images of RSL zones - one image per frame....'"
Write(19,*) " " 
!
! ---- Here we consider 4 zones (10, 11, and 20, 21) Some of these may be empty... this depends 
!      on ice chronology, mantle viscosity and the range of harmonic degrees. Default input files 
!      are rslzones-**.dat and lonlat-**.dat. 
!
	DO I=1, NREGIONS 
!
! --- xy plot of RSL 
!
        Write(19,*) " psbasemap -U'SELEN 3.2' -X3 -Y10.5 -Ba2f1:'time (ka)':/a50f50WSen:'RSL (m)': ", & 
	              r_option(i), "-JX14/9  -K >  plot.ps" 
 	Write(19,*) " psxy rslzones-"//label(i)//".dat", " -M -H -B -R -JX -W0.05/0 -O -K >> plot.ps"
	Write(19,*) " psxy rslzones-"//label(i)//"-ave.dat", " -M -H -B -R -JX ", W_OPTION(I), " -O -K >> plot.ps"	
!
! --- first sphere 
 	Write (19,*) " psbasemap -Y0 -X15 -Ba90/a85f90WSEN -R0/360/-90/90  -JG-20/70/9 -O -K  >> plot.ps"
 	Write (19,*) " psxy lonlat-"//label(i)//".dat", " -B -R -JG -Sh0.32 ", g_option(i), " -O -K >> plot.ps"
 	Write (19,*) " pscoast -G255 -B -R -JG -Dc -W2 -A10000 -O -K >> plot.ps"
!
! --- second sphere 
 	Write (19,*) " psbasemap -Y-10 -X0 -B -R -JG-20/-70/9 -O -K  >> plot.ps"
 	Write (19,*) " psxy lonlat-"//label(i)//".dat", " -B -R -JG -Sh0.32 ", g_option(i), " -O -K >> plot.ps"
 	Write (19,*) " pscoast -G255 -B -R -JG -Dc -W2 -A10000 -O -K >> plot.ps"
!
! --- a new basemap for titles... 
	Write (19,*) " psbasemap -X-15 -Y10 -Bf1000wesn   -R0/10/0/10     -JX14/9  -O -K >>  plot.ps"
!
	If(cd=='-1')then 
 	open (8,file='tmpz'//"."//label(i),status='unknown')	
		Write(8,*) ".5 -5  26 0 0 BL  RSL curves for zone ", label(i)
		Write(8,*) ".5 -6 16  0 0 BL  Ice model: ", trim(adjustl(titlice))
		Write(8,*) ".5 -7 16  0 0 BL  Repository: ", "./depot-/", trim(adjustl(run))
		Write(8,*) ".5 -8 16 0 0 BL  -LMAX =", trim(adjustl(DEGREE)), & 
					       " -RES =", trim(adjustl(RESOLUTION)), & 
					       " -ITER =", trim(adjustl(ITER)), & 
					       " -Mode =", trim(adjustl(MODE)) 
		Write(8,*) ".5 -9 16 0 0 BL -ALMA model =", trim(adjustl(SHORT_VISCO_FILENAME))   
 	close(8) 
	Else 
 	open (8,file='tmpz'//"."//label(i),status='unknown')	
		Write(8,*) ".5 -5  26 0 0 BL  RSL curves for zone ", label(i)
		Write(8,*) ".5 -6 16  0 0 BL  Ice model: ", trim(adjustl(titlice))
		Write(8,*) ".5 -7 16  0 0 BL  Repository: ", "./depot-/", trim(adjustl(run))
		Write(8,*) ".5 -8 16  0 0 BL  Viscosity profile: ", trim(adjustl(VSTRING))
		Write(8,*) ".5 -9 16 0 0 BL  -LMAX =", trim(adjustl(DEGREE)), & 
					       "-RES =", trim(adjustl(RESOLUTION)), & 
					       " -ITER =", trim(adjustl(ITER)), & 
					       " -Mode =", trim(adjustl(MODE)), & 
					       " -NV =", trim(adjustl(NV)), & 
					       " -CODE =", trim(adjustl(CD))			  	         
 	close(8) 	
	Endif 		
!		
	Write (19,*) " pstext -N ", 'tmpz'//"."//label(i), " -JX -R -O >> plot.ps"
!
	Write (19,*) " ps2pdf plot.ps"	
!
	Write (19,*) " mv plot.pdf rslzones-"//label(i)//".pdf"	
!
	Write (19,*) " mv plot.ps rslzones-"//label(i)//".ps"	
!
	Write (19,*) " "	
	
	ENDDO 
!
!
! ======== PART 2 - a frame with all zones
!
! --- Target GMT file for plotting <<RSL zones>> 
!
 	open(19,file=file2_gmt,status='unknown')

	Write(19,*) " gmtset PAPER_MEDIA A4+" 
	Write(19,*) " gmtset FRAME_WIDTH 0.1c"
	Write(19,*) " gmtset ANOT_FONT_SIZE 12p"
	Write(19,*) " gmtset LABEL_FONT_SIZE 18p"
	Write(19,*) " " 
!Write(19,*) & 
!"echo '     - ", trim(adjustl(file2_gmt))//":", " Creating ps images of RSL zones - one images in one frame....'"
Write(19,*) " " 
!
! --- xy plots of RSL zones  
!
        Write(19,*) " psbasemap -X3 -Y14 -Ba2f1:'time (ka)':/a50f50Wsen:'RSL (m)':  ", r_option(1), "-JX6.5/5  -K >  plot.ps" 
 	Write(19,*) " psxy rslzones-"//label(1)//".dat", " -M -H -B -R -JX -Sc0.01 -G0/0/255 -O -K >> plot.ps"	
	Write(19,*) " psxy rslzones-"//label(1)//"-ave.dat", " -M -H -B -R -JX ", W_OPTION(1), " -O -K >> plot.ps"
 	open (8,file='tmpzu1.dat',status='unknown')	
		Write(8,*) "1 -180 18 0 2 BL  zone 10"  	         
 	close(8) 
	Write (19,*) " pstext -N tmpzu1.dat -JX -R -O -K >> plot.ps"
!
        Write(19,*) " psbasemap -X8.25 -Y0 -Ba2f1:'time (ka)':/a50f50Wsen  ", r_option(2), "-JX6.5/5  -O -K >>  plot.ps" 
 	Write(19,*) " psxy rslzones-"//label(2)//".dat", " -M -H -B -R -JX -Sc0.01 -G107/199/231 -O -K >> plot.ps"
	Write(19,*) " psxy rslzones-"//label(2)//"-ave.dat", " -M -H -B -R -JX ", W_OPTION(2), " -O -K >> plot.ps"
 	open (8,file='tmpzu2.dat',status='unknown')	
		Write(8,*) "1 -180 18 0 2 BL  11"  	         
 	close(8) 
	Write (19,*) " pstext -N tmpzu2.dat -JX -R -O -K >> plot.ps"	
!
        Write(19,*) " psbasemap -U'SELEN 3.2' -X-8.25 -Y-6  -Ba2f1:'time (ka)':/a100f100WSen:'RSL (m)':  ", & 
	              r_option(3), "-JX6.5/5  -O -K >>  plot.ps" 
 	Write(19,*) " psxy rslzones-"//label(3)//".dat", " -M -H -B -R -JX -Sc0.01 -G226/0/122 -O -K >> plot.ps"
	Write(19,*) " psxy rslzones-"//label(3)//"-ave.dat", " -M -H -B -R -JX ", W_OPTION(3), " -O -K >> plot.ps"
 	open (8,file='tmpzu3.dat',status='unknown')	
		Write(8,*) "1 400 18 0 2 BL  20"  	         
 	close(8) 
	Write (19,*) " pstext -N tmpzu3.dat -JX -R -O -K >> plot.ps"		
!
        Write(19,*) " psbasemap -X8.25 -Y0    -Ba2f1:'time (ka)':/a50f50WSen  ", r_option(4), "-JX6.5/5  -O -K >>  plot.ps" 
 	Write(19,*) " psxy rslzones-"//label(4)//".dat", " -M -H -B -R -JX -Sc0.01 -G105/62/142  -O -K  >> plot.ps"
	Write(19,*) " psxy rslzones-"//label(4)//"-ave.dat", " -M -H -B -R -JX ", W_OPTION(4), " -O -K  >> plot.ps"
 	open (8,file='tmpzu4.dat',status='unknown')	
		Write(8,*) "1 115 18 0 2 BL  21"  	         
 	close(8) 
	Write (19,*) " pstext -N tmpzu4.dat -JX -R -O -K >> plot.ps"		
!
! --- first sphere 
 	Write (19,*) " psbasemap -X7.5 -Y3 -Ba90/a85f90WSEN -R0/360/-90/90  -JG-20/70/9 -O -K  >> plot.ps"
 	Write (19,*) " psxy lonlat-"//label(1)//".dat", " -B -R -JG -Sh0.32 ", g_option(1), " -O -K >> plot.ps"
 	Write (19,*) " psxy lonlat-"//label(2)//".dat", " -B -R -JG -Sh0.32 ", g_option(2), " -O -K >> plot.ps"
 	Write (19,*) " psxy lonlat-"//label(3)//".dat", " -B -R -JG -Sh0.32 ", g_option(3), " -O -K >> plot.ps"
 	Write (19,*) " psxy lonlat-"//label(4)//".dat", " -B -R -JG -Sh0.32 ", g_option(4), " -O -K >> plot.ps"
 	Write (19,*) " pscoast -G255 -B -R -JG -Dc -W2 -A10000 -O -K >> plot.ps"
!
! --- second sphere 
 	Write (19,*) " psbasemap -Y-10 -X0 -B -R -JG-20/-70/9 -O -K  >> plot.ps"
 	Write (19,*) " psxy lonlat-"//label(1)//".dat", " -B -R -JG -Sh0.32 ", g_option(1), " -O -K >> plot.ps"
 	Write (19,*) " psxy lonlat-"//label(2)//".dat", " -B -R -JG -Sh0.32 ", g_option(2), " -O -K >> plot.ps"
 	Write (19,*) " psxy lonlat-"//label(3)//".dat", " -B -R -JG -Sh0.32 ", g_option(3), " -O -K >> plot.ps"
 	Write (19,*) " psxy lonlat-"//label(4)//".dat", " -B -R -JG -Sh0.32 ", g_option(4), " -O -K >> plot.ps"
 	Write (19,*) " pscoast -G255 -B -R -JG -Dc -W2 -A10000 -O -K >> plot.ps"
!
! --- a new basemap for titles... 
	If(cd=='-1')then 
	Write (19,*) " psbasemap -X-7.5 -Y7 -Bf1000wesn   -R0/10/0/10     -JX6.5/5  -O -K >>  plot.ps"
 	open (8,file='tmpzu.dat',status='unknown')	
		Write(8,*) "-9 -7  18  0  0 BL  RSL zones "
		Write(8,*) "-9 -9  14  0 0 BL  Ice model: ", trim(adjustl(titlice))
		Write(8,*) "-9 -10.5 14  0 0 BL  Repository: ./depot-", trim(adjustl(run)) 
		Write(8,*) "-9 -12.0 14  0 0 BL  -LMAX =", trim(adjustl(DEGREE)), & 
					       " -RES =", trim(adjustl(RESOLUTION)), & 
					       " -ITER =", trim(adjustl(ITER)), & 
					       " -Mode =", trim(adjustl(MODE))
		Write(8,*) "-9 -13.5 14 0 0  BL -ALMA model =", trim(adjustl(SHORT_VISCO_FILENAME))   
 	close(8) 	
	Else 
	Write (19,*) " psbasemap -X-7.5 -Y7 -Bf1000wesn   -R0/10/0/10     -JX6.5/5  -O -K >>  plot.ps"
 	open (8,file='tmpzu.dat',status='unknown')	

		Write(8,*) "-9 -7  18 0  0 BL  RSL zones "
		Write(8,*) "-9 -9  14  0 0 BL  Ice model: ", trim(adjustl(titlice))
		Write(8,*) "-9 -10.5 14  0 0 BL  Repository: ./depot-", trim(adjustl(run)) 
		Write(8,*) "-9 -12.0 14  0 0 BL  Viscosity profile: ", trim(adjustl(VSTRING))
		Write(8,*) "-9 -13.5 14  0 0 BL  -LMAX =", trim(adjustl(DEGREE)), & 
					       " -RES =", trim(adjustl(RESOLUTION)), & 
					       " -ITER =", trim(adjustl(ITER)), & 
					       " -Mode =", trim(adjustl(MODE)), & 
					       " -NV =", trim(adjustl(NV)), & 
					       " -CODE =", trim(adjustl(CD))			  	         
 	close(8) 	
	Endif 		
!
!
	Write (19,*) " pstext -N tmpzu.dat -JX -R -O >> plot.ps"
!	
	Write (19,*) " ps2pdf plot.ps"	
!
	Write (19,*) " mv plot.ps rslzones-all.ps"	
!
	Write (19,*) " mv plot.pdf rslzones-all.pdf"	
!
	close(19) 
!	
 END SUBROUTINE MAKE_RSLZONES
!
!
!
!
!
!
 SUBROUTINE MAKE_RSLMIS (NV, CD, RUN, NINC, NRSL, TITLICE, RESOLUTION, ITER, & 
 		      MODE, DEGREE, VSTRING, FILE_GMT, SHORT_VISCO_FILENAME)
!
! --- Prepares a GMT script for plotting a simple histogram
!     showing the RSL misfit distribution - GS January 23 2008 
!     Revised July 2008 for version 2.6 of SELEN
!     Revised April 2010 GS for version g95 
!
 IMPLICIT NONE
 CHARACTER*30 VSTRING
 CHARACTER*20 FILE_GMT
 CHARACTER*10 DEGREE, RESOLUTION, TITLICE   
 CHARACTER*4 RUN  
 CHARACTER*3 NV, CD 
 CHARACTER*3 NINC
 CHARACTER*1 ITER, MODE
 CHARACTER*100 SHORT_VISCO_FILENAME
 INTEGER NRSL 
!
!
! --- Target GMT file for misfit histogram   
!
 	open(19,file=file_gmt,status='unknown')

	Write(19,*) "gmtset PAPER_MEDIA A4+" 
	Write(19,*) "gmtset FRAME_WIDTH 0.1c"
!	Write(19,*) "echo       - Creating a postscript image of the misfit histogram..."	
!
  	Write(19,*) "psbasemap -X3 -Y3 -U'SELEN 3.2' -JX16/10  -G255 -R0/100/0/50 -K > rsl-misfit.ps" 
	Write(19,*) "awk '{print $1, $2}' mis.dat > mis.tmp"
   	Write(19,*) "pshistogram mis.tmp -JX -R -B100a20:misfit:/a10:frequency:WSen -G140 -W2 -O -K >> rsl-misfit.ps" 
	Write(19,*) "/bin/rm mis.tmp"
!
 	open (8,file='tmph.dat',status='unknown')		  	         
	If(cd=='-1')then 
	Write(8,*) "45 40 12 0 0 BL  Ice model: ", trim(adjustl(titlice))
	Write(8,*) "45 37 12 0 0 BL  Repository: ./depot-", trim(adjustl(run))
	Write(8,*) "45 34 10 0 0 BL  -LMAX =", trim(adjustl(DEGREE)), & 
				   " -RES =",  trim(adjustl(RESOLUTION)), & 
				   " -ITER =", trim(adjustl(ITER)), & 
				   " -Mode =", trim(adjustl(MODE))					  
	Write(8,*) "45 31 10 0 0 BL -ALMA model =", trim(adjustl(SHORT_VISCO_FILENAME))
	Write(8,*) "45 28 10 0 2 BL  Number of sites =", NRSL	
 	close(8) 
!
	Else 
!
	Write(8,*) "45 40 12 0 0 BL  Ice model: ", trim(adjustl(titlice))
	Write(8,*) "45 37 12 0 0 BL  Repository: ./depot-", trim(adjustl(run))
	Write(8,*) "45 34 12 0 0 BL  Viscosity profile: ", trim(adjustl(VSTRING))
	Write(8,*) "45 31 10 0 0 BL  -LMAX =", trim(adjustl(DEGREE)), & 
				   " -RES =",  trim(adjustl(RESOLUTION)), & 
				   " -ITER =", trim(adjustl(ITER)), & 
				   " -Mode =", trim(adjustl(MODE))					  
	Write(8,*) "45 28 10 0 0 BL    -NV =", trim(adjustl(NV)),     & 
		       	           " -CODE =", trim(adjustl(CD))

	Write(8,*) "45 25 10 0 2 BL  Number of sites =", NRSL	
 	close(8) 
!	
	Endif
!
	Write(19,*) "pstext tmph.dat -O -JX -R >> rsl-misfit.ps"	
	close(19) 
!	
 END SUBROUTINE MAKE_RSLMIS
!
!
!
!
!
!
 SUBROUTINE MAKE_RSLDB (DATAFILE, RSL_FORMAT, NRSL, FILE_GMT)
 IMPLICIT NONE
!
! Writes a GMT script for plotting a map of RSL sites distribution 
! The output filename is "File_GMT"   --- GS October 27 2007 ---
! Also modified in August 2009 for the implementation of SELEN 3
! *** Revised GS July 2010 - g95 - Double precision implementation 
! *** Revised GS Dec 2010 - RSL database of type "3" (for data from Dorit Sivan) 
!
 CHARACTER*1   CJUNK, RSL_FORMAT 
 CHARACTER*100 SS(100)
 CHARACTER*10  LATSC10, LONSC10 
 CHARACTER*50  LINE
 CHARACTER*20  FILE_GMT
 CHARACTER*30  DATAFILE 
 INTEGER I, J, NRSL, NOUT, CODE, NDATA
 REAL*8 LONS, LATS
 CHARACTER*200 LINEP 
! 
! new new new new new new new new new
 CHARACTER*100 ANOTHERROW
 INTEGER, PARAMETER :: VERY_LARGE_INTEGER = 100000
 INTEGER K, HEADER_LINES 
! new new new new new new new new new
!
!
! --- Open and reads the RSL database for getting info about lon/lat of data... 
!
!
! --- Target file for lon-lat of RSL sites 
 	open(37,file='lon-lat-rsl.dat',status='unknown')
!
! --- An header 
	Write(37,*) "Longitudes and latitudes of RSL sites (degrees)"
!
!
!
	If(RSL_FORMAT=='0') then 
!
        OPEN(1,FILE=DATAFILE,STATUS='unknown')
   	do i=1, nrsl 	
		READ (1, 6) code, lats, lons, ndata
		if(lons<=0.) lons=lons+360. 			 			 		 			 
		write(37,'(F6.1, 1X, F6.1)') lons, lats
		do j=1, ndata
			read(1,*) LINE			
		enddo
	enddo
	CLOSE(1) 
	CLOSE(37) 
6       FORMAT(1X,I3,1X,F5.1,1X,F6.1,1X,I2,1X,A22)
!
! NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW 
! NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW 
!
	ELSEIF(RSL_FORMAT=='3') then
!
! ...   Counting the header lines 
!
        open(1,file=DATAFILE,status='old') 
        HEADER_LINES=0 
        do k=1, VERY_LARGE_INTEGER
			read(1,'(a80)',end=23132)  ANOTHERROW
			if(ANOTHERROW(1:1)=='!') HEADER_LINES = HEADER_LINES + 1
		enddo
23132   close(1) 
!
        OPEN(1,FILE=DATAFILE,STATUS='old')
        do k=1, HEADER_LINES
			read(1,'(a80)')  ANOTHERROW	
        enddo
   	do i=1, nrsl 	
		READ (1, *) code
		READ(1,'(A80)')  ANOTHERROW
	        READ(1,'(a80)')  ANOTHERROW	  		       
		READ(1, *) lons, lats, ndata
		if(lons<=0.) lons=lons+360. 			 			 		 			 
		write(37,'(F6.1, 1X, F6.1)') lons, lats
		do j=1, ndata
			read(1,'(a80)')  ANOTHERROW			
		enddo
	enddo
	CLOSE(1) 
	CLOSE(37) 
!
! NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW 
! NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW 
!
	ELSEIF(RSL_FORMAT=='1') then
!
        OPEN(1,FILE=DATAFILE,STATUS='unknown')
!	
   	DO I=1, NRSL 
			READ(1,'(A1)')CJUNK 
			READ(1,'(A1)')CJUNK 
			READ(1,'(A1)')CJUNK 
			        READ(1,*) LATS, LONS  
                        	IF(LONS<=0.) LONS=LONS+360. 	
			READ(1,'(A1)')CJUNK 
			READ(1,'(A1)')CJUNK 
 			WRITE(37,'(F9.4, 1X, F9.4)') LONS, LATS
	ENDDO 
	CLOSE(1) 
	CLOSE(37) 
!
	ELSEIF(RSL_FORMAT=='2') then
!
        OPEN(1,FILE=DATAFILE,STATUS='unknown')
!	
   	DO I=1, NRSL 
			READ(1,'(A1)')CJUNK 
			READ(1,'(A1)')CJUNK 
			READ(1,'(A1)')CJUNK 			
			read(1,'(a200)') linep 	
!		
			call scan_string (linep, 2, ss, nout)				
!
			LATSC10=trim(adjustl(ss(1))) 
			LONSC10=trim(adjustl(ss(2)))			
!
 			call FROMDDMMSS_2_DEGREES (LATSC10, LATS) 
 			call FROMDDMMSS_2_DEGREES (LONSC10, LONS) 
!
                        IF(LONS<=0.) LONS=LONS+360. 	
			READ(1,'(A1)')CJUNK 
			READ(1,'(A1)')CJUNK 
 			WRITE(37,'(F9.4, 1X, F9.4)') LONS, LATS
	ENDDO 
	CLOSE(1) 
	CLOSE(37) 



        Endif 
!
!
! -------------------------------------------------------------
! --- Creates a GMT script for plotting a map of RSL sites  --- 
! -------------------------------------------------------------
!
! --- Target GMT file for a map of RSL sites 
!
 	open(29,file=file_gmt,status='unknown')
!
	Write(29,*) "gmtset PAPER_MEDIA A4+" 
	Write(29,*) "gmtset HEADER_FONT_SIZE 24p"
	Write(29,*) "gmtset FRAME_WIDTH 0.1c"
        Write(29,*) "gmtset ANOT_FONT_SIZE 12p"
	Write(29,*) ""
	Write(29,*) ""

	open(22,file='tmptitle',status='unknown') 
	Write(22,*) "180 87.4 22 0 1 BC Relative Sea Level (RSL) sites"
	Write(22,'(a29,a12,a7,i4)') "180 86.2 16 0 2 BC from file ", trim(adjustl(datafile)), " -NRSL=", NRSL 			
	close(22) 
!
! ------ Main frame (Mercator projection) 
!
	Write(29,*) "psbasemap -X4 -Y5 -Ba180/a85f90WSEn -R0/360/-85/85  -JM12 -K > maprsl.ps" 	
	Write(29,*) "pscoast -G120 -S0/0/220 -B -R -O -K  -JM -Dc  -A10000 >> maprsl.ps" 	
	Write(29,*) "psxy -H1 lon-lat-rsl.dat -B -R -JM -Sc0.1 -G220 -O -K >> maprsl.ps"  
	Write(29,*) "pstext -N tmptitle", " -G0 -O -K -JM -R >> maprsl.ps"
	Write(29,*) ""
!
! ------ Northern emisphere 
!	     	
 	Write(29,*) "psbasemap -Y6 -X14 -Ba90/a85f90WSEN -R0/360/0/90  -JG-45/90/8 -O -K  >> maprsl.ps"  
	Write(29,*) "pscoast -G120 -S0/0/220 -B -R -O -K -JG -Dc  -A10000 >> maprsl.ps" 	
	Write(29,*) "psxy -H1 lon-lat-rsl.dat -B -R -JG -Sc0.1 -G220 -O -K >> maprsl.ps"  
	Write(29,*) ""
!
! ------Southern emisphere 
!	     	
 	Write(29,*) "psbasemap -Y-9 -X0 -Ba90/a85f90WSEN -R0/360/-90/0  -JG-90/-90/8 -U/-10/1/'SELEN 3.2' -O -K  >> maprsl.ps"  
	Write(29,*) "pscoast -G120 -S0/0/220 -B -R -O -K -JG -Dc  -A10000 >> maprsl.ps" 	
	Write(29,*) "psxy -H1 lon-lat-rsl.dat -B -R -JG -Sc0.1 -G220 -O -K >> maprsl.ps"  
	Write(29,*) ""
!
	close(29) 	
!	
 END SUBROUTINE MAKE_RSLDB
!
!
!
!
!
!
 SUBROUTINE MAKE_RSL (NV, CD, DATAFILE, RSL_FORMAT, RUN, NINC, NRSL, TITLICE, & 
 	              RESOLUTION, ITER, MODE, DEGREE, VSTRING, FILE_GMT,      & 
		      SHORT_VISCO_FILENAME)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! Writes a GMT script for plotting the observed RSL and predictions 
! --- GS February 21, 2008 ---
! Revised July 2008 for v 2.6 
! Revised August 2008 for v 2.7 
! Also revised August 2009 for v 3
! Updated on April 2010 by GS for v. 3.1 (g95)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! *** Revised GS July 2010 - g95 - Double precision implementation
! *** Revised GS Dec 2010 - RSL database of type "3" (for data from Dorit Sivan) 
!
 IMPLICIT NONE
 CHARACTER*20 R_OPTION
 CHARACTER*50 B_OPTION
 CHARACTER*30 VSTRING
 CHARACTER*20 FILE_GMT, YMINC, YMAXC, DELTAC
 CHARACTER*30 DATAFILE
 CHARACTER*10 RESOLUTION, DEGREE, TITLICE 
 CHARACTER*4 RUN, STRING, CJUNK     
 CHARACTER*3 NV, CD
 CHARACTER*3 NINC
 CHARACTER*1 ITER, MODE, RSL_FORMAT 
 INTEGER I, J, K, NRSL 
 INTEGER, PARAMETER :: MAXSITES=1000, MAXDATA=1000
 CHARACTER*50 TITRE(MAXSITES), TITREA(MAXSITES), TITREB(MAXSITES)
 INTEGER CODE(MAXSITES), NDATA(MAXSITES)
 REAL*8    TIME(MAXSITES,0:MAXDATA)    	! times bp for which data are available at site#k
 REAL*8   DTIME(MAXSITES,0:MAXDATA)    	! uncertainties on the above times...
 REAL*8     RSL(MAXSITES,0:MAXDATA)    	! rsl datum for site k ad the times above
 REAL*8    DRSL(MAXSITES,0:MAXDATA)    	! uncertainties on the rsl datum
 REAL*8 AJUNK, LONS(MAXSITES), LATS(MAXSITES)
 REAL*8 DMIN, DMAX, RANGE, DELTA 
 CHARACTER*3 CODEC(MAXSITES)
 CHARACTER*10 LATSC10, LONSC10 
 CHARACTER*100 SS(2)
 CHARACTER*100 SHORT_VISCO_FILENAME  
 CHARACTER*200 LINEP 
 INTEGER NOUT
!
! new new new new new new new new new new new new 
!
 REAL*8,  PARAMETER :: time_today=2011. 
 INTEGER, PARAMETER :: VERY_LARGE_INTEGER=10000
 INTEGER III, KKK, N_HEADER_LINES  
 CHARACTER*80  ANOTHERROW
 CHARACTER*40  TITRE1, TITRE2
!
! new new new new new new new new new new new new 
!
!
!
! ------ Open and reads the sealevel database for getting info about data... 
!

	If    (RSL_FORMAT=='0') then 
!  
        do k=1, 2 
 		OPEN(1, FILE=DATAFILE,STATUS='unknown')
   	do i=1, nrsl 	
 		if(k==1) then 
			 READ (1, 6) code(i), lats(i), lons(i), ndata(i)
6  			 FORMAT(1X,I3,1X,F5.1,1X,F6.1,1X,I2,1X,A22)
			 if(lons(i)<=0.) lons(i)=lons(i)+360.
!			 
			 open(33,file='junk.dat',status='unknown'); write(33,'(i3)')  code(i); close(33) 
			 open(33,file='junk.dat',status='unknown'); read (33,'(a3)') codec(i); close(33) 			 			 
!			 			 
			 else
			 READ (1, '(A50)') TITRE(I)
			 endif 
!			 
		do j=1, ndata(i)
			read(1,*) time(i,j), dtime(i,j), rsl(i,j), drsl(i,j)
		enddo
	enddo
!
	CLOSE(1) 
	enddo
!
!
	Elseif(RSL_FORMAT=='3') then
!
!
! NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW DEC 2010
! NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW DEC 2010  
! 
        OPEN(1, FILE=DATAFILE,STATUS='unknown')
!
        N_HEADER_LINES=0 
        do k=1, VERY_LARGE_INTEGER
			read(1,'(a80)',end=23132)  ANOTHERROW
			if(ANOTHERROW(1:1)=='!') N_HEADER_LINES = N_HEADER_LINES + 1
		enddo
23132   close(1) 
!
        OPEN(1, FILE=DATAFILE,STATUS='unknown')
	do k=1, N_HEADER_LINES 
		read(1,'(a80)') ANOTHERROW
	enddo
!
	do 77788 i=1, NRSL 
!
		READ (1,       *) code(i)
        	open(33,file='xjunk.dat',status='unknown'); write(33,'(i3)')  code(i); close(33) 
        	open(33,file='xjunk.dat',status='unknown'); read (33,'(a3)') codec(i); close(33) 			 			 
!
		READ (1, '(a40)') TITRE1 
!
	        READ (1, '(a80)') ANOTHERROW	  		       
!
		READ (1, '(a40)') TITRE2 
        	open(33,file='yjunk.dat',status='unknown'); write(33,'(a40)') TITRE2                    ; close(33) 
        	open(33,file='yjunk.dat',status='unknown'); read (33, *     ) lons(i), lats(i), ndata(i); close(33) 			 			 
!
		TITRE(I) = codec(i)//" "//trim(adjustl(TITRE1))//" "//trim(adjustl(TITRE2)) 
!	
        	do j=1, ndata(i)
            		read(1,*) iii, kkk, time(i,j), dtime(i,j), rsl(i,j), drsl(i,j)
			if(iii==1) time(i,j) = time_today - time(i,j) 
			if(iii==2) time(i,j) =              time(i,j)		        
        	enddo
!
77788   CONTINUE

        CLOSE(1)                                                                               
!
! NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW DEC 2010
! NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW DEC 2010
!
!
	Elseif(RSL_FORMAT=='1') then
!
 		OPEN(1, FILE=DATAFILE,STATUS='unknown')
!
   		DO I = 1, NRSL 
		    NDATA(I)=1
		    READ(1,*)CODE(I)
			 open(33,file='junk.dat',status='unknown'); write(33,'(i3)')  code(i); close(33) 
			 open(33,file='junk.dat',status='unknown'); read (33,'(a3)') codec(i); close(33)			
!
	            READ(1,'(A50)')TITREA(I) 
		    READ(1,'(A50)')TITREB(I)
		         TITRE(I)=TRIM(ADJUSTL(TITREA(I)))//" - "//TRIM(ADJUSTL(TITREB(I))) 
	            READ(1,*) LATS(I), LONS(I)  
                         IF(LONS(I)<=0.) LONS(I)=LONS(I)+360. 					
		    READ(1,*) TIME(I,1), DTIME(I,1), RSL(I,1), DRSL(I,1)
		    READ(1,'(A1)')CJUNK 
		ENDDO 
		CLOSE(1) 
!
	Elseif(RSL_FORMAT=='2') then
!
 		OPEN(1, FILE=DATAFILE,STATUS='unknown')
!
   		DO I = 1, NRSL 
		    NDATA(I)=1
		    READ(1,*)CODE(I)
			 open(33,file='junk.dat',status='unknown'); write(33,'(i3)')  code(i); close(33) 
			 open(33,file='junk.dat',status='unknown'); read (33,'(a3)') codec(i); close(33)			
!
	            READ(1,'(A50)')TITREA(I) 
		    READ(1,'(A50)')TITREB(I)
		         TITRE(I)=TRIM(ADJUSTL(TITREA(I)))//" - "//TRIM(ADJUSTL(TITREB(I))) 
!
			READ(1,'(A200)') LINEP 			
			CALL SCAN_STRING (LINEP, 2, SS, NOUT)
!	
			LATSC10=trim(adjustl(SS(1))) 
			LONSC10=trim(adjustl(SS(2)))				
! 
 			call FROMDDMMSS_2_DEGREES (LATSC10, LATS(I)) 
 			call FROMDDMMSS_2_DEGREES (LONSC10, LONS(I)) 
!
                         IF(LONS(I)<=0.) LONS(I)=LONS(I)+360. 					
		    READ(1,*) TIME(I,1), DTIME(I,1), RSL(I,1), DRSL(I,1)
		    READ(1,'(A1)')CJUNK 
		ENDDO 
		CLOSE(1) 
!
	ENDIF
!
!
!
! ------ Part 1: Creates a GMT script for RSL curves  
!
! ------ Target GMT file for RSL curves  
!
 	open(19,file=file_gmt,status='unknown')

	Write(19,*) "gmtset PAPER_MEDIA A4+" 
	Write(19,*) "gmtset FRAME_WIDTH 0.1c"
!
!Write(19,*) & 
!"echo '--- ", trim(adjustl(file_gmt))//":", " Creating ps images of RSL data and predictions'"
!
	DO 11111 I=1, NRSL
!
! ------- Based on the sealevel database, finds a range for plotting
!
	dmax=-99999.
	dmin=+99999.
	do j=1, ndata(i) 
		if(rsl(i,j)<=dmin) dmin=rsl(i,j) 
		if(rsl(i,j)>=dmax) dmax=rsl(i,j)  	
	enddo
	if(dmin<=0.and.dmax<=0) dmax=0 
	if(dmin>=0.and.dmax>=0) dmin=0 
!	 	
	range=dmax-dmin 
!	
	if(range>=0  .and.range<=10)  delta=1
	if(range>=10 .and.range<=50)  delta=5	  
	if(range>=50 .and.range<100) delta=10	
	if(range>=100.and.range<200) delta=20
	if(range>=200.and.range<300) delta=50
	if(range>=300.and.range<400) delta=100
	if(range>=400) 		  delta=150 
!
   	OPEN  (10,FILE='junky'//trim(adjustl(codec(i))),STATUS='unknown')
	WRITE (10,'(i10)') int(dmin-delta)
	WRITE (10,'(i10)') int(dmax+delta)
	WRITE (10,'(i10)') int(delta) 	
	CLOSE(10)
   	OPEN  (10,FILE='junky'//trim(adjustl(codec(i))),STATUS='unknown')
	READ (10,'(a10)') yminc
	READ (10,'(a10)') ymaxc
	READ (10,'(a10)') deltac 
	close(10)		
!
! ------ GMT options 
!
	r_option = & 
	"-R0/"//trim(adjustl(ninc))//"/"//trim(adjustl(yminc))//"/"//trim(adjustl(ymaxc))
	b_option = & 
	"-Ba2f1:'time (ka)':/a"//trim(adjustl(deltac))//"f"//trim(adjustl(deltac))//"WSen:'RSL (m)':"
!	
	Write(19,*) "psbasemap -U'SELEN 3.2' -X4 -Y14 ", b_option, r_option," -P -JX12.8/8  -K >  plot.ps" 
	Write(19,*) "psxy -H2 ", "rslp"//"-"//trim(adjustl(CODEC(i)))//".dat", " -B -R -JX -W4 -O -K >> plot.ps" 
	Write(19,*) "psxy -H2 ", "rsld"//"-"//trim(adjustl(CODEC(i)))//".dat", " -B -R -JX -Ss0.25 -G0 -O -K >> plot.ps" 
	if(RSL_FORMAT/='3')& 
	Write(19,*) "psxy -H2 ", "rsld"//"-"//trim(adjustl(CODEC(i)))//".dat", "  -Ey0.25/4 -B -R -JX -O -K >> plot.ps" 
	if(RSL_FORMAT=='3')& 
	Write(19,*) "psxy -H2 ", "rsld"//"-"//trim(adjustl(CODEC(i)))//".dat", " -Exy0.25/4 -B -R -JX -O -K >> plot.ps" 

!
! --- A new basemap for titles... 
	Write(19,*) "psbasemap -X0 -Y0 -Ba1000wesn -R0/18/0/10 -P -JX12.8/8 -O -K >>  plot.ps" 
!
	If(cd=='-1') then 
 	open (8,file='tmpb'//"."//trim(adjustl(CODEC(i))),status='unknown')
	Write(8,*) "8 9.0 9 0 0 BL  ", titre(i)			  	         
	Write(8,*) "8 8.5 9 0 0 BL  Ice model: ", trim(adjustl(titlice))
	Write(8,*) "8 8.0 9 0 0 BL  Repository: ./depot-", trim(adjustl(run))
	Write(8,*) "8 7.5 9 0 0 BL  Viscosity profile: ", trim(adjustl(VSTRING))
	Write(8,*) "8 7.0 8 0 0 BL      -LMAX =", trim(adjustl(DEGREE)),     & 
		       	                "-RES =", trim(adjustl(RESOLUTION)), & 
				      " -ITER =", trim(adjustl(ITER)),       & 
				      " -Mode =", trim(adjustl(MODE))	  
	Write(8,*) "8 6.5 8 0 0 BL -ALMA model =", trim(adjustl(SHORT_VISCO_FILENAME))
 	close(8) 	
	Else 
 	open (8,file='tmpb'//"."//trim(adjustl(CODEC(i))),status='unknown')
	Write(8,*) "8 9.0 9 0 0 BL  ", titre(i)			  	         
	Write(8,*) "8 8.5 9 0 0 BL  Ice model: ", trim(adjustl(titlice))
	Write(8,*) "8 8.0 9 0 0 BL  Repository label: ", trim(adjustl(run))
	Write(8,*) "8 7.5 9 0 0 BL  Viscosity profile: ", trim(adjustl(VSTRING))
	Write(8,*) "8 7.0 8 0 0 BL     -LMAX =", trim(adjustl(DEGREE)),     & 
		       	               "-RES =", trim(adjustl(RESOLUTION)), & 
				     " -ITER =", trim(adjustl(ITER)),       & 
				     " -Mode =", trim(adjustl(MODE))	  
	Write(8,*) "8 6.5 8 0 0 BL       -NV =", trim(adjustl(NV)),     & 
		       	             " -CODE =", trim(adjustl(CD))				     	         
 	close(8) 		
	Endif
!		
	Write(19,*) "pstext ", 'tmpb'//"."//trim(adjustl(CODEC(i))), " -O -JX -R >> plot.ps"	
!
	Write(19,*) "ps2pdf plot.ps"
!
	Write(19,*) "mv plot.pdf ", "rslp-"//trim(adjustl(CODEC(i)))//".pdf"
!
	Write(19,*) "mv plot.ps ", "rslp-"//trim(adjustl(CODEC(i)))//".ps"
!
	close(11) 
!
11111  CONTINUE
!
!Write(19,*) "echo       - RSL curves reported on files rslp-yyy.ps and rslp-yyy.pdf"		
	close(19)	
!	
 END SUBROUTINE MAKE_RSL
!
!
!
!
!
!
 SUBROUTINE MAKE_PALEO_TOPO_MAP (NINC,FILE)    		
 IMPLICIT NONE
!
! --- Creates a GMT script for plotting the paleotopography 
!     GS & FC July 29 2009

 CHARACTER*20 FILE, FNAME, OUT_FILE
 CHARACTER*2  LABCHAR, CHARJUNK
 CHARACTER*3  NINC 
 INTEGER I, N, NN
 CHARACTER*20 J_OPTION, B_OPTION, R_OPTION 
 CHARACTER*10 TITLICE  
 CHARACTER*10 C_TABLE_1, C_TABLE_2  
 CHARACTER*20 T_OPT_1, T_OPT_2,     & 
 	      N_TABLE_1, N_TABLE_2, & 
	      D_OPT_1, D_OPT_2  

!
!
!
!
 Call CHAR3_2_INT(NINC,NN)
!
 open(17,file=file,status='unknown')
!
!
! # Some settings ... 
!
 Write(17,*) "gmtset PAPER_MEDIA A4+"
 Write(17,*) "gmtset HEADER_FONT_SIZE 24p" 
 Write(17,*) "gmtset FRAME_WIDTH 0.1c"
 Write(17,*) "gmtset LABEL_FONT_SIZE 10p"
 Write(17,*) "gmtset ANOT_FONT_SIZE 10p"
 Write(17,*) "gmtset COLOR_BACKGROUND  -"
 Write(17,*) "gmtset COLOR_FOREGROUND -" 
 Write(17,*) "gmtset COLOR_NAN -"  
!
!
 Write(17,*) ""
!
 R_OPTION="-R0/360/-90/90"
 B_OPTION="-Ba180/a60WESn"
 J_OPTION="-JQ0/13"
 out_file="map.ps"
!
!
!---- a palette for the rock topography 
 c_table_1='-Crelief'
 n_table_1='pale_topo.cpt'
 t_opt_1  ='-T-10/10/1'
 d_opt_1="-D3.25/-1/6/0.6h"
 Write(17,*) "makecpt", " ", c_table_1, " ", t_opt_1, " > ", n_table_1 
 
!
!---- a palette for the ice topography 
 c_table_2='-Cgray'
 n_table_2='pale_ice.cpt'
 t_opt_2  ='-T0.0001/5.5/0.5'
 d_opt_2="-D9.75/-1/6/0.6h"

 Write(17,*) "makecpt ", " ", c_table_2, " -I ", t_opt_2, " > ", n_table_2 
!
!
!
 DO 1 I=0, NN 
! 
 Write(17,*) ""
 Write(17,*) ""
!
 If (i==0)    Write(17,*) "echo '     - until'", nn, "ka "
 If (i.eq.nn) Write(17,*) "echo '     - between 1 ka and today'"
!
 open(3,file='junk.dat',status='unknown') 
 if(i<=9) write(3,'(a1,i1)') '0',i  
 if(i> 9) write(3,'(i2)')        i
 close(3)
 open(3,file='junk.dat',status='unknown')
 read(3,'(a2)') labchar
 close(3)
!
 fname='ptopo-'//labchar//'.dat'
!
 Write(17,*) "awk '{print $1, $2, $3/1000}' ",  "ptopo-"//trim(adjustl(labchar))//".dat",  "> ", " ptopo.xyz"      
 Write(17,*) "pscontour -X3 -Y3 ptopo.xyz -I ", trim(adjustl(R_OPTION)), " ", trim(adjustl(J_OPTION)),       " ", & 
 				        trim(adjustl(B_OPTION)), " ", "-C"//trim(adjustl(n_table_1))," ", " -K > map.ps"
 Write(17,*) "psscale -Bf1a2:elevation'(km)': ", " ", trim(adjustl(d_opt_1)), " ", "-C"//trim(adjustl(n_table_1)), & 
             " -K -O >> map.ps"
!
 Write(17,*) ""
!
 Write(17,*) "awk '{print $1, $2, $4/1000}' ",  "imask-"//trim(adjustl(labchar))//".dat", "> ", " imask.xyz" 
 Write(17,*) "pscontour imask.xyz -I ", trim(adjustl(R_OPTION)), " ", trim(adjustl(J_OPTION)),       " ", & 
 				        trim(adjustl(B_OPTION)), " ", "-C"//trim(adjustl(n_table_2))," ", " -O -K >> map.ps" 
 Write(17,*) "psscale -Bf1a1:elevation'(km)': ", " ", trim(adjustl(d_opt_2)), " ", "-C"//trim(adjustl(n_table_2)), & 
             " -K -O >> map.ps"
 Write(17,*) "pscoast -Di -U/0.5/0.5/'SELEN 3.2' -R -A1000 ", trim(adjustl(J_OPTION)), " -W1/0/0/0  -O -K >> map.ps" 
!
! ------ Prepares a small title...

      open(8,file='ptmp-'//labchar//'.dat',status='unknown')
      if(i==0) then 
      Write(8,'(a15,1x,a23,i3,1x,a3)')    "0 99 18 0 1 BC ", "Paleo-topography until ", nn, "ka"
	       elseif (1<=i.and.i<=nn)     then       
      Write(8,'(a18,1x,a17,i3,a1,i3,a3)') "0 99 18 0 1 BC ", "Paleo-topography ", nn-i+1,"-",nn-i, "ka"
      endif
      close(8) 
      Write(17,*) "pstext -N ", "ptmp-"//labchar//".dat", " -O -JQ -R -G0 >> map.ps"



!
! ------ A pdf image and a name for each frame  
!
 Write(17,*) ""
	Write(17,*) "ps2pdf map.ps"
	Write(17,*) "mv map.pdf ", "map-topo-"//labchar//".pdf"
	Write(17,*) "mv map.ps ", "map-topo-"//labchar//".ps"
	Write(17,*) "/bin/rm ptopo.xyz imask.xyz"
!
!
1 CONTINUE 
!
 Close(17) 
!
!
 end subroutine make_paleo_topo_map 
!
!
!
!
!
!
 SUBROUTINE MAKE_ICEMAP(NINC, TITLICE, OPTION_ROF, FILECAP, FILE)    
 IMPLICIT NONE
 CHARACTER*1  OPTION_ROF
 CHARACTER*3  LABCHAR, CHARJUNK
 CHARACTER*3 NINC 
 CHARACTER*10 TITLICE  
 CHARACTER*20 FILE, FILECAP, FILENAME
 CHARACTER*20 C_TABLE, C_OPTION, B_OPTION 
 INTEGER I, N, NN
!
! --- Creates a GMT script for plotting the original ice sheets
!     distributions using various projections - 
!  GS November 14 2007
!  Revised February 21, 2008- 
!  Revised June 14, 2008. 
!  Revised (IJ05) July 5, 2008 
!  Also revised July 14, 2008 
!  Revised July 26, 2008 
!  Revised GS & GG May 2011 
!    *** Revised GS & GG June 2011 - Implementation of the "ANTA" ice model...
!    *** Revised GS July 2011 - Implementation of the "RECT" ice model...

 open(17,file=file,status='unknown')
!
!
! # Some settings ... 
!
 Write(17,*) "gmtset PAPER_MEDIA A4+"
 Write(17,*) "gmtset HEADER_FONT_SIZE 24p" 
 Write(17,*) "gmtset FRAME_WIDTH 0.1c" 
!
!
! ------ 
 C_TABLE="pice.cpt"
 C_OPTION=" -C"//trim(adjustl(C_TABLE))
 If(titlice=="IJ05MOD")then 
   Write(17,*) "makecpt -Crainbow  -T-100/1000/100 -D > ", trim(adjustl(C_TABLE)) 
   B_OPTION = " -Bf1000a100g100/:: "			      
 endif 
 If(titlice/="IJ05MOD")then 
   Write(17,*) "makecpt -Crainbow  -T-500/4500/500 -D > ", trim(adjustl(C_TABLE)) 
   B_OPTION = " -Bf4500a1000g500/:: "
 endif
!
 Write(17,*) ""
!
 open(10,file='junk.dat',status='unknown'); Write(10,'(a3)') NINC; close(10) 
 open(10,file='junk.dat',status='unknown'); Read(10,*) NN; close(10) 
!
! Write(17,*) "echo '--- ", trim(adjustl(file))//":", " Creating ps images of original ice sheets distributions'"
!
 do i=0, nn+1 
! 
 Write(17,*) ""
 If (i==0)  then 
 	Write(17,*) "echo '     - until'", nn, "ka "
!elseif (1<=i.and.i<=nn) then
!Write(17,*) "echo '     - Between '", nn-i+1, " and ", nn-i, "ka "	    
	    elseif (i.eq.nn+1) 	        then
 	Write(17,*) "echo '     - today'"
 Endif 
!
	open(3,file='junk.dat',status='unknown') 
!	if(i<=9) write(3,'(a2,i1)') '00',i  
!	if(i> 9) write(3,'(i3)')        i
	if(i<=9) write(3,'(a1,i1)') '0',i  
	if(i> 9) write(3,'(i2)')        i
 	close(3)
	open(3,file='junk.dat',status='unknown'); read(3,'(a3)') labchar; close(3)
!
 	IF(titlice=="ICE5G26")  filename='msg5-'//trim(adjustl(labchar))//'.dat'
 	IF(titlice=="ICE5G")    filename='msg5-'//trim(adjustl(labchar))//'.dat'
 	IF(titlice=="ICE3G")    filename='msg3-'//trim(adjustl(labchar))//'.dat'
 	IF(titlice=="IMED1")    filename='msgM-'//trim(adjustl(labchar))//'.dat'
 	IF(titlice=="DISK")     filename='msgD-'//trim(adjustl(labchar))//'.dat'
 	IF(titlice=="ICAP")     filename='msgC-'//trim(adjustl(labchar))//'.dat'
 	IF(titlice=="IJ05MOD")  filename='msgI-'//trim(adjustl(labchar))//'.dat'
 	IF(titlice=="ANU05")    filename='msgU-'//trim(adjustl(labchar))//'.dat'
 	IF(titlice=="ICE1")     filename='msg1-'//trim(adjustl(labchar))//'.dat'
 	IF(titlice=="ALPS")     filename='msgA-'//trim(adjustl(labchar))//'.dat'
 	IF(titlice=="FLOR")     filename='msgF-'//trim(adjustl(labchar))//'.dat'
 	IF(titlice=="GREEN")    filename='msgG-'//trim(adjustl(labchar))//'.dat'
 	IF(titlice=="GLAC")     filename='msgL-'//trim(adjustl(labchar))//'.dat'
 	IF(titlice=="ANTA")     filename='msgR-'//trim(adjustl(labchar))//'.dat'
 	IF(titlice=="TRANS")    filename='msgT-'//trim(adjustl(labchar))//'.dat'
 	IF(titlice=="RETT")     filename='msgY-'//trim(adjustl(labchar))//'.dat'
 	IF(titlice=="GRD")      filename='msgZ-'//trim(adjustl(labchar))//'.dat'
 	IF(titlice=="BENO")     filename='msgB-'//trim(adjustl(labchar))//'.dat'
	
!
!
! ------ Main frame (Mercator projection) 
!
	Write(17,*) "psbasemap -X4 -Y6 -Ba180/a85f90WSEn -R0/360/-85/85  -JM12 -U/4/13.5/'SELEN 3.2' -K > map.ps" 	
	Write(17,*) "psxy ", trim(adjustl(C_OPTION)), " -A -L -M ", trim(adjustl(filename)), " -O -K -JM -R >> map.ps" 
	If(option_rof=='r')then
		           Write(17,*) "pscoast -R -JM -Dc  -B -W1/0   -A10000 -O -K >> map.ps" 
		           else
		           Write(17,*) "pscoast -R -JM -Dc  -B -W1/240 -A10000 -O -K >> map.ps" 
	Endif			   
	Write(17,*) "psscale -E ", trim(adjustl(C_OPTION)), " ", trim(adjustl(B_OPTION)), " -D6/-1/11/1h -O  -K >> map.ps"
!
! ------ Prepares a small title 
!
 	open(8,file='rtmp0'//labchar//'.dat',status='unknown')
	if(i==0) then 
	Write(8,'(a18,1x,a5,a18,i2,1x,a8)')    "180 86.5 18 0 1 BC ", & 
	      				     trim(adjustl(titlice)), " thickness until ", nn, " ka "
 	      	 elseif (1<=i.and.i<=nn)     then	
	Write(8,'(a18,1x,a5,a11,i2,a1,i2,a8)') "180 86.5 18 0 1 BC ", & 
		   	                     trim(adjustl(titlice)), " thickness ", nn-i+1,"-",nn-i, " ka "
	         elseif (i==nn+1) 	     then
	Write(8,'(a18,1x,a5,a16)')             "180 86.5 18 0 1 BC ", & 
				             trim(adjustl(titlice)), " thickness today"
	endif
	Write(8,*) "180 -89.2 14 0 2 BC ", adjustl(trim(titlice)), " thickness (m)"		
	close(8) 
	Write(17,*) "pstext -N ", "rtmp0"//labchar//".dat", " -O -K -JM -R -G0 >> map.ps"

        if(option_rof=='z')& 
	Write(17,*) "psxy -R -JM ", trim(adjustl(filecap)), " -M -W1/0 -B -A -O >> map.ps"
!
!
	if(option_rof=='r') then 
!
! ------ Northern emisphere 
!	     	
 	Write(17,*) "psbasemap -Y6 -X14.5 -Ba90/a85f90WSEN -R0/360/0/90 -JG-45/90/8 -O -K >> map.ps"  
	Write(17,*) "psxy ", trim(adjustl(C_OPTION)), " -A -L -M ", trim(adjustl(filename)), " -O  -K -JG -R >> map.ps"
	If(option_rof=='r')then
		           Write(17,*) "pscoast -R -JG -Dc  -B -W1/0   -A10000 -O -K >> map.ps" 
		           else
		           Write(17,*) "pscoast -R -JG -Dc  -B -W1/240 -A10000 -O -K >> map.ps" 
	Endif	
!
! ------ Southern emisphere 
!
 	Write(17,*) "psbasemap -Y-9 -X0 -Ba90/a85f90WSEN -R0/360/-90/0  -JG-90/-90/8 -O -K  >> map.ps"  
	Write(17,*) "psxy ", trim(adjustl(C_OPTION)), " -A -L -M ", trim(adjustl(filename)), " -O  -K -JG -R >> map.ps"
	If(option_rof=='r')then
		           Write(17,*) "pscoast -R -JG -Dc  -B -W1/0   -A10000 -O >> map.ps" 
		           else
		           Write(17,*) "pscoast -R -JG -Dc  -B -W1/240 -A10000 -O >> map.ps" 
	Endif	
!
	Endif
!	
! ------ A pdf image and a name for each frame  
!
	Write(17,*) "ps2pdf map.ps"
	Write(17,*) "mv map.pdf ", "mapice"//labchar//".pdf"
	Write(17,*) "mv map.ps ", "mapice"//labchar//".ps"
!
enddo
!
 close(17)
!
 END SUBROUTINE MAKE_ICEMAP
!
!
!
!
!
!
 SUBROUTINE MAKE_RECICEMAP(NINC, DEGREE, TITLICE, OPTION_ROF, FILECAP, FILE)
 IMPLICIT NONE
 CHARACTER*1 OPTION_ROF
 CHARACTER*2 LABCHAR, CHARJUNK
 CHARACTER*3 NINC  
 CHARACTER*10 DEGREE, TITLICE   
 CHARACTER*20 FILE, FILECAP, B_OPTION, C_TABLE, C_OPTION 
 CHARACTER*12 FILENAME
 INTEGER I, N, NN
!
!
! --- Creates a GMT script for plotting the reconstructed ice sheets
!     distributions using various projections - GS November 07
! 
!     Revised July 05, 2008 - For INTEL version (2.6)  
!     Also revised July 14, 2008 
!     Revised July 26, 2008 
!
!
 open(1,file=file,status='unknown')
!
!
! # Some settings ... 
!
 Write(1,*) "gmtset PAPER_MEDIA A4+"
 Write(1,*) "gmtset HEADER_FONT_SIZE 24p" 
 Write(1,*) "gmtset FRAME_WIDTH 0.1c" 
!
 Write(1,*) ""
!
!
 open(10,file='junk.dat',status='unknown'); Write(10,'(a3)') NINC; close(10) 
 open(10,file='junk.dat',status='unknown'); Read(10,*) NN; close(10) 
!
! Write(1,*) "echo '--- ", trim(adjustl(file))//":", " Creating ps images of reconstructed ice sheets'"
!
 do i=0, nn+1 
! 
 Write(1,*) ""
 If (i==0)  then 
 	Write(1,*) "echo '     - until'", nn, "ka "
!elseif (1<=i.and.i<=nn) then
!Write(1,*) "echo '     - Between '", nn-i+1, " and ", nn-i, "ka "	    
	    elseif (i.eq.nn+1) 	        then
 	Write(1,*) "echo '     - today'"
 Endif 
!
	open(3,file='junk.dat',status='unknown') 
	if(i<=9) write(3,'(a1,i1)') '0',i  
	if(i> 9) write(3,'(i2)')        i
 	close(3)
	open(3,file='junk.dat',status='unknown'); read(3,'(a2)') labchar; close(3)
!
	filename='rect'//labchar//'.dat'
!
!
! ------ [NEW] Main frame (Mercator projection + PScontour) 
 C_TABLE="pice.cpt"
 C_OPTION=" -C"//trim(adjustl(C_TABLE))
 If(titlice=="IJ05MOD")then 
   Write(1,*) "makecpt -Crainbow  -T-100/1000/100 -D > ", trim(adjustl(C_TABLE)) 
   B_OPTION = " -Bf1000a100g100/:: "			      
 endif 
 If(titlice/="IJ05MOD")then 
   Write(1,*) "makecpt -Crainbow  -T-500/4500/500 -D > ", trim(adjustl(C_TABLE))  
   B_OPTION = " -Bf4500a1000g500/:: "
 endif
!
        Write(1,*) "pscontour -X3 -Y5  -I -JM12 -Ba180/a85f90WSEn -R0/360/-85/85 ", filename, " ", & 
	            trim(adjustl(C_OPTION)), " -K > map.ps"
!
	If(option_rof=='r')then
		   Write(1,*) "pscoast -U/4/13.5/'SELEN 3.2' -B -R -O -K -W2/255 -JM -Dc -A10000 >> map.ps" 
		           else
		   Write(1,*) "pscoast -U/4/13.5/'SELEN 3.2' -B -R -O -K -W2/0/0/255 -JM -Dc -A10000 >> map.ps" 		   
	Endif
!
        if(option_rof=='z')& 
	Write(1,*) "psxy -R -JM ", trim(adjustl(filecap)), " -M -W1/0 -B -A -O -K >> map.ps"
!        
        Write(1,*) "psscale -E ", trim(adjustl(C_OPTION)), " ", trim(adjustl(B_OPTION)), & 
	           " -D6/-1/11/1h -O -K >> map.ps"
        Write(1,*) "pstext -N ", " tmp0"//labchar//".dat",  " -JM -R -G0 -O -K >> map.ps"
!
! ------ Prepares a small title 
!
 	open(8,file='tmp0'//labchar//'.dat',status='unknown')
	if(i==0) then 
	Write(8,'(a18,1x,a7,a21,a4,a3,a7,i2,a9)')    "180 86.5 18 0 1 BC  ", trim(adjustl(titlice)), & 
		                                  " thickness to degree ", trim(adjustl(degree)), & 
					          " / ", " until ", nn, " ka  "
 	      	 elseif (1<=i.and.i<=nn) then	
	Write(8,'(a18,1x,a7,a21,a4,a3,i2,a1,i2,a9)') "180 86.5 18 0 1 BC  ", trim(adjustl(titlice)), & 
	        			          " thickness to degree ", trim(adjustl(degree)), & 
					          " / ", nn-i+1,"-",nn-i, " ka  "
	         elseif (i.eq.nn+1) then
	Write(8,'(a20,1x,a7,a21,a4,a3,a6)')          "180 86.5 18 0 1 BC  ", trim(adjustl(titlice)), & 
	        				  " thickness to degree ", trim(adjustl(degree)), & 
						  " / ", " today"
	endif	
	Write(8,*) "180 -89.2 14 0 2 BC ", trim(adjustl(titlice)), " thickness (m)"		
	close(8) 
!
!--- Southern emisphere 
!
 	Write(1,*) "pscontour -I -Y1 -X14.5 -JA0/-90/10 -R0/360/-90/-50 ", filename, " ", & 
	            trim(adjustl(C_OPTION)), " -O -K >> map.ps"
! 	
	If(option_rof=='r')then
		   Write(1,*) "pscoast -R -JA -Di -Ba45f45 -W2/255     -A0 -O >> map.ps"
		           else
		   Write(1,*) "pscoast -R -JA -Di -Ba45f45 -W2/0/0/255 -A0 -O >> map.ps"		   
	Endif
!
!
! ------ A pdf image and a name for each frame  
!
	Write(1,*) "ps2pdf map.ps"
	Write(1,*) "mv map.pdf ", "recice"//labchar//".pdf"
	Write(1,*) "mv map.ps ",  "recice"//labchar//".ps"
!
enddo
!
 close(1)
!
 END SUBROUTINE MAKE_RECICEMAP
!
!
!
!
!
!
 SUBROUTINE MAKE_OFMAP(DEG, OPT, AMP, FILECAP, FILE)
 IMPLICIT NONE
 CHARACTER*1 OPT
 CHARACTER*10 DEG
 CHARACTER*10 AMP 
 CHARACTER*20 FILE, FILECAP 
!
! Prepares a GMT script for plotting the reconstructed and 
! the original Ocean Function - 
! ---> Last changes GS Nov 14 2007
! ---> Revised June 25, 2008 for implementation of "Selen 2.6"
! ---> Revised July XX, 2008 for implementation of the "Zonal" OF
!
 open(1,file=file,status='unknown')
!
 If(opt=='r')then 
!
! --- Realistic OF 
!
   Write(1,*)"gmtset PAPER_MEDIA A4+"
   Write(1,*)"gmtset HEADER_FONT_SIZE 24p"
   Write(1,*)"gmtset FRAME_WIDTH 0.1c"
   Write(1,*)"gmtset ANOT_FONT_SIZE 12p"
   Write(1,*)" "
!
   Write(1,*)"# Map of the original ocean function"
   Write(1,*)"psbasemap -X3 -Y18 -P -Bf180a90/f90a60WSEN -R0/360/-85/85 -JQ180/16 -K > of.ps"
   Write(1,*)"pscoast -R -JQ -Di -B -W1/255/255/255 -A1000 -S255/0/0 -G0/0/255 -O -K >> of.ps"
   Write(1,*)"echo '180 120 18 0 2 CM OCEAN FUNCTION 0 ka' | pstext -U/0.5/0.5/'SELEN 3.2' -N -R -JQ -O -K >> of.ps"
   Write(1,*)" "
!
   Write(1,*)"# Map of reconstructed ocean function"
   Write(1,*)"makecpt  -Cno_green -T0/1/0.1 > pale.cpt"
   Write(1,*)"psbasemap -X0 -Y-11 -P -B -R -JQ -K -O >> of.ps"
   Write(1,*)"pscontour -I -JQ -R recof.dat -Cpale.cpt -O -K >> of.ps"
   Write(1,*)"pscoast -R -JQ -Di -W1/255/255/255 -B -A1000 -O -K >> of.ps"
   Write(1,*)"echo '180 120 18 0 2 CM OCEAN FUNCTION 0 ka to degree ", trim(adjustl(deg)),"'", & 
   " | pstext -U/0.5/0.5/'SELEN 3.2' -N -R -JQ -O -K >> of.ps"
   Write(1,*)" "
!
   Write(1,*)"# A color table"
   Write(1,*)"psscale -E -Cpale.cpt -B1f0.1a0.5g0.5/:O.F.: -D8/-2/8/1h -O >> of.ps"
!
 Else
!
! --- "Zonal" OF 
!
   Write(1,*)"gmtset PAPER_MEDIA A4+"
   Write(1,*)"gmtset HEADER_FONT_SIZE 24p"
   Write(1,*)"gmtset FRAME_WIDTH 0.1c"
   Write(1,*)"gmtset ANOT_FONT_SIZE 12p"
   Write(1,*)" "
!
   Write(1,*)"# Map of the original ocean function"
   Write(1,*)"psbasemap -X3 -Y18 -P -Bf180a90/f90a60WSEN -Gred -R0/360/-85/85 -JQ180/16 -K > of.ps"
   Write(1,*)"psxy -R -JQ ", trim(adjustl(FILECAP)), " -M -W1/255 -Gblue -L -B -A -O -K >> of.ps"
   Write(1,*)"echo '180 120 18 0 2 CM OCEAN FUNCTION' | pstext -U/0.5/0.5/'SELEN 3.2' -N -R -JQ -O -K >> of.ps"
   Write(1,*)" "
!
   Write(1,*)"# Map of reconstructed ocean function"
   Write(1,*)"makecpt  -Cno_green -T0/1/0.1 > pale.cpt"
   Write(1,*)"psbasemap -X0 -Y-11 -P -B -R -JQ -K -O >> of.ps"
   Write(1,*)"pscontour -I -JQ -R recof.dat -Cpale.cpt -O -K >> of.ps" 
   Write(1,*)"psxy -R -JQ ", trim(adjustl(FILECAP)), " -M -W1/255 -B -A -O -K >> of.ps"
   Write(1,*)"echo '180 120 18 0 2 CM OCEAN FUNCTION to degree ", trim(adjustl(deg)),"'", & 
   " | pstext -U/0.5/0.5/'SELEN 3.2' -N -R -JQ -O -K >> of.ps"
   Write(1,*)" "
!
   Write(1,*)"# A color table"
   Write(1,*)"psscale -E -Cpale.cpt -B1f0.1a0.5g0.5/:O.F.: -D8/-2/8/1h -O >> of.ps"
!
 Endif
!
 Close(1) 
!
 END SUBROUTINE MAKE_OFMAP
!
!
!
!
!
!
 SUBROUTINE MAKE_PXMAP(NRES, FILE)
 IMPLICIT NONE
 INTEGER NP, NRES
 CHARACTER*20 FILE
 CHARACTER*20, PARAMETER :: SIZE_OF_PIXELS="-Sc0.04" 
! 
!
! Prepares a GMT script for plotting  wet, dry, and global 
! pixels distribution  ---- Last changes GS Feb 21 2008 ----
!
! Last reviewed July 1, 2008 - GS for port to SELEN 2.6 
! Now the script generated TWO plots: 1) wet/dry distributions
! and 2) a global "spherical" view of pixels - 
! 
 NP=2*nres*(nres-1)*20+12
!
Open(11,file=file,status='unknown') 

Write(11,*)"gmtset PAPER_MEDIA A4+" 
Write(11,*)"gmtset HEADER_FONT_SIZE 24p"
Write(11,*)"gmtset FRAME_WIDTH 0.1c"
Write(11,*)"gmtset ANOT_FONT_SIZE 12p"
!
Write(11,*)"echo '     - wet pixels'"
Write(11,*)""
Write(11,*)"# Map of wet pixels distribution"	
Write(11,*)""	
Write(11,*)"psbasemap -Y18 -P -Bf180a90/f90a60WSEN -R0/360/-85/85 -JQ180/16 -K > px.ps"
Write(11,*)"pscoast -B -R -O -K -W1/255/0/0 -JQ -Di -A1000 >> px.ps"
Write(11,*)"psxy -G255/0/0 -H4 weta.dat -O ", trim(adjustl(size_of_pixels)), " -JQ -R -K >> px.ps"
Write(11,*)"echo '180 120 18 0 0 CM WET pixels' | pstext -N -R -JQ -O -K >> px.ps"
Write(11,'(a27,i3,a4,i6,a1,a52)') "echo '10 -70 14 0 2 LM RES=", nres, ", N=", NP, "'", & 
" | pstext -N -R -JQ -O -K -G255/0/0 -W255 >> px.ps"
!
Write(11,*)"echo '     - dry pixels'"
Write(11,*)""
Write(11,*)"# Map of dry pixels distribution"	
Write(11,*)""	
Write(11,*)"psbasemap -X0 -Y-11 -P -Bf180a90/f90a60WSEN -R0/360/-85/85 -JQ180/16 -K -O >>  px.ps"
Write(11,*)"pscoast -B -R -O -K -W1/0/0/255 -JQ -Di -A1000 >> px.ps"
Write(11,*)"psxy -G0/0/255 -H4 drya.dat -O ", trim(adjustl(size_of_pixels)), " -JQ -R -K >> px.ps"
Write(11,*)"psbasemap -U/6/-2/'SELEN 3.2' -B -R -P -JQ -K -O >> px.ps"
Write(11,*)"echo '180 120 18 0 0 CM DRY pixels' | pstext -N -R -JQ -O -K >> px.ps"
Write(11,'(a27,i3,a4,i6,a1,a52)') "echo '10 -70 14 0 2 LM RES=", nres, ", N=", NP, "'", & 
" | pstext -N -R -JQ -O -G0/0/255 -W255 >> px.ps"
!
Write(11,*)"echo '     - Spherical map of wet and dry pixels'"
Write(11,*)""
Write(11,*)"# Spherical map of dry & wet pixels distribution"	
Write(11,*)""
Write(11,*)"psbasemap -X4 -Y12 -U/0/-3/'SELEN 3.2' -Bg60f60/g60f60 -R0/360/-80/80 -P -JG-20/-44/14 -K > px-sphere.ps"
Write(11,*)"pscoast -B -R -O -K -W1/220/220/220 -JG -Di -S0/60/255 -G100/250/100 -A1000 >> px-sphere.ps"
Write(11,*)"psxy  -G255 -H4 pxa.dat -O -K ", trim(adjustl(size_of_pixels)),  " -JG -R >> px-sphere.ps"
Write(11,'(a27,i3,a4,i6,a1,a49)') "echo '0 -110 14 0 2 LM RES=", nres, ", N=", NP, "'", & 
" | pstext -N -R -JQ -O -G0 -W255 >> px-sphere.ps"
!
 close(11)
! 
  end subroutine make_pxmap
!
!
!
!
!
!
 SUBROUTINE MAKE_RSLSCA (NV, CODE, RSL_FILE, RUN, NINC, NRSL, TITLICE, RESOLUTION, & 
			          ITER, MODE, DEGREE, VSTRING, FILE_GMT, VISCO_FILENAME)
!
! Prepares a simple GMT script for drawing a scatterplot of 
! global RSL predictions vs. observations  - GS 31.10-.2007 
!
! Deeply Revised on July 2008 for the implementation of SELEN 2.6
! Also revised on April 2010 GS 
!
 IMPLICIT NONE
 INTEGER, PARAMETER :: MAXDATA=100000
 REAL TIME, RSL, DMIN, DMAX, TMIN, TMAX, RANGE, DELTA 
 CHARACTER*100 VSTRING, VISCO_FILENAME
 CHARACTER*20 FILE_GMT, YMINC, YMAXC, TMINC, TMAXC, DELTAC, R_OPTION
 CHARACTER*30 RSL_FILE 
 CHARACTER*10 RESOLUTION, DEGREE, TITLICE, STRING 
 CHARACTER*50  B_OPTION
 CHARACTER*3 RUN
 CHARACTER*3 NV, CODE 
 CHARACTER*3 NINC
 CHARACTER*1 ITER, MODE
 INTEGER I, J, K, NRSL 
!
!
! ==== Creates a GMT script for a scatterplot of RSL data 
! 	 
! --- Determines min and max of cumulative RSL data 
open(11,file='scatter-data.dat',status='unknown') 
dmax=-99999.
dmin=+99999.
tmax=-99999.
tmin=+99999.
do j=1, maxdata
       read(11,*,end=88) time, rsl 
       if(rsl<=dmin)  dmin=rsl 
       if(rsl>=dmax)  dmax=rsl 
       if(time<=tmin) tmin=time
       if(time>=tmax) tmax=time        
enddo
88 continue
   close(11)  
!
	range=dmax-dmin  
!	
if(range>=0  .and.range<=50) delta=10  
if(range>=50 .and.range<100) delta=10	
if(range>=100.and.range<200) delta=20
if(range>=200.and.range<300) delta=50
if(range>=300.and.range<400) delta=100
if(range>400) 		     delta=150 
!
OPEN  (10,FILE='junk.dat',STATUS='unknown')
WRITE (10,'(i10)') int(dmin-delta/2.)
WRITE (10,'(i10)') int(dmax+delta/2.)
WRITE (10,'(i10)') int(delta) 
WRITE (10,'(i10)') int(tmin)-1		
WRITE (10,'(i10)') int(tmax)+1	
 CLOSE(10)
OPEN  (10,FILE='junk.dat',STATUS='unknown')
READ (10,'(a10)') yminc
READ (10,'(a10)') ymaxc
READ (10,'(a10)') deltac
READ (10,'(a10)') tminc 
READ (10,'(a10)') tmaxc 	
 CLOSE(10)
!
! --- Target GMT file for scatterplot of RSL data 
!   
 open (19,file=file_gmt,status='unknown')
 Write(19,*)"gmtset PAPER_MEDIA A4+"
 Write(19,*)"gmtset HEADER_FONT_SIZE 24p"
 Write(19,*)"gmtset FRAME_WIDTH 0.1c"
 Write(19,*)"gmtset ANOT_FONT_SIZE 16p"
 Write(19,*)"gmtset LABEL_FONT_SIZE 16p"
!
 r_option = & 
 "-R"//trim(adjustl(tminc))//"/"//trim(adjustl(tmaxc))//"/"//trim(adjustl(yminc))//"/"//trim(adjustl(ymaxc))
!
 b_option = & 
 "-Ba2f1/a"//trim(adjustl(deltac))//"f"//trim(adjustl(deltac))//"WSen:'RSL (m)':"
!
 Write(19,*)"psbasemap -X4.5 -Y18 ", b_option, r_option," -P -JX12/7 -K > plot.ps" 
 Write(19,*)"psxy scatter-data.dat ", " -B -R -JX -Ss0.1 -G0 -O -K >> plot.ps" 
 Write(19,*)"psxy scatter-data.dat ", "-Ey0.25/2 -B -R -JX -O -K >> plot.ps" 
 Write(19,*)"echo '", int(tmin)+(int(tmax)-int(tmin))/2., (11./10.)*int(dmax+delta/2.), " 22 0 0 BC RSL data from file: ", & 
 trim(adjustl(RSL_FILE)), "'", " | pstext -JX -R -O -K -N >> plot.ps "  
!
 b_option = & 
 "-Ba2f1:'time (ka)':/a"//trim(adjustl(deltac))//"f"//trim(adjustl(deltac))//"WSen:'RSL (m)':"
!
 Write(19,*)"psbasemap -U/-2/-3/'SELEN 3.2' -Y-11 ", b_option, r_option, " -P -JX -O -K >> plot.ps"
 Write(19,*)"psxy scatter-pred.dat -B -R -JX -Sc0.12 -G0 -O -K >> plot.ps"
 Write(19,*)"echo '", int(tmin)+(int(tmax)-int(tmin))/2., (13./10.)*int(dmax+delta/2.), " 22 0 2 BC RSL predictions for: ", & 
 trim(adjustl(TITLICE)), "'", " | pstext -JX -R -O -K -N >> plot.ps "   
!
 open (8,file='title.tmp',status='unknown')  
 If(code=='-1')then 
 Write(8,*) int(tmin)+(int(tmax)-int(tmin))/2., (11./10.)*int(dmax+delta/2.), " 13 0 0 BC ",& 
    " -LMAX=", trim(adjustl(DEGREE)),& 
    " -RES=", trim(adjustl(RESOLUTION)),& 
    " -MODE=", trim(adjustl(MODE)),&   
    " -ITER=", trim(adjustl(ITER)),& 
    " -ALMA model=", trim(adjustl(VISCO_FILENAME))
 close(8) 
 Else 
 Write(8,*) int(tmin)+(int(tmax)-int(tmin))/2., (11./10.)*int(dmax+delta/2.), " 13 0 0 BC ",& 
    " -LMAX=", trim(adjustl(DEGREE)),& 
    " -RES=", trim(adjustl(RESOLUTION)),& 
    " -NV=", trim(adjustl(NV)),& 
    " -CODE=", trim(adjustl(CODE)),& 
    " -MODE=", trim(adjustl(MODE)),&   
    " -ITER=", trim(adjustl(ITER)),& 
    "  -Visco: ", trim(adjustl(VSTRING))  
 close(8) 
 Endif 
 Write(19,*)"pstext title.tmp -N -JX -R -O >> plot.ps"
 Write(19,*)"mv plot.ps scatter-plot.ps"
!
 Close(19) 
!
 END SUBROUTINE MAKE_RSLSCA 
!
!
!
!
!
!
 SUBROUTINE MAKE_PM (RUN,      & 
                     NV,       & 
		     CODE,     & 
		     TITLICE,  & 
		     RESOLUTION, & 
 	             ITER,     & 
		     MODE,     & 
		     DEGREE,   & 
		     VSTRING,  & 
		     SHORT,    & 
		     NINC,     & 
		     FILE_GMT)
!
! -----------------------------------------------------------------
! --- Creates a GMT script for plotting the polar motion and the 
!     rate of polar motion as a function of time  -GS Aug 15 2010
! -----------------------------------------------------------------
!
 IMPLICIT NONE
! 
 CHARACTER*10   TITLICE,  DEGREE, RESOLUTION     
 CHARACTER*20   R_OPTION, G_OPTION, W_OPTION, H_OPTION, J_OPTION 
 CHARACTER*200  B_OPTION 
 CHARACTER*20   FILE_GMT   
 CHARACTER*30   VSTRING
 CHARACTER*60   TITRE 
 CHARACTER*3    RUN, STRING
 CHARACTER*3    NV, NINC, CODE
 CHARACTER*1    ITER, MODE 
 CHARACTER*100  SHORT
 CHARACTER*7    MIDDLE
!
 open(9,file=file_gmt,status='unknown')
!
 Write(9,*) "" 
 Write(9,*) "# -------------------------------------- "
 Write(9,*) "# ------ Settings for BOTH plots  ------ " 		      
 Write(9,*) "# -------------------------------------- "
 Write(9,*) "" 
!
 Write(9,*) "gmtset PAPER_MEDIA     A4+"
 Write(9,*) "gmtset ANOT_FONT_SIZE  16p"
 Write(9,*) "gmtset LABEL_FONT_SIZE 16p"
 Write(9,*) "gmtset TICK_LENGTH   -0.2c"
!
! --- Header lines in files "m.dat" and "m.dot" 
 H_OPTION = " -H7 " 
!
! --- Header lines in files "m.dat" and "m.dot" 
 J_OPTION = " -JX12/8 " 
!
!
!
!
 Write(9,*) "" 
 Write(9,*) "# ----------------------------------------- "
 Write(9,*) "# ------ Polar motion since the LGM  ------ " 		      
 Write(9,*) "# ----------------------------------------- "
!
! --- Recommended y-range for plot: the range for polar motion is [-0.1,0.1] degrees 
!
                R_OPTION=" -R0/20.0/-0.1/0.1 "      
 if(ninc=='18') R_OPTION=" -R0/18.0/-0.1/0.1 "      
 if(ninc=='21') R_OPTION=" -R0/10.5/-0.1/0.1 "      
 if(ninc=='26') R_OPTION=" -R0/13.0/-0.1/0.1 "      
! 
                MIDDLE= " 10.00 "
 if(ninc=='18') MIDDLE= "  9.00 "      
 if(ninc=='21') MIDDLE= " 10.50 "      
 if(ninc=='26') MIDDLE= " 13.00 "     
!
! --- B-option 
 B_OPTION="-Ba2f1:'time since LGM (ka)':/f0.025a0.05:'polar motion (deg)'"//":WSne"
! 
 Write(9,*) "" 
 Write(9,*) "#--- Base of plot for polar motion"
 Write(9,*) "psbasemap -X6 -Y4 -U'SELEN 3.2' ", trim(adjustl(R_OPTION)), " ", & 
                                                trim(adjustl(B_OPTION)), " ", & 
						trim(adjustl(J_OPTION)), " -K > m.ps "
 Write(9,*) "" 
 Write(9,*) "#--- Extracts polar motion data (x, y, and modulus)"
 Write(9,*) "awk '{print $1, $3}' m.dat >   mx.tmp  " 
 Write(9,*) "awk '{print $1, $4}' m.dat >   my.tmp  "
 Write(9,*) "awk '{print $1, $6}' m.dat >  mod.tmp "
!
 Write(9,*) "" 
 Write(9,*) "#--- x component of polar motion"
 Write(9,*) "psxy mx.tmp -B -JX -R -W3/.. ", trim(adjustl(H_OPTION)), " -K -O >> m.ps"
 Write(9,*) "" 
 Write(9,*) "#--- x legend"
 Write(9,*) "psxy -R -JX -W3/.. -K -O  <<END  >> m.ps"
 Write(9,'(A9)') "1  -0.050" 
 Write(9,'(A9)') "3  -0.050" 
 Write(9,'(A3)') "END"
 Write(9,*) "echo '3.5 -0.050 10 0 2 LM x' | pstext  -N -R -JX -O -K >> m.ps"  
!
 Write(9,*) "" 
 Write(9,*) "#--- y component of polar motion"
 Write(9,*) "psxy my.tmp -B -JX -R -W3/- ", trim(adjustl(H_OPTION)), " -K -O >> m.ps"
 Write(9,*) "" 
 Write(9,*) "#--- x legend"
 Write(9,*) "psxy -R -JX -W3/- -K -O  <<END  >> m.ps"
 Write(9,'(A9)') "1  -0.065" 
 Write(9,'(A9)') "3  -0.065" 
 Write(9,'(A3)') "END"
 Write(9,*) "echo '3.5 -0.065 10 0 2 LM y' | pstext  -N -R -JX -O -K >> m.ps"  
!
 Write(9,*) "" 
 Write(9,*) "#--- Modulus of polar motion"
 Write(9,*) "psxy mod.tmp -B -JX -R -W6 ", trim(adjustl(H_OPTION)), " -K -O >> m.ps"
 Write(9,*) "" 
 Write(9,*) "#--- legend"
 Write(9,*) "psxy -R -JX -W6 -K -O  <<END  >> m.ps"
 Write(9,'(A9)') "1  -0.080" 
 Write(9,'(A9)') "3  -0.080" 
 Write(9,'(A3)') "END"
 Write(9,*) "echo '3.5 -0.080 10 0 2 LM Mod' | pstext  -N -R -JX -O -K >> m.ps"  
!
 Write(9,*) "" 
 Write(9,*) " ##### Zero line" 
 Write(9,*) " psxy -B -R -JX -O -W1/... -K <<END >> m.ps " 
 Write(9,'(A4)')"-2 0" 
 Write(9,'(A5)')"100 0" 
 Write(9,'(A3)')"END"
!
 If(code/='-1') then 
 open(4,file='pm-1.tmp',status='unknown') 
     Write(4,*) MIDDLE, "0.130 16 0 0 BC -Ice model: ", trim(adjustl(TITLICE)), & 
                             " -Viscosity profile: ", trim(adjustl(VSTRING)) 
     Write(4,*) MIDDLE, "0.115 11 0 2 BC ", &
                " -LMAX=", trim(adjustl(DEGREE)),  " -RES=",  trim(adjustl(RESOLUTION)),& 	  
                " -NV=",   trim(adjustl(NV)),      " -CODE=", trim(adjustl(CODE)),& 
                " -MODE=", trim(adjustl(MODE)),    " -ITER=", trim(adjustl(ITER))    
 close(4) 
 else
 open(4,file='pm-1.tmp',status='unknown') 
     Write(4,*) MIDDLE, "0.130 16 0 0 BC -Ice model: ", trim(adjustl(TITLICE)), & 
                                  " -ALMA rheology:", trim(adjustl(SHORT)) 
     Write(4,*) MIDDLE, "0.115 11 0 2 BC ", &
                " -LMAX=", trim(adjustl(DEGREE)),  " -RES=",  trim(adjustl(RESOLUTION)),& 	  
                " -MODE=", trim(adjustl(MODE)),    " -ITER=", trim(adjustl(ITER))    
 close(4) 
 endif
!
 Write(9,*) " "
 Write(9,*) "pstext pm-1.tmp -N -JX -R -G0 -O -K >> m.ps"
!
 Write(9,*) "psbasemap -J -R -B::/::EWSN -O >> m.ps " 
!
!
!
!
!
 Write(9,*) "" 
 Write(9,*) "# ------------------------------------------------- "
 Write(9,*) "# ------ Rate of polar motion since the LGM  ------ " 		      
 Write(9,*) "# ------------------------------------------------- "
!
!
! --- Recommended y-range for plot: the range for polar motion is [-0.1,0.1] degrees 
!
                R_OPTION=" -R0/20.0/-7.5/7.5 "      
 if(ninc=='18') R_OPTION=" -R0/18.0/-7.5/7.5 "      
 if(ninc=='21') R_OPTION=" -R0/10.5/-7.5/7.5 "      
 if(ninc=='26') R_OPTION=" -R0/13.0/-7.5/7.5 "      
! 
                MIDDLE= " 10.00 "
 if(ninc=='18') MIDDLE= "  9.00 "      
 if(ninc=='21') MIDDLE= " 10.50 "      
 if(ninc=='26') MIDDLE= " 13.00 "     
!
! --- B-option 
 B_OPTION="-Ba2f1:'time since LGM (ka)':/f0.5a2:'rate of polar motion (deg/Ma)'"//":WSne"
! 
 Write(9,*) "" 
 Write(9,*) "#--- Base of plot for polar motion"
 Write(9,*) "psbasemap -X6 -Y4 -U'SELEN 3.2' ", trim(adjustl(R_OPTION)), " ", & 
                                                trim(adjustl(B_OPTION)), " ", & 
						trim(adjustl(J_OPTION)), " -K > m-dot.ps "
 Write(9,*) "" 
 Write(9,*) "#--- Extracts polar motion data (x, y, and modulus)"
 Write(9,*) "awk '{print $1, $3}' m.dot >   mdot-x.tmp  " 
 Write(9,*) "awk '{print $1, $4}' m.dot >   mdot-y.tmp  "
 Write(9,*) "awk '{print $1, $6}' m.dot > mdot-mod.tmp  "
!
 Write(9,*) "" 
 Write(9,*) "#--- x component of polar motion"
 Write(9,*) "psxy mdot-x.tmp -B -JX -R -W3/.. ", trim(adjustl(H_OPTION)), " -K -O >> m-dot.ps"
 Write(9,*) "" 
 Write(9,*) "#--- x legend"
 Write(9,*) "psxy -R -JX -W3/.. -K -O  <<END  >> m-dot.ps"
 Write(9,'(A9)') "1  -3" 
 Write(9,'(A9)') "3  -3" 
 Write(9,'(A3)') "END"
 Write(9,*) "echo '3.5 -3 10 0 2 LM x' | pstext  -N -R -JX -O -K >> m-dot.ps"  
!
 Write(9,*) "" 
 Write(9,*) "#--- y component of polar motion"
 Write(9,*) "psxy mdot-y.tmp -B -JX -R -W3/- ", trim(adjustl(H_OPTION)), " -K -O >> m-dot.ps"
 Write(9,*) "" 
 Write(9,*) "#--- x legend"
 Write(9,*) "psxy -R -JX -W3/- -K -O  <<END  >> m-dot.ps"
 Write(9,'(A9)') "1  -4" 
 Write(9,'(A9)') "3  -4" 
 Write(9,'(A3)') "END"
 Write(9,*) "echo '3.5 -4 10 0 2 LM y' | pstext  -N -R -JX -O -K >> m-dot.ps"  
!
 Write(9,*) "" 
 Write(9,*) "#--- Modulus of polar motion"
 Write(9,*) "psxy mdot-mod.tmp -B -JX -R -W6 ", trim(adjustl(H_OPTION)), " -K -O >> m-dot.ps"
 Write(9,*) "" 
 Write(9,*) "#--- legend"
 Write(9,*) "psxy -R -JX -W6 -K -O  <<END  >> m-dot.ps"
 Write(9,'(A9)') "1  -5" 
 Write(9,'(A9)') "3  -5" 
 Write(9,'(A3)') "END"
 Write(9,*) "echo '3.5 -5 10 0 2 LM Mod' | pstext  -N -R -JX -O -K >> m-dot.ps"  
!
 Write(9,*) "" 
 Write(9,*) " ##### Zero line" 
 Write(9,*) " psxy -B -R -JX -O -W1/... -K <<END >> m-dot.ps " 
 Write(9,'(A4)')"-2 0" 
 Write(9,'(A5)')"100 0" 
 Write(9,'(A3)')"END"
!
 If(code/='-1') then 
 open(4,file='pm-2.tmp',status='unknown') 
     Write(4,*) MIDDLE, "9.50 16 0 0 BC -Ice model: ", trim(adjustl(TITLICE)), & 
                             " -Viscosity profile: ", trim(adjustl(VSTRING)) 
     Write(4,*) MIDDLE, "8.25 11 0 2 BC ", &
                " -LMAX=", trim(adjustl(DEGREE)),  " -RES=",  trim(adjustl(RESOLUTION)),& 	  
                " -NV=",   trim(adjustl(NV)),      " -CODE=", trim(adjustl(CODE)),& 
                " -MODE=", trim(adjustl(MODE)),    " -ITER=", trim(adjustl(ITER))    
 close(4) 
 else
 open(4,file='pm-2.tmp',status='unknown') 
     Write(4,*) MIDDLE, "9.50 16 0 0 BC -Ice model: ", trim(adjustl(TITLICE)), & 
                                  " -ALMA rheology:", trim(adjustl(SHORT)) 
     Write(4,*) MIDDLE, "8.25 11 0 2 BC ", &
                " -LMAX=", trim(adjustl(DEGREE)),  " -RES=",  trim(adjustl(RESOLUTION)),& 	  
                " -MODE=", trim(adjustl(MODE)),    " -ITER=", trim(adjustl(ITER))    
 close(4) 
 endif
!
 Write(9,*) " "
 Write(9,*) "pstext pm-2.tmp -N -JX -R -G0 -O -K >> m-dot.ps"
!
 Write(9,*) "psbasemap -J -R -B::/::EWSN -O >> m-dot.ps " 

 
 
 
 
 
 
 
 
!
 CLOSE(9) 
!
 END SUBROUTINE MAKE_PM

!
!
!
!
!
!
 SUBROUTINE MAKE_STOKES (RUN, NV, CODE, TITLICE, RESOLUTION, & 
 			ITER, MODE, DEGREE, VSTRING, SHORT, FILE_GMT)
!
! --- Creates a GMT script for plotting the rates - of - variation
!     of cosine and sine stokes coeffciients -GS November 12 2007
!
!     Revised JULY 2008 for version 2.6 of SELEN - 
!     Revised Luky 2010 - Fully normalized coefficients---
!
 IMPLICIT NONE
 CHARACTER*10 TITLICE,  DEGREE, RESOLUTION     
 CHARACTER*20 R_OPTION, G_OPTION, W_OPTION, H_OPTION, FILE_GMT 
 CHARACTER*150 B_OPTION
 CHARACTER*30 VSTRING
 CHARACTER*60 TITRE 
 CHARACTER*3 RUN, STRING
 CHARACTER*3 NV, CODE
 CHARACTER*1 ITER, MODE 
 CHARACTER*100 SHORT
 INTEGER L, M, J_INDEX 
!
 open(9,file=file_gmt,status='unknown')
!
 Write(9,*) "gmtset PAPER_MEDIA A4+"
 Write(9,*) "gmtset ANOT_FONT_SIZE 16p"
 Write(9,*) "gmtset LABEL_FONT_SIZE 20p"
 Write(9,*) " "
 Write(9,*) "# ------ A plot of dot(c_lm,s_lm) vs degree (at present time)  ----" 		      
 Write(9,*) " "
!
! --- Recommended y-range for plot 
 R_OPTION="-R0/48/-2/2"      
!
! --- Header lines in file "stokes.dat" 
 H_OPTION="-H15" 
!
! --- Plot title  
 TITRE="'Rate of change of fully-normalized Stokes coefficients'"
!
! --- B-option 
 B_OPTION="-Ba4f2:'j(l,m)=l(l+1)/2+m+1':/f0.1a0.5:'d/dt (c@-lm@-, s@-lm@-) x 10@+-11@+/yr'"//":WSne"
! 
!
!
 Write(9,*) "" ; Write(9,*) "#--- Base of plot"
 Write(9,*) "psbasemap -X6 -Y4 -U'SELEN 3.2' "//trim(adjustl(R_OPTION))//" "//B_OPTION//"-JX18/14 -K > stokes.ps"
!
! A title (revised) 
 open(4,file='title_stokes.tmp',status='unknown') 
   write(4,*) "24 2.4 22 0 1 CB Rate of change of the fully-normalized Stokes coefficients"
 close(4) 
 Write(9,*) "pstext title_stokes.tmp -N -JX -R -G0 -O -K >> stokes.ps"
!
 Write(9,*) "" ; Write(9,*) "#--- Extracts cosine and sine coefficients"
 Write(9,*) "awk '{print $1, $4}' stokes.dat > cosine.tmp"
 Write(9,*) "awk '{print $1, $5}' stokes.dat >  sine.tmp"
!
 Write(9,*) "" ; Write(9,*) "#--- Plots labels for zonal degrees..."
!
!
  open(4,file='stokes0.tmp',status='unknown') 
  do 2 l=2, 48
  do 2 m=0, l
 	if(m==0) then 
	Write(4,'(i4,1x,a16,i2,a1,i1)') j_index(l,m), " -1.8 12 0 1 CM ", l, "-", m 
	Write(4,'(i4,1x,a16,i2,a1,i1)') j_index(l,m), " +1.8 12 0 1 CM ", l, "-", m 
	Endif
2 continue 
  close(4) 
 Write(9,*) "pstext stokes0.tmp -JX -R -G0 -O -W220 -K >> stokes.ps" 
!
 Write(9,*) "" ; Write(9,*) "#--- Plots open squares for sine components"
 Write(9,*) "psxy  sine.tmp ", trim(adjustl(H_OPTION)), " -B -R -JX  -W4 -K -O >> stokes.ps"
 Write(9,*) "psxy  sine.tmp ", trim(adjustl(H_OPTION)), " -B -R -JX -Ss0.35 -G0 -K -O >> stokes.ps"
 Write(9,*) "psxy  sine.tmp ", trim(adjustl(H_OPTION)), " -B -R -JX -Ss0.25 -G255 -K -O >> stokes.ps"
!
 Write(9,*) "" ; Write(9,*) "#--- Plots filled squares for cosine components"
 Write(9,*) "psxy cosine.tmp ", trim(adjustl(H_OPTION)), " -B -R -JX -Ss0.35 -G0 -K -O >> stokes.ps"
 Write(9,*) "psxy cosine.tmp ", trim(adjustl(H_OPTION)), " -B -R -JX  -W3 -K -O >> stokes.ps" 
!

 If(code=='-1') then 
 open(4,file='stokes1.tmp',status='unknown') 
 Write(4,*) "24 2.2 14 0 2 BC -Ice model: ", trim(adjustl(TITLICE)), & 
 " -LMAX=", trim(adjustl(DEGREE)),  " -RES=",  trim(adjustl(RESOLUTION)),& 	  
 " -ALMA rheology:", trim(adjustl(SHORT)), &  
 " -MODE=", trim(adjustl(MODE)),    " -ITER=", trim(adjustl(ITER)) 
 close(4) 
 else
 open(4,file='stokes1.tmp',status='unknown') 
 Write(4,*) "24 2.2 14 0 2 BC -Ice model: ", trim(adjustl(TITLICE)), & 
 " -Viscosity profile: ", trim(adjustl(VSTRING)),&
 " -LMAX=", trim(adjustl(DEGREE)),  " -RES=",  trim(adjustl(RESOLUTION)),& 	  
 " -NV=",   trim(adjustl(NV)),      " -CODE=", trim(adjustl(CODE)),& 
 " -MODE=", trim(adjustl(MODE)),    " -ITER=", trim(adjustl(ITER))    
 close(4) 
 endif
!
!
 Write(9,*) "" ; Write(9,*) "#--- A subtitle with some parameters"
 Write(9,*) "pstext stokes1.tmp -N -JX -R -G0 -O -K >> stokes.ps" 
!
 open(4,file='stokes2.tmp',status='unknown') ; write(4,*) "32 -1.00" ; close(4) 
 open(4,file='stokes3.tmp',status='unknown') ; write(4,*) "32 -1.25" ; close(4)  
 open(4,file='stokes4.tmp',status='unknown') 
 write(4,*) "33 -1.00  12 0 1 ML cosine" 
 write(4,*) "33 -1.25  12 0 1 ML  sine" 
 close(4) 
! 
 Write(9,*) "" ; Write(9,*) "#--- Plots a legend"
 Write(9,*) "psxy stokes2.tmp -N -JX -R -G0 -Ss0.35 -K -O >> stokes.ps" 
 Write(9,*) "psxy stokes3.tmp -N -JX -R -G0 -Ss0.35 -K -O >> stokes.ps" 
 Write(9,*) "psxy stokes3.tmp -N -JX -R -G255 -Ss0.25 -K -O >> stokes.ps" 
 Write(9,*) "pstext stokes4.tmp -N -JX -R -G0 -O -K >> stokes.ps"
!
 CLOSE(9) 
!
 END SUBROUTINE MAKE_STOKES 
!
!
!
!
!
!
 	SUBROUTINE SCAN_3D_REGIONS (FILE_IN, N_3D_REGIONS,       & 
					     NAME_3D_REGIONS,    &     
				 	     NP_3D_REGIONS,      & 
					     PIX_3D_FILENAMES)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
! # The routine scans the 'file_in' file to get information about the regions
!   across which the 3D velocity will be plotted. GS October 2008 for v. 2.7-
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! > Modified on August 2009 for the implementation in SELEN 3
! 
!
	IMPLICIT NONE 	
	INTEGER, PARAMETER :: LARGE_INTEGER = 300, & 
	 		      NMAX_3D_REGIONS=10,  & 
			      RES_MAX=48,          & 
			      NP_MAX=2*RES_MAX*(RES_MAX-1)*20+12
	INTEGER I, J, K, L, M, N, NOUT, IUNIT, N_3D_REGIONS, NP_3D_REGIONS(NMAX_3D_REGIONS)
	CHARACTER*30 FILE_IN, & 
		     NAME_3D_REGIONS (NMAX_3D_REGIONS), & 
		     PIX_3D_FILENAMES(NMAX_3D_REGIONS), & 
		     GRD_SIZE_KM_C(NMAX_3D_REGIONS)
	CHARACTER*20 	LON_MIN_C(NMAX_3D_REGIONS), &  
			LON_MAX_C(NMAX_3D_REGIONS), & 
			LAT_MIN_C(NMAX_3D_REGIONS), & 
			LAT_MAX_C(NMAX_3D_REGIONS) 
	INTEGER RES_3D_MAP(NMAX_3D_REGIONS) 
	REAL*4 LONXX(NP_MAX), LATXX(NP_MAX), GRD_SIZE_KM(NMAX_3D_REGIONS) 
	REAL*4 LON_3D(NMAX_3D_REGIONS,NP_MAX), LAT_3D(NMAX_3D_REGIONS,NP_MAX)
	REAL*4		LON_MIN(NMAX_3D_REGIONS), &  
			LON_MAX(NMAX_3D_REGIONS), & 
			LAT_MIN(NMAX_3D_REGIONS), & 
			LAT_MAX(NMAX_3D_REGIONS) 
        CHARACTER*200 LINE	
        CHARACTER*20 date, timc  	
	CHARACTER*4, PARAMETER :: SUFFIX = ".pix"
        INTEGER, PARAMETER :: NIMP=100
	CHARACTER*100 SS(NIMP)
	CHARACTER*1 YES_OR_NOT 

!				
!
	open(101,file=file_in,status='unknown') 
!
 	N_3D_REGIONS=0  
	i=0       
!
	do 10 j=1, large_integer 
!
        read(101,'(a200)',end=20) line 
!		
		call scan_string_soft (line, 7, ss, nout)
!	
		yes_or_not=trim(adjustl(ss(1))) 
!		
	        if(yes_or_not=='y') THEN 
!
		yes_or_not='n'

			N_3D_REGIONS = N_3D_REGIONS + 1 
!
			i=i+1 
!
!		# Write(*,'(i4,i4,a100)') i, N_3D_REGIONS, line 	

!
!		# Reads the region name & grid size 		
			name_3d_regions(i) = ss(2) 
			grd_size_km_c(i)   = ss(3) 
!
!		# Grid sixze in km (floating point) 
			call CHAR10_2_REAL(grd_size_km_c(i),grd_size_km(i))				
! 
!		# Lon-lat bounds for each region  			
			lon_min_c(i) = ss(4) 
			lon_max_c(i) = ss(5) 
			lat_min_c(i) = ss(6) 
			lat_max_c(i) = ss(7) 
!						 
!		# Converts in degrees (floating point)  			
			call CHAR10_2_REAL(lon_min_c(i),lon_min(i))				
			call CHAR10_2_REAL(lon_max_c(i),lon_max(i))				
			call CHAR10_2_REAL(lat_min_c(i),lat_min(i))				
			call CHAR10_2_REAL(lat_max_c(i),lat_max(i))				
		
!
!		# Tegmark resolution 			
			RES_3D_MAP(i)=INT(0.4*(6.371E3/grd_size_km(i)))				
			If(RES_3D_MAP(i)>RES_MAX) then 
				Write(*,*) "WARNING/ Maximum resolution exceeded in Sbr SCAN_3D_REGIONS"
				Write(*,*) "WARNING/ Resolution is set to the maximum allowed, i. e., ", RES_MAX
				RES_3D_MAP(i)=RES_MAX 				
			endif
!
!		# A name for the pixels file 			
			pix_3d_filenames(i)=trim(adjustl(name_3d_regions(i)))//SUFFIX
!
!		# Pixelization of the sphere 			
			call FINDPX(RES_3D_MAP(i), NP_3D_REGIONS(i), DBLE(LONXX), DBLE(LATXX))	
			LON_3D(i,:)=LONXX(:)
			LAT_3D(i,:)=LATXX(:)
!
!               # Logging the basic data
		        write(88,'(a30,1x,i3,1x,a14,i3,i6)') &
		        ' Region / Resolution / Pixels: ', i, & 
							  name_3d_regions(i), &
						          RES_3D_MAP(i), & 
							  NP_3D_REGIONS(i)
!
 	ENDIF
!
10      continue 	
20      close(101)  
!	
	write(88,*) 'Number of regions: ', N_3D_REGIONS
!
!
        call DATE_AND_TIME (date,timc)      
!
	do 30 i=1, N_3D_REGIONS 
!	
		iunit = 100+i 
		open(iunit, file=pix_3d_filenames(i), status='unknown') 
		Write(iunit,'(a1)')'#'
		Write(iunit,'(a1,1x,a27,1x,a20)')'#', 'Pixels lon/lat for region: ', &
				trim(adjustl(NAME_3D_REGIONS(I)))
		Write(iunit,'(a1,1x,a21,1x,a20)')'#', 'Grid resolution (km):', &
				trim(adjustl(grd_size_km_c(i)))	
                Write(iunit,'(a31,1x,a4,a1,a2,a1,a2,a1,a2,a1,a2,a1,a2,a1)') & 
		'# File created by config.f90 on ', & 
                date(1:4), '.', date(5:6), '.', date(7:8), ' ', & 
	        timc(1:2), '.', timc(3:4), '.', timc(5:6)			
		Write(iunit,'(a1)')'#'		
!
		do 40 j=1, NP_3D_REGIONS(i) 
!		
		if ((lon_3d(i,j).ge.lon_min(i).and.lon_3d(i,j).le.lon_max(i)).and. & 
		    (lat_3d(i,j).ge.lat_min(i).and.lat_3d(i,j).le.lat_max(i)))     &		 	   
		   Write(iunit, *) lon_3d(i,j), lat_3d(i,j)    
!
40	continue  
!
	close(iunit) 
!	
30      continue 
!		
 	END SUBROUTINE SCAN_3D_REGIONS 
!
!
!
!
!
!
 SUBROUTINE ALMA_2_CONFIG (wdir,            & 
                           visco_model,     & 
 		           n_mantle_layers, & 
			   lmax,            & 
			   option_d1)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
! **** This routine writes a configuration for "OPTION 2" of ALMA 
!      GS & FC July 21 2009 for v. 3.1 -
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
!      Re-touched for the implementation of degree one...  
!      April 13, 2010 by FC & GS
!       
 IMPLICIT NONE 
 CHARACTER*100  visco_model
 INTEGER n_mantle_layers
 CHARACTER*3 lmax 
 CHARACTER*100 WDIR 
 CHARACTER*1 OPTION_D1
!
! ======== Relevant parameters ========
 CHARACTER*1, PARAMETER :: TH='h'
 CHARACTER*1, PARAMETER :: ILOAD='1'
 CHARACTER*1, PARAMETER :: LMIN='1'
 CHARACTER*1, PARAMETER :: LSTP='1'
 CHARACTER*1, PARAMETER :: IH='1'
 CHARACTER*1, PARAMETER :: IL='1'
 CHARACTER*1, PARAMETER :: IK='1'
 CHARACTER*2, PARAMETER :: NSD='64'
 CHARACTER*2, PARAMETER :: NG='8'
 CHARACTER*1, PARAMETER :: ISALZ='1'
 CHARACTER*3, PARAMETER :: SCALE='lin'
! ======== Relevant parameters ========
!
write(*,*) trim(adjustl(wdir))//"/ALMA/alma.inc"

!OPEN(120,file=trim(adjustl(wdir))//"/ALMA/alma.inc",status='new')
!

OPEN(120,file="alma.inc",status='new')


Write(120,*)"! - - - - - - - - - - - - - - - - - - - - - - - -"	   
Write(120,*)"! This is file 'alma.inc' - SELEN port July 2009 "
Write(120,*)"! - - - - - - - - - - - - - - - - - - - - - - - -" 
Write(120,*)"USE FMZM"
Write(120,*)"integer, parameter :: imode=2   !  Mode 2: Multi-layer user-supplied model"    
Write(120,*)"! ~~~~~~~~~~~~~~~~~~~~~~~~~~"
Write(120,*)"! # Parameters for imode=2  "
Write(120,*)"! ~~~~~~~~~~~~~~~~~~~~~~~~~~"
Write(120,*)"character*100, parameter :: mod2 = ", trim(adjustl(visco_model)) 
Write(120,*)"! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
Write(120,*)"! # Parameters for imode=1 & imode=3 (declared but not used) "
Write(120,*)"! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
Write(120,*)"integer, parameter ::      kv=0       ! DUMMY PARAMETER Type of viscosity profile"
Write(120,*)"integer, parameter ::      ioc=0      ! DUMMY PARAMETER PREM switch (1 oceanic, 0 continental)"
Write(120,*)"real, parameter::          lth=0      ! DUMMY PARAMETER Lithospheric thickness (m)"
Write(120,*)"character*30, parameter :: mod3='0'   ! DUMMY PARAMETER nla=2, user supplied model"
Write(120,*)"integer, parameter ::      ire=0      ! DUMMY PARAMETER Mantle rheology (Table 7)" 
Write(120,*)"!~~~~~~~~~~~~~~~~~~~~~"
Write(120,*)"!# General parameters "
Write(120,*)"!~~~~~~~~~~~~~~~~~~~~~"
Write(120,*)"integer, parameter :: nla=    ", n_mantle_layers,  "     ! Number of mantle layers"
Write(120,*)"character*1, parameter :: th= ", "'"//th//"'", 	        "     ! Stands for 'Heaviside' Love numbers "
Write(120,*)"integer, parameter :: iload=  ", ILOAD,            "     ! Loading (1) or tidal (0) Love numbers "          
Write(120,*)"integer, parameter ::    L1=  ", LMIN,               "     ! Minimum degree"
Write(120,*)"integer, parameter ::    L2=  ", LMAX,              "     ! Maximum degree"
Write(120,*)"integer, parameter ::    LS=  ", LSTP,              "     ! Step degree"
Write(120,*)"integer, parameter ::    IH=  ", IH,                 "     ! h Love number (0/1)"
Write(120,*)"integer, parameter ::    IL=  ", IL,                 "     ! l Love number (0/1)"
Write(120,*)"integer, parameter ::    IK=  ", IK,                 "     ! k Love number (0/1)"
Write(120,*)"integer, parameter ::   NSD=  ", NSD,                "     ! Significant digits"
Write(120,*)"integer, parameter ::    NG=  ", NG,                "     ! Order of the Gaver method"
Write(120,*)"integer, parameter :: ISALZ=  ", ISALZ,            "      !Salzer acceleration"
Write(120,*)"character*3, parameter :: tscale= ", "'"//scale//"'",       "   ! toggle between lin/log time scale"
!
! Note: Unit 120 will be closed within the body of program "config.f90" - GS & FC July 2009 
!
END SUBROUTINE ALMA_2_CONFIG  	
!
!
!  
!
!
!
 SUBROUTINE MAKE_ELAREB_MAPS_ANTARCTICA (option_reb_ag, & 
	  			         option_reb_ar, & 
 	  			         resolution,    & 
	   				 nv,            & 
	  				 code,          & 
	  				 iter,          & 
	  				 mode,          &
	  				 degree,        & 
	  				 vstring,       & 
					 TITLICE,       & 
					 FILE_GMT) 
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! Prepares a GMT script for *** GLOBAL and REGIONAL maps *** showing variations of 
! of S, U & N in response to melting of the Antarctic ice sheet *** GS May 2011 ** 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Copied and adapted from "MAKE_ELAREB_MAPS_GREENLAND" and philosophycally 
! based on "SUBROUTINE MAKE_GMAPS"  ******* Last change GS May 03 2011  ******
!
!
 IMPLICIT NONE
 CHARACTER*20  FILE_GMT
 CHARACTER*1   OPTION_REB_AG, OPTION_REB_AR
 CHARACTER*30  NAMEIN, NAMEOUT, NAMEF   
 CHARACTER*80  TITRE
 INTEGER I 
!
 CHARACTER*100 SHORT_VISCO 
 CHARACTER*10  TITLICE, RESOLUTION, DEGREE
 CHARACTER*30  VSTRING
 CHARACTER*20  R_OPTION, T_OPTION, J_OPTION
 CHARACTER*3   NV, CODE
 CHARACTER*1   OPTION_ROF
 CHARACTER*1   ITER, MODE 
!
 INTEGER K
!
!
!
!
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
! A unique file for global and regional analyses 
!
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
 OPEN (9,FILE=FILE_GMT,STATUS='UNKNOWN')
! 
 Write(9,*)"gmtset PAPER_MEDIA A4+"
 Write(9,*)"gmtset LABEL_FONT_SIZE 24p"
 Write(9,*)"gmtset ANOT_FONT_SIZE 16p"
 Write(9,*)"gmtset FRAME_WIDTH 0.1c"
 Write(9,*)" "
!
!
!




!
! ///////////////////////////////////////// ---------- Global Antarctica!  Global!
!        Global map for Antarctica 
! ///////////////////////////////////////// ---------- Global Antarctica!  Global!
! 
 IF(OPTION_REB_AG=='y') THEN 
! 
  T_OPTION="-T-0.4/0.4/0.05"         ! Range of the palette 
  R_OPTION="-R0/360/-80/80"     ! Range of plot  
! 
 Write(9,*) "# "
 Write(9,*) "# "
 Write(9,*) "# "
 Write(9,*) "# GLOBAL MAP - GLOBAL MAP - GLOBAL MAP - GLOBAL MAP - GLOBAL MAP "
 Write(9,*) "# GLOBAL MAP - GLOBAL MAP - GLOBAL MAP - GLOBAL MAP - GLOBAL MAP "
 Write(9,*) "# GLOBAL MAP - GLOBAL MAP - GLOBAL MAP - GLOBAL MAP - GLOBAL MAP "
 Write(9,*) "# "
!
 DO 10 I=1, 3 
!
 if(i==1)namein ="smap_anta_glob.dat" 
 if(i==2)namein ="umap_anta_glob.dat" 
 if(i==3)namein ="nmap_anta_glob.dat" 
 if(i==1)nameout="smap_anta_glob.ps" 
 if(i==2)nameout="umap_anta_glob.ps" 
 if(i==3)nameout="nmap_anta_glob.ps" 
 if(i==1)namef="gtmpfs.dat"
 if(i==2)namef="gtmpfu.dat" 
 if(i==3)namef="gtmpfn.dat"  
 if(i==1)TITRE="Rate of relative SL change (\dot S)"
 if(i==2)TITRE="Rate of vertical uplift (\dot U)"
 if(i==3)TITRE="Rate of absolute SL change (\dot N)"
!
 Write(9,*) " "
 if(i==1)Write(9,*)"# ---- Global map of S at present time for ELASTIC REBOUND in ANTARCTICA ----" 		      	   
 if(i==2)Write(9,*)"# ---- Global map of U at present time for ELASTIC REBOUND in ANTARCTICA ----" 		      	   
 if(i==3)Write(9,*)"# ---- Global map of N at present time for ELASTIC REBOUND in ANTARCTICA ----" 		      	   
!
 if(i==1)Write(9,*) "echo", "    - S for Antarctica - global" 
 if(i==2)Write(9,*) "echo", "    - U for Antarctica - global" 
 if(i==3)Write(9,*) "echo", "    - N for Antarctica - global" 
!
 Write(9,*) " "
 Write(9,*) "makecpt -Cno_green ", trim(adjustl(T_OPTION)), " > pale.cpt" 
 Write(9,*) "psbasemap -X3 -Y5 -Ba180/a40WSEn -Jm0.018i ", & 
             trim(adjustl(R_OPTION)), " -K > ",            & 
	     trim(adjustl(nameout))  
!
 Write(9,*) "pscontour -I -Jm -O -K ", & 
             trim(adjustl(R_OPTION)), " ", & 
             trim(adjustl(namein)), " -Cpale.cpt  >> ", &   
	     trim(adjustl(nameout)) 
!
 Write(9,*) "pscontour -O -K -G8 -W2/0/100/0 -Jm -A+g255 ",         & 
             trim(adjustl(R_OPTION)), " ", & 
             trim(adjustl(namein)), " -Cpale.cpt  >> ", & 
	     trim(adjustl(nameout)) 
!
 Write(9,*) "pscoast -Jm -Dc -B -W2/0  -A1000 -O -K ",  &
             trim(adjustl(R_OPTION)), " >> ",           & 
             trim(adjustl(nameout))  
!
! A file for "titles" 
 open (4,file=NAMEF,status='unknown')  
!
! Alma 
 If(CODE=='-1')  THEN 
		Write(4,*) "180 +84.2 26 0 3 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	 	  "  -Elasticity: ALMA ", trim(adjustl(SHORT_VISCO))  
		Write(4,*) "180 -87.5   13  0 0 BC ",          & 
 	   	 	   " -LMAX=", trim(adjustl(DEGREE)),   & 
 	    		   " -RES=", trim(adjustl(RESOLUTION)),& 
	    	 	   " -NV=", trim(adjustl(NV)),         & 
	    		   " -MODE=", trim(adjustl(MODE)),     &   
	    	 	   " -ITER=", trim(adjustl(ITER))
		close(4) 
!
! Taboo 
 ELSEIF(CODE/='-1'.and.CODE/='-2')THEN 
		Write(4,*) "180 +84.2 26 0 3 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	   	"  -Elasticity: TABOO "  
		Write(4,*) "180 -87.5   13  0 0 BC ",          & 
 	    		   " -LMAX=", trim(adjustl(DEGREE)),   & 
 	    		   " -RES=", trim(adjustl(RESOLUTION)),& 
	    		   " -NV=", trim(adjustl(NV)),         & 
	    		   " -CODE=", trim(adjustl(CODE)),     & 
	    		   " -MODE=", trim(adjustl(MODE)),     &   
	    		   " -ITER=", trim(adjustl(ITER))
		close(4) 
! External 
 ELSEIF(CODE=='-2')THEN 
		Write(4,*) "180 +84.2 26 0 3 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	   	"  -Elasticity: External "  
		Write(4,*) "180 -87.5   13  0 0 BC ",          & 
 	    		   " -LMAX=", trim(adjustl(DEGREE)),   & 
 	    		   " -RES=", trim(adjustl(RESOLUTION)),& 
	    		   " -CODE=", trim(adjustl(CODE)),     & 
	    		   " -MODE=", trim(adjustl(MODE)),     &   
	    		   " -ITER=", trim(adjustl(ITER))
		close(4) 
 ENDIF
		Write(9,*) "pstext -N -R -Jm -B ", trim(adjustl(NAMEF)),    & 
		           " -G0 -O -K >> ", trim(adjustl(nameout)) 
		Write(9,*) "psscale -U/0.5/0.5/'SELEN 3.2' -E -Cpale.cpt ", & 
		           "-B1f0.1a0.2/:mm/yr: -D8.25/-1/10/1h -O -K >> ",   & 
	                    trim(adjustl(nameout)) 
 		Write(9,*) "psbasemap -J -R -B::/::EWSN -O >> ", trim(adjustl(nameout)) 
 10 CONTINUE
!
!
!
! \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/
!   A script for the Normalized Sea Level change 
! \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/
! \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/
!   A script for the Normalized Sea Level change 
! \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/
! \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/
!   A script for the Normalized Sea Level change 
! \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/ \/\/\/

!
!
 Write(9,*) "# "
 Write(9,*) "# "
 Write(9,*) "# "
 Write(9,*) "# GLOBAL NORMALIZED SEA LEVEL CHANGE MAP - GLOBAL NORMALIZED SEA LEVEL CHANGE MAP "
 Write(9,*) "# GLOBAL NORMALIZED SEA LEVEL CHANGE MAP - GLOBAL NORMALIZED SEA LEVEL CHANGE MAP "
 Write(9,*) "# GLOBAL NORMALIZED SEA LEVEL CHANGE MAP - GLOBAL NORMALIZED SEA LEVEL CHANGE MAP "
 Write(9,*) "# "
!
 namein  = "smap_anta_glob_norm.dat" 
 nameout = "smap_anta_glob_norm.ps" 
 namef   = "ngtmpfs.dat"  
 TITRE   = "Sea Level Change / Eustatic Change"
!
 Write(9,*)" "
 Write(9,*)"gmtset PAPER_MEDIA A4+"
 Write(9,*)"gmtset LABEL_FONT_SIZE 24p"
 Write(9,*)"gmtset ANOT_FONT_SIZE 16p"
 Write(9,*)"gmtset FRAME_WIDTH 0.1c"
 Write(9,*)" "
!
 T_OPTION="-T-0.1/1.3/0.1"         ! Range of the palette 
 R_OPTION="-R0/360/-80/85"         ! Range of plot  
!
 Write(9,*) "echo", "    - S normalized for Antarctica - global" 
!
 Write(9,*) " "
 Write(9,*) "makecpt -Cno_green ", trim(adjustl(T_OPTION)), " > pale_norm.cpt" 
 Write(9,*) "psbasemap -X3 -Y5 -JQ0/18 -Ba90/a30WSEn  ", & 
             trim(adjustl(R_OPTION)), " -K > ",            & 
	     trim(adjustl(nameout))  
 Write(9,*) "pscontour -I -O -K  -JQ ",         & 
             trim(adjustl(R_OPTION)), " ",      & 
             trim(adjustl(namein)), " -Cpale_norm.cpt  >> ", &   
	     trim(adjustl(nameout)) 
 Write(9,*) "pscoast -Dc -B -W2/0 -JQ -A10000 -G200 -O -K ", & 
             trim(adjustl(R_OPTION)), " >> ",      & 
             trim(adjustl(nameout))  
!
!
! A file for "titles" 
 open (4,file=NAMEF,status='unknown')  
!
! Alma 
 If(CODE=='-1')  THEN 
		Write(4,*) "0 +100     22 0 0 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "0 -175     14 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	 	  "  -Elasticity: ALMA ", trim(adjustl(SHORT_VISCO))  
		Write(4,*) "0 -160.5   13  0 0 BC  ",          & 
 	   	 	   " -LMAX=", trim(adjustl(DEGREE)),   & 
 	    		   " -RES=", trim(adjustl(RESOLUTION)),& 
	    	 	   " -NV=", trim(adjustl(NV)),         & 
	    		   " -MODE=", trim(adjustl(MODE)),     &   
	    	 	   " -ITER=", trim(adjustl(ITER))
		close(4) 
!
! Taboo 
 ELSEIF(CODE/='-1'.and.CODE/='-2')THEN 
		Write(4,*) "0 +100     22 0 0 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "0 -175     14 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	   	"  -Elasticity: TABOO "  
		Write(4,*) "0 -160.5   13  0 0 BC ",          & 
 	    		   " -LMAX=", trim(adjustl(DEGREE)),   & 
 	    		   " -RES=", trim(adjustl(RESOLUTION)),& 
	    		   " -NV=", trim(adjustl(NV)),         & 
	    		   " -CODE=", trim(adjustl(CODE)),     & 
	    		   " -MODE=", trim(adjustl(MODE)),     &   
	    		   " -ITER=", trim(adjustl(ITER))
		close(4) 
! External 
 ELSEIF(CODE=='-2')THEN 
		Write(4,*) "0 +100     22 0 0 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "0 -175     14 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	   	"  -Elasticity: External "  
		Write(4,*) "0 -160.5   13  0 0 BC ",          & 
 	    		   " -LMAX=", trim(adjustl(DEGREE)),   & 
 	    		   " -RES=", trim(adjustl(RESOLUTION)),& 
	    		   " -CODE=", trim(adjustl(CODE)),     & 
	    		   " -MODE=", trim(adjustl(MODE)),     &   
	    		   " -ITER=", trim(adjustl(ITER))
		close(4) 
 ENDIF
!
!
!
 Write(9,*) "pstext -N -R -JQ -B ", trim(adjustl(NAMEF)),    & 
            " -G0 -O -K >> ", trim(adjustl(nameout)) 
 Write(9,*) "psscale -U/0.5/0.5/'SELEN 3.2' -E -Cpale_norm.cpt ", & 
		           "-B1f0.1a0.1/:S/S@-E@-: -D9/-1.5/-16/1h -O -K >> ",   &    
	                    trim(adjustl(nameout)) 
 Write(9,*) "psbasemap -JQ -R -B::/::EWSN -O >> ", trim(adjustl(nameout)) 


! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ---------- Global Antarctica!
 ENDIF    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ---------- Global Antarctica!    
!
!
!
 
 
!
! ///////////////////////////////////////// ---------- Regional Antarctica!
!        Regional map for Antarctica 
! ///////////////////////////////////////// ---------- Regional Antarctica!
! 
 IF(OPTION_REB_AR=='y') THEN 
!
 Write(9,*) "# "
 Write(9,*) "# "
 Write(9,*) "# "
 Write(9,*) "# REGIONAL MAP - REGIONAL MAP - REGIONAL MAP - REGIONAL MAP - REGIONAL MAP# "
 Write(9,*) "# REGIONAL MAP - REGIONAL MAP - REGIONAL MAP - REGIONAL MAP - REGIONAL MAP# "
 Write(9,*) "# REGIONAL MAP - REGIONAL MAP - REGIONAL MAP - REGIONAL MAP - REGIONAL MAP# "
 Write(9,*) "# "
!
 Write(9,*)"gmtset PAPER_MEDIA A4+"
 Write(9,*)"gmtset LABEL_FONT_SIZE 24p"
 Write(9,*)"gmtset ANOT_FONT_SIZE 16p"
 Write(9,*)"gmtset FRAME_WIDTH 0.1c"
 Write(9,*)" "
!
  T_OPTION="-T-1/1/0.1"         ! Range of the palette 
  R_OPTION="-R0/360/-90/-52"    ! Range of plot   
  J_OPTION="-JE0/-90/12"
!
 DO 20 I=1, 3 
!
 if(i==1)namein ="smap_anta_reg.dat" 
 if(i==2)namein ="umap_anta_reg.dat" 
 if(i==3)namein ="nmap_anta_reg.dat" 
 if(i==1)nameout="smap_anta_reg.ps" 
 if(i==2)nameout="umap_anta_reg.ps" 
 if(i==3)nameout="nmap_anta_reg.ps" 
 if(i==1)namef="rrtmpfs.dat"
 if(i==2)namef="rrtmpfu.dat" 
 if(i==3)namef="rrtmpfn.dat"  
 if(i==1)TITRE="Rate of sea level change"
 if(i==2)TITRE="Vertical velocity"
 if(i==3)TITRE="Rate of geoid change"
!
 Write(9,*) " "
 if(i==1)Write(9,*)"# ---- Regional map of S at present time for ELASTIC REBOUND in Antarctica ----" 		      	   
 if(i==2)Write(9,*)"# ---- Regional map of U at present time for ELASTIC REBOUND in Antarctica ----" 		      	   
 if(i==3)Write(9,*)"# ---- Regional map of N at present time for ELASTIC REBOUND in Antarctica ----" 		      	   
!
!
!
 if(i==1)Write(9,*) "echo", "    - S for Antarctica - regional" 
 if(i==2)Write(9,*) "echo", "    - U for Antarctica - regional" 
 if(i==3)Write(9,*) "echo", "    - N for Antarctica - regional" 
!
!
 Write(9,*) " "
 Write(9,*) "makecpt -Cno_green ", trim(adjustl(T_OPTION)), " > pale.cpt" 
 Write(9,*) "psbasemap -P -X1 -Y6 -Ba90g45/a90g45WSEN ",     trim(adjustl(R_OPTION)), " ", & 
             trim(adjustl(J_OPTION)), " -K > ", " ", trim(adjustl(nameout))  
 Write(9,*) "pscontour -I -O -K ",         & 
             trim(adjustl(R_OPTION)), " ", & 
   	     trim(adjustl(J_OPTION)), " ", &
             trim(adjustl(namein)), " -Cpale.cpt  >> ", & 
	     trim(adjustl(nameout)) 
!
 Write(9,*) "pscontour -O -K -G8 -W2/0/100/0 -A+g255 ",         & 
             trim(adjustl(R_OPTION)), " ", & 
   	     trim(adjustl(J_OPTION)), " ", &
             trim(adjustl(namein)), " -Cpale.cpt  >> ", & 
	     trim(adjustl(nameout)) 
!	     
 Write(9,*) "pscoast -Dh -B -W2/0 -A1000 -O -K ", & 
             trim(adjustl(R_OPTION)), " ",         & 
             trim(adjustl(J_OPTION)), " ", " >> ", & 
	     trim(adjustl(nameout))  
!
! A file for "titles" 
 open (4,file=NAMEF,status='unknown')  
!
! Alma 
 If(CODE=='-1')  THEN 
		Write(4,*) "0   -42      24  0 2 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "180 -42      16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	 	  "  -Elasticity: ALMA ", trim(adjustl(SHORT_VISCO))  
		Write(4,*) "319 +52.5   13  0 0 BC ",& 
 	   	 	   " -LMAX=", trim(adjustl(DEGREE)),& 
 	    		   " -RES=", trim(adjustl(RESOLUTION)),& 
	    		   " -NV=", trim(adjustl(NV)),& 
	    		   " -MODE=", trim(adjustl(MODE)),&   
	    		   " -ITER=", trim(adjustl(ITER))
		close(4) 
!
! Taboo 
 ELSEIF(CODE/='-1'.and.CODE/='-2')THEN 
 		Write(4,*) "0	-42  	 24  0 2 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "180 -42  	 16 0 0 BC   -Ice model: ", trim(adjustl(TITLICE)),& 
 	   	"  -Elasticity: TABOO "  
		Write(4,*) "180 -38   13  0 0 BC ",& 
 	    		   " -LMAX=", trim(adjustl(DEGREE)),& 
 	    		   " -RES=", trim(adjustl(RESOLUTION)),& 
	    	 	   " -NV=", trim(adjustl(NV)),& 
	    		   " -CODE=", trim(adjustl(CODE)),& 
	    		   " -MODE=", trim(adjustl(MODE)),&   
	    		   " -ITER=", trim(adjustl(ITER))
		close(4) 
!
 ELSEIF(CODE=='-2')THEN 
 		Write(4,*) "0	-42  	 24  0 2 BC ", trim(adjustl(TITRE)) 
		Write(4,*) "180 -42  	 16 0 0 BC   -Ice model: ", trim(adjustl(TITLICE)),& 
 	   	"  -Elasticity: External "  
		Write(4,*) "180 -38   13  0 0 BC ",& 
 	    		   " -LMAX=", trim(adjustl(DEGREE)),& 
 	    		   " -RES=", trim(adjustl(RESOLUTION)),& 
	    		   " -CODE=", trim(adjustl(CODE)),& 
	    		   " -MODE=", trim(adjustl(MODE)),&   
	    		   " -ITER=", trim(adjustl(ITER))
		close(4) 
 ENDIF
!
		Write(9,*) "pstext -N -R -B ", trim(adjustl(NAMEF)), " ", &
		                               trim(adjustl(J_OPTION)), " -G0 -O -K >> ", & 
					       trim(adjustl(nameout)) 
		Write(9,*) "psscale -U/0.5/0.5/'SELEN 3.2' -E -Cpale.cpt -B1f0.1a0.5/:mm/yr: -D14/6/6/1  -O -K >> ", & 
	                    trim(adjustl(nameout)) 
 		Write(9,*) "psbasemap ", trim(adjustl(J_OPTION)), " ", & 
		           " -R -B::/::EWSN -O >> ", trim(adjustl(nameout)) 
!
!
20 CONTINUE 
!
!

!
!
 IF(OPTION_REB_AG=='y'.or.OPTION_REB_AR=='y') CLOSE(9)
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ---------- On Regional Antarctica!
 ENDIF    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ---------- On Regional Antarctica!

!
!
!
!
 END SUBROUTINE MAKE_ELAREB_MAPS_ANTARCTICA
!
!
!
!
!
!
