!
! This is program "TGAUGES_ELA_MS.F90" 
! 
! Created as a cpy of "TGAUGES.F90" by GS on February 6, 2012 
! (for the implementaion of the multi-step rebound at tide-gauges)
! GS March 08 2012 - REAL -> DBLE and expressions rearranged to preserve precision in
!                    the computation of time derivatives. 
! GS March 2012- Added information on vertical displacement  
! Option GIA on April 19 (Urbino Hospital...) 
!
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
! Copyright (C) 2008 Giorgio Spada, Florence Colleoni, and Paolo Stocchi 
!
! This file is part of SELEN. 
!  
! SELEN is free software: you can redistribute it and/or modify it under the 
! terms of the GNU General Public License as published by the Free Software 
! Foundation, either version 3 of the License, or at your option) any later 
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
! -------------------------------------------------------------------------
!  Computes the present-day rate of sealevel change at the PSMSL sites ... 
! -------------------------------------------------------------------------
!
! Input files:
!	- rlr-trends.txt  
! 	- shs.bin
! 	- shn.bin
! 	- shu.bin
! 	- shtidegauges.bin
!
! Output files:
!       - psmsl*.dat
!
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! INCLUDE "harmonics.f90"
  PROGRAM PSMSL
  IMPLICIT NONE
  INCLUDE "data.inc"
  CHARACTER*1, PARAMETER :: NORTH='N', SOUTH='S', WEST='W', EAST='E'
  INTEGER, PARAMETER :: MAXN=10000  ! A large number   
  CHARACTER*100 RIGA		    ! A row 
  CHARACTER*8 PCODE                 ! PSMSL code  
  CHARACTER*3 GLOSS_CODE            ! GLOSS code 
  CHARACTER*3 YEARS                 ! Number of years 
  CHARACTER*11 RANGE_YEARS          ! Range of years 
  CHARACTER*9 DATUM 	            ! Trend 
  CHARACTER*3 PLUS_MINUS 	    ! +/- 
  CHARACTER*5 ERROR	            ! Error 
  CHARACTER*7 STDV 	            ! Std. deviation of residues 							
  CHARACTER*2 LAT1, LAT2            ! Latitide degrees and minutes 	
  CHARACTER*3 LON1	            ! Longitude degrees 
  CHARACTER*2 LON2 		    ! Longitude minutes 
  CHARACTER*30 NAME                 ! Name of the station 
  CHARACTER*1 WLAT, WLON            ! Where is the station (e, w, s, n)  
  CHARACTER*20 DATE, TIMC           ! DATE AND TIME
  REAL*8 LAT, NLAT, DLAT            ! Latitude  
  REAL*8 LON, NLON, DLON            ! Longitude 
  REAL*8 PSDOT, PNDOT, PUDOT        ! Predicted rate 
  INTEGER I, J, N, NH, DOM
  REAL*8, PARAMETER :: DDELTA=DELTA 
!
! Updated February 2012 with "nn+1"
!
  COMPLEX*16 Y(JMAX), SA(JMAX,0:NN+1), & 
  		      NA(JMAX,0:NN+1), & 
		      UA(JMAX,0:NN+1)  
!
! ---- Declarations for the file format '1' (new as of may 13, 2011) 
  INTEGER IJ, STAT_ID, TMIN, TMAX, NVAL
  REAL*8 STAT_LON, STAT_LAT, TREND, DTREND
  CHARACTER*30 STAT_NAME
  CHARACTER*15 INDIVIDUAL_FILENAME 
  CHARACTER*10 FAKE_NAME  
  CHARACTER*1 CA(6) 
  INTEGER K
!
! ----
  INTEGER KT
  REAL*8 SVAR, NVAR, UVAR, SDOT, NDOT, UDOT 
!
!
! PATH to the GIA corrections and misc variables [New as of April 21, 2012]
 CHARACTER*30, PARAMETER :: GIA_PATH="./DATA/GIA-corrections/" 
 CHARACTER*40, PARAMETER :: file_gia_corr = trim(adjustl(GIA_PATH))//"p3.dat"
 INTEGER, PARAMETER :: NH_FILE_GIA_CORR=17 
 REAL*8 SDOT_GIA_CORR, UDOT_GIA_CORR, NDOT_GIA_CORR  
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!
!
!
!
! --- SH file for PSMSL sites (input) 
        open(17,file='shtidegauges.bin',status='unknown',form='unformatted')  
!
! --- Writes an header for the output file 
!     This file "ptidegauges.dat" contains predictions for ALL TGs in database 	
	OPEN (99,FILE='ptidegauges.dat',STATUS='unknown')
! 
!
! --- Header information about the SELEN settings (cumulative file)
!
        call DATE_AND_TIME (date,timc)      
	Write (99,*)  "# **************************************************"
        Write (99,*)  "# ", date(1:4), '.', date(5:6), '.', date(7:8), & 
	             ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6)
        Write (99,*)  "# -----> ELASTIC REBOUND in MULTI-STEPS <-----"
	Write (99,*)  "# Ice model: ", trim(adjustl(ice_model))
	Write (99,*)  "# SLE iterations: ", SMAX 
	Write (99,*)  "# SLE mode of solution: ", IMODE
	Write (99,*)  "# Maximum harmonic degree: ", LMAX 
	Write (99,*)  "# Tegmark resolution: ", RES 
	Write (99,*)  "# Elastic model code: ", CDE
	IF(GIA_CORR==1) & 
	Write (99,*)  "# ====== a GIA correction IS applied ======"
	IF(GIA_CORR==0) & 
	Write (99,*)  "# ====== NO GIA corrections applied ======"
	Write (99,*)  "# **************************************************"
!
!
! --- tgauges.f: reading sealevel CSH coefficients from file', filename 
        open(3,file='shs.bin',status='unknown',form='unformatted') ; read(3) SA ; close(3)
        open(3,file='shn.bin',status='unknown',form='unformatted') ; read(3) NA ; close(3)
        open(3,file='shu.bin',status='unknown',form='unformatted') ; read(3) UA ; close(3)
!
! --- Counting the header lines beginning with '#' and the data lines
        open(10,file=TGAUGES_DATABASE,status='unknown')
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
        open(10,file=TGAUGES_DATABASE,status='unknown')
!
! --- Reading again the header lines  
        do i=1, nh 
	    read(10,'(a100)',end=1) riga 
        enddo  
!
!
! GIA Correction GIA correction GIA Correction GIA correction GIA Correction GIA correction 
! GIA Correction GIA correction GIA Correction GIA correction GIA Correction GIA correction 
! GIA Correction GIA correction GIA Correction GIA correction GIA Correction GIA correction 
!
!
! --- If the GIA correction is requested, reads the header of the corrections filename 
!     The number of header lines is determined by the file format (see TGAUGES.F90).  
!
        IF    (GIA_CORR==1) THEN  
!
		OPEN(111, file=file_gia_corr,status='unknown')  
!
		DO I=1, NH_FILE_GIA_CORR 
	            read(111,'(a100)',end=1) riga 
		ENDDO
!
        ENDIF
!
!
!
! === Reading all lines of TGAUGES_DATABASE GS = October 5 2007 - 
!
! *************************************************************************************************
! *************************************************************************************************
!
!
!
         IF(TG_DATABASE_FORMAT=='1')   THEN  !  ------- Input TG file has format '1'
!
!
!
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
         DO 6 I=1, N   ! Loop on the tide-gauges 
!
! //////////////////////////////////////////////////////
!
!
! GIA Correction GIA correction GIA Correction GIA correction GIA Correction GIA correction 
! GIA Correction GIA correction GIA Correction GIA correction GIA Correction GIA correction 
! GIA Correction GIA correction GIA Correction GIA correction GIA Correction GIA correction 
!
! --- GIA corrections from the external file 
!
         IF(GIA_CORR==1) THEN 
	 READ(111,202)  IJ, STAT_ID, STAT_LON, STAT_LAT, & 
		            TMIN,    TMAX,  NVAL,        & 
		            TREND,    DTREND,	         & 		   
                            SDOT_GIA_CORR, &
                            UDOT_GIA_CORR, &
                            NDOT_GIA_CORR, &  		  			   
		            STAT_NAME
         ELSE 
			    SDOT_GIA_CORR=0D0
			    UDOT_GIA_CORR=0D0
			    NDOT_GIA_CORR=0D0	 
         ENDIF
!
! GIA Correction GIA correction GIA Correction GIA correction GIA Correction GIA correction 
! GIA Correction GIA correction GIA Correction GIA correction GIA Correction GIA correction 
! GIA Correction GIA correction GIA Correction GIA correction GIA Correction GIA correction 
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
	 FAKE_NAME=CA(1)//CA(2)//CA(3)//CA(4)//CA(5)//CA(6)//".tga"



         open(3,file='junk.dat',status='unknown') 
             if(              i.le.9)     write(3,'(A3,i1,A1,A10)') "000",i,"-",FAKE_NAME 
             if(10  .le.i.and.i.le.99)    write(3,'(A2,i2,A1,A10)')  "00",i,"-",FAKE_NAME 
             if(100 .le.i.and.i.le.999)   write(3,'(A1,i3,A1,A10)')   "0",i,"-",FAKE_NAME  
             if(1000.le.i.and.i.le.9999)  write(3,'(   i4,A1,A10)')       i,"-",FAKE_NAME
         close(3)
	     open(3,file='junk.dat',status='unknown')
	     read(3,'(a15)') individual_filename
	 close(3)    
!	 
! --- write(*,*) i,  individual_filename     
!
         OPEN(200,file=individual_filename,status='UNKNOWN')
!
! --- Header for each individual TG site in the cumulative file
         Write(99,*) " "
         Write(99,'(A2,A10,1X, A30)')           '# ','Stat name:',             STAT_NAME
         Write(99,'(A2,A19,1X, 2(i5, 1x))')     '# ','Stat number and ID:',    IJ, STAT_ID
         Write(99,'(A2,A14,1X, 2(f10.4, 1x))')  '# ','Lon/lat (deg):',         STAT_LON, STAT_LAT
         Write(99,'(A2,A17,1X, 3(i5, 1x))')     '# ','tmin, tmax, nval:',      TMIN, TMAX,  NVAL	 
         Write(99,'(A2,A22,1X, 2(F7.2, 1x))')   '# ','trend & sigma (mm/yr):', TREND, DTREND 
	 Write(99,'(A2,A12)')                   '# ','Predictions:'
!
! --- The same header is also reported in each individual filename 
         Write(200,'(A2,A10,1X, A30)')           '# ','Stat name:',             STAT_NAME
         Write(200,'(A2,A19,1X, 2(i5, 1x))')     '# ','Stat number and ID:',    IJ, STAT_ID
         Write(200,'(A2,A14,1X, 2(f10.4, 1x))')  '# ','Lon/lat (deg):',         STAT_LON, STAT_LAT
         Write(200,'(A2,A17,1X, 3(i5, 1x))')     '# ','tmin, tmax, nval:',      TMIN, TMAX,  NVAL	 
         Write(200,'(A2,A22,1X, 2(F7.2, 1x))')   '# ','trend & sigma (mm/yr):', TREND, DTREND 
	 Write(200,'(A2,A12)')                   '# ','Predictions:'

!
! --- SH at the current PSMSL site...
      	        read(17) y
!
! --- Predicted rate at the site 'i' 
		If(n>=50)                then 
      	        If(mod(i,50)==0.or.i==1) write(*,*) '    - ', i, 'of', n, stat_name 
		elseif(10<=n.and.n<=50)  then
      	        If(mod(i,10)==0.or.i==1) write(*,*) '    - ', i, 'of', n, stat_name 		
		elseif(n<10)             then
		write(*,*) '    - ', i, 'of', n, stat_name 			       
					 endif 		
!
! Header for the predictions on the ***cumulative*** file
!
       Write (99,'(A2,A81)')   "# ", "k     time       S           N           U        dot S       dot N       dot U "
       Write (99,'(A2,A81)')   "# ", "      (yr)      (m)         (m)         (m)      (mm/yr)     (mm/yr)     (mm/yr)"
!
!
! The same header for the predictions on the ***individual*** file
!
!
       Write(200,'(A2,A81)')   "# ", "k     time       S           N           U        dot S       dot N       dot U "
       Write(200,'(A2,A81)')   "# ", "      (yr)      (m)         (m)         (m)      (mm/yr)     (mm/yr)     (mm/yr)"
!
!
!
       DO 100 KT=0, NN 
!		
		if(KT.GE.0.and.KT.le.NN) then 
		   SVAR=0d0 
		   NVAR=0d0
		   UVAR=0D0
		Endif 
!
		if(KT.GE.1.and.KT.le.NN-1)then 
		   SDOT=0d0 
		   NDOT=0d0
		   UDOT=0D0
		Endif
!
		do j=1, jmax
!
			if(KT.GE.0.and.KT.le.NN) then 
			   SVAR = SVAR + DBLE(SA(J,KT)*Y(J))  
			   NVAR = NVAR + DBLE(NA(J,KT)*Y(J))  
			   UVAR = UVAR + DBLE(UA(J,KT)*Y(J))  
			Endif
!
			if(KT.GE.1.and.KT.le.NN-1) then 
			   SDOT = SDOT + DBLE(((SA(J,KT+1)-SA(J,KT-1)))*Y(J)) 
			   NDOT = NDOT + DBLE(((NA(J,KT+1)-NA(J,KT-1)))*Y(J)) 
			   UDOT = UDOT + DBLE(((UA(J,KT+1)-UA(J,KT-1)))*Y(J)) 				   		  			   			   
			Endif 
		enddo
!			
!
! ---- Writes on the ***cumulative*** file
!
		        if(KT==0.or.KT==NN)        & 
			write(99,'(I4,2X,F7.3,3(F10.3,2X))') KT, DELTA*KT*1000., & 
			SVAR +  SDOT_GIA_CORR*DELTA*FLOAT(KT), &
			NVAR +  NDOT_GIA_CORR*DELTA*FLOAT(KT), &
			UVAR +  UDOT_GIA_CORR*DELTA*FLOAT(KT)  
			if(KT.GE.1.and.KT.le.NN-1) & 
			write(99,'(I4,2X,F7.3,6(F10.3,2X))') KT, DELTA*KT*1000., & 
			SVAR + SDOT_GIA_CORR*DELTA*FLOAT(KT), &
			NVAR + NDOT_GIA_CORR*DELTA*FLOAT(KT), &
			UVAR + UDOT_GIA_CORR*DELTA*FLOAT(KT), &   
			SDOT/2D0/DDELTA  + SDOT_GIA_CORR, & 
			NDOT/2D0/DDELTA  + NDOT_GIA_CORR, & 
			UDOT/2D0/DDELTA  + UDOT_GIA_CORR  
!
! ---- Writes on the ***individual*** file
!
		        if(KT==0.or.KT==NN)        & 
			write(200,'(I4,2X,F7.3,3(F10.3,2X))') KT, DELTA*KT*1000., & 
			SVAR +  SDOT_GIA_CORR*DELTA*FLOAT(KT), &
			NVAR +  NDOT_GIA_CORR*DELTA*FLOAT(KT), &
			UVAR +  UDOT_GIA_CORR*DELTA*FLOAT(KT)  
			if(KT.GE.1.and.KT.le.NN-1) & 
			write(200,'(I4,2X,F7.3,6(F10.3,2X))') KT, DELTA*KT*1000., & 
			SVAR +  SDOT_GIA_CORR*DELTA*FLOAT(KT), &
			NVAR +  NDOT_GIA_CORR*DELTA*FLOAT(KT), &
			UVAR +  UDOT_GIA_CORR*DELTA*FLOAT(KT), &     
			SDOT/2D0/DDELTA  + SDOT_GIA_CORR, & 
			NDOT/2D0/DDELTA  + NDOT_GIA_CORR, & 
			UDOT/2D0/DDELTA  + UDOT_GIA_CORR  
!
 100   CONTINUE
!
       CLOSE(200) 
!
!
  6   CONTINUE   ! --- End of Loop on the tide-gauges 	  
!
!
!     
      ELSEIF(TG_DATABASE_FORMAT/='1')   THEN
!
!
! ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... 
! ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... 
! ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... 
! ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... 
! ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... 
! ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... 
!
!
!
      ENDIF  !  Endif on the type of input file
!
!
!
!
!
! --- Writing format for type '0'  
201 format(a8,2x,a3,2x,a11,1x,a8,2x,a5,2x,f8.3,2x,f8.3,2x,f8.3,2x,f8.3,2x,f8.3,4x,a30) 
! --- Writing format for type '1'  
202  FORMAT (2(i5, 1x), 2(f10.4, 1x), 3(i5, 1x), F7.2, 1x, F7.2, 4x, 3(f8.3, 2x), A30)
!
!
! --- Reading format for type '0'
200    format(a8,1x,a3,1x,a3,2x,a11,1x,a8,1x,      & 
              a3,1x,a5,1x,a7,3x,a2,1x,a2,1x,a1,1x, & 
              a3,1x,a2,1x,a1,3x,a30)
! --- Reading format for type '1'
119  FORMAT (2(i5, 1x), 2(f10.4, 1x), 3(i5, 1x), F7.2, 1x, F7.2, 4x, A30)
!
 close(17)
 close(10)
 close(99)  
!
 END PROGRAM PSMSL
!
!
!
!
!
!
