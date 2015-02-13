!
! This is program "TGAUGES.F90" 
! 
! Last modified GS 04-11-2008 "Intel port"
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
! Reviewed again by GS on April 2010 Porting under g95 
! *** Revised GS July 2010 - g95 - Double precision implementation
! === Revised GS & GG May 2011 - New format ('1') for the tide gauge analysis   
! === Revised GS & FC May 21 2011 - DELTA parameter in ice history     
! Feb 2012: Implementation of the numerical derivative "in the future"    
! Apr 2012: The format '1' had a misprint in the header. Corrected  while in the Hospital... 
! April 2012: Added a "#" in the first column of the header lines (while in the Hospital)
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
!
  REAL*8, PARAMETER :: DDELTA=DELTA 
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
	OPEN (99,FILE='ptidegauges.dat',STATUS='unknown')
! 
!
! --- Header information about the SELEN settings 
!
!
        call DATE_AND_TIME (date,timc)      
        Write(99,*) date(1:4), '.', date(5:6), '.', date(7:8), & 
	          ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6)
	Write(99,*) "Ice model: ", trim(adjustl(ice_model))
	If(NV.NE.1.AND.CDE.NE.-1) THEN
		Write(99,*) "Number of mantle layers: ", NV 
		Write(99,*) "Model code (see TABOO User guide): ", CDE
	ELSE
		Write(99,*) "Computations by ALMA - See log files for model details" 
	ENDIF 
	Write(99,*) "Viscosity model: ", visco_model
!
	open(32,file=visco_model,status='unknown') 
	do i=1,NV+1 
	read(32,'(a100)')riga 
	If(NV.NE.1.AND.CDE.NE.-1) THEN
	if(i==1)  then 
		Write(99,*) "Thickness of the lithosphere (km): ", trim(adjustl(riga))
		  else
		Write(99,*) "Viscosity (bottom-to-top, Haskell units) ", trim(adjustl(riga))
	          Endif	
	ENDIF	
	enddo ; close(32) 
	Write (99,*) "SLE iterations: ", SMAX 
	Write (99,*) "SLE mode of solution: ", IMODE
	Write (99,*) "Maximum harmonic degree: ", LMAX 
	Write (99,*) "Tegmark resolution: ", RES 
!
	If   (TG_DATABASE_FORMAT=='0')  THEN   ! ----------------  Header info for 'type 0' TG input files... 
!
	WRITE (99,*) ''
	WRITE (99,*) &
	' PSMSL    yrs    range       trend  error     lon       lat       s-dot     u-dot     n-dot    PSMSL station name' 
	WRITE (99,*) &
	' code             yrs            mm/yr        deg       deg       mm/yr     mm/yr     mm/yr    =================='
	WRITE (99,*) ''
!
        ELSEIF(TG_DATABASE_FORMAT=='1') THEN   ! ----------------  Header info for 'type 1' TG input files... 
	WRITE (99,*) ''
	WRITE (99,*) &
        '   #   Id   Longitude   Latitude  tmin  tmay  nval   trend   error       s-dot     u-dot     n-dot  PSMSL station name'
	WRITE (99,*) &
        '            (degrees)  (degrees)  ---- (years) ---   -- (mm/yr) --       -------- (mm/yr) --------  =================='
	WRITE (99,*) ''
!
	ENDIF                                 ! ----------------  End of database-sensitive header info ... 
!
!
! --- tgauges.f: reading sealevel CSH coefficients from file', filename 
       open(3,file='shs.bin',status='unknown',form='unformatted') ; read(3) SA ; close(3)
       open(3,file='shn.bin',status='unknown',form='unformatted') ; read(3) NA ; close(3)
       open(3,file='shu.bin',status='unknown',form='unformatted') ; read(3) UA ; close(3)
!
! --- Counting the header lines (, beginning with '#'), and the data lines
       open(10,file=TGAUGES_DATABASE,status='unknown')
       nh=0 
       n=0
       do i=1, maxn 
       		read(10,'(a100)',end=1) riga
       		if(riga(1:1)=='#') nh=nh+1
       		if(riga(1:1)/='#') n=n+1
       enddo 
 1     close(10)  
!
! --- Opening the PSMSL data file... 
       open(10,file=TGAUGES_DATABASE,status='unknown')
!
! --- Reading again the header lines  
       do i=1, nh ; read(10,'(a100)',end=1) riga ; enddo  
!
! === Reading all lines. GS = October 5 2007 - 
!
       write(*,*) '    - Rate of sealevel change at site'
!



!
! *************************************************************************************************
! *************************************************************************************************
!
       If(TG_DATABASE_FORMAT=='0')       THEN   !  -------- Input TG file has format '0'
!
       do 5 i=1, n 
!
       		read  (10,200) pcode, gloss_code, years, range_years, & 
   		      	       datum , plus_minus, error, stdv, & 
			       lat1, lat2, wlat, lon1, lon2, wlon, name        
!
! --- Extracts lon and lat of the station from the texst strings... 
                call  find_lon_lat (lat1, lat2, wlat, lon1, lon2, wlon, lon, lat) 
!
! --- SH st the current PSMSL site...
      	        read(17) y
!    
! --- Predicted rate at the site 'i' 
		If(n>=50)                then 
      	        If(mod(i,50)==0.or.i==1) write(*,*) '    - ', i, 'of', n, name 
		elseif(10<=n.and.n<=50)  then
      	        If(mod(i,10)==0.or.i==1) write(*,*) '    - ', i, 'of', n, name 		
		elseif(n<10)             then
		write(*,*) '    - ', i, 'of', n, name 			       
					 endif 
!
! --- Computing Sdot, Udot and Ndot at tide gauge		
! --- Computing Sdot, Udot and Ndot at tide gauge		
		PSDOT = 0.
		PNDOT = 0. 
		PUDOT = 0. 
		do j=1, jmax
                IF(IDER==1)     THEN    ! Updated on February 2012 
			PSDOT = PSDOT + real(((sa(j,nn+1)-sa(j,nn-1))/DDELTA)*y(j)) 
   		        PNDOT = PNDOT + real(((na(j,nn+1)-na(j,nn-1))/DDELTA)*y(j)) 
			PUDOT = PUDOT + real(((ua(j,nn+1)-ua(j,nn-1))/DDELTA)*y(j)) 
		ELSEIF(IDER==2) THEN 
			PSDOT = PSDOT + real(((sa(j,nn+1)-sa(j,nn-1))/2D0/DDELTA)*y(j)) 
   		        PNDOT = PNDOT + real(((na(j,nn+1)-na(j,nn-1))/2D0/DDELTA)*y(j)) 
			PUDOT = PUDOT + real(((ua(j,nn+1)-ua(j,nn-1))/2D0/DDELTA)*y(j)) 
                ENDIF
		enddo
!
  		WRITE (99,201)  pcode, & 
				years, & 
				range_years, & 
				datum, & 
				error, & 
				lon, lat, & 
				psdot, pudot, pndot, & 
				name 
!
 5    continue 
!
	ELSEIF(TG_DATABASE_FORMAT=='1')   THEN  !  ------- Input TG file has format '1'
!
        Do 6 i=1, n 
!
! --- Reading one line      		     	    
        Read(10,119)  IJ, STAT_ID, STAT_LON, STAT_LAT, & 
		          TMIN,    TMAX,  NVAL,        & 
		          TREND,    DTREND,	       & 
		          STAT_NAME

119  FORMAT (2(i5, 1x), 2(f10.4, 1x), 3(i5, 1x), F7.2, 1x, F7.2, 4x, A30)

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
! --- Computing Sdot, Udot and Ndot at tide gauge		
		PSDOT = 0.
		PNDOT = 0. 
		PUDOT = 0. 
		do j=1, jmax
                IF(IDER==1)     THEN    ! Updated on February 2012 
			PSDOT = PSDOT + real(((sa(j,nn)-sa(j,nn-1))/DDELTA)*y(j)) 
   		        PNDOT = PNDOT + real(((na(j,nn)-na(j,nn-1))/DDELTA)*y(j)) 
			PUDOT = PUDOT + real(((ua(j,nn)-ua(j,nn-1))/DDELTA)*y(j)) 
		ELSEIF(IDER==2) THEN 
			PSDOT = PSDOT + real(((sa(j,nn+1)-sa(j,nn-1))/2D0/DDELTA)*y(j)) 
   		        PNDOT = PNDOT + real(((na(j,nn+1)-na(j,nn-1))/2D0/DDELTA)*y(j)) 
			PUDOT = PUDOT + real(((ua(j,nn+1)-ua(j,nn-1))/2D0/DDELTA)*y(j)) 
                ENDIF
		enddo
!
        WRITE(99,202)  IJ, STAT_ID, STAT_LON, STAT_LAT, & 
		           TMIN,    TMAX,  NVAL,        & 
		           TREND,    DTREND,	        & 		   
                           PSDOT, PUDOT, PNDOT,         & 		  			   
		           STAT_NAME
!
! --- End of Loop
  6     continue 	  
!
!
        ENDIF                              !  -----------------    Endif on the type of input file


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
!119  FORMAT (2(i5, 1x), 2(f10.4, 1x), 3(i5, 1x), F7.2, 1x, F7.2, 4x, A30)
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
