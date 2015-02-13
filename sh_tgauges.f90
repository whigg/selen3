!
! This is program "SH_TGAUGES.F90" 
!
! Last modified GS 04-07-2008 "Intel port..."
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
! ***  Revised GS August 2010 - g95 - Double precision implementation
! === Revised GS & GG May 2011 - New format ('1') for the tide gauge analysis   
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
!  Computes the SHs at the sites of the PSMSL database. The Shs to degree
!  'lmax' are sequentially written on binary file "shtidegauges.bin". 
!  Information about the sea level trends PSMSL stations comes from file 
!  "rlr-trends.txt", obtained from the PSMSL web site. Current version dates 
!  January 22, 2007. 
!
! Input files:
!	- a tide-gauges database
!
! Output files:
!	- shtidegauges.bin
!
!
!
! INCLUDE "harmonics.f90"
  PROGRAM SHPSMSL
  IMPLICIT NONE 
  INCLUDE "data.inc"
  INTEGER I, J, N, NH, DOM
  INTEGER, PARAMETER :: MAXN=10000  ! A large integer    
  CHARACTER*200 RIGA    	    ! A given row 
  COMPLEX*16 Y(JMAX) 		    ! Shs
  CHARACTER*56 STR		    ! First part of the string 				
  CHARACTER*2 LAT1, LAT2            ! Latitide degrees and minutes 	
  CHARACTER*3 LON1	            ! Longitude degrees 
  CHARACTER*2 LON2 		    ! Longitude minutes 
  CHARACTER*30 NAME                 ! Name of the station 
  CHARACTER*1 WLAT, WLON            ! Where is the station (e, w, s, n)  
  REAL*8 LON, LAT 		    ! Longitude and latitude 
!
!
!
! --- SH file for tide gauges sites (output) 
       open(17,file='shtidegauges.bin',status='unknown',form='unformatted') 
!
! --- Counting the header lines (beginning with '#'), and the data lines
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
       write(*,*) '    - There are', n, 'sites in file ', trim(adjustl(TGAUGES_DATABASE)) 
!
! === Reading all lines to extract longitude and latitude of the PSMSL stations. 
!     Other information will be extracted by PSMSL.F. GS 05-10-2007 - 
!
! --- Opening the PSMSL file... 
       open(10,file=TGAUGES_DATABASE,status='unknown')
!
! --- Reading the header lines  
       do i=1, nh ; read(10,'(a100)')riga ; enddo  
!
       do 5 i=1, n
!
        IF     (TG_DATABASE_FORMAT=='0') then 
! 
! --- Reading one line   
             read  (10,200) str, lat1, lat2, wlat, lon1, lon2, wlon, name       
!
! --- Extracts lon and lat of the station from the texst strings... 
 	     call  find_lon_lat (lat1, lat2, wlat, lon1, lon2, wlon, lon, lat) 
!	     
!
        ELSEIF (TG_DATABASE_FORMAT=='1') then        ! New format as of May 2011 ****** 
!
!
! --- Reading one line      		     	    
             Read(10,119)  j, j, lon, lat, j, j, j
!
	ENDIF                                        ! New format as of May 2011 ******


!
! --- Computing the SHs at the station coordinates 
 	     call harmo(lmax, lon, lat, y) 
!
! --- Adjusts 
	     do j=1,jmax ; y(j)=y(j)*(2-dom(j)) ; enddo
!
! --- Reports on shtidegauges.bin 
 	     write(17) y
!
! --- End of Loop
 5 continue 
!
!
! --- Reading format(s) 
!
! traditional 
200  format(a56,3x, a2,1x,a2,1x,a1,1x,a3,1x,a2,1x,a1,3x,a30)              
!
! new
119  FORMAT (2(i5, 1x), 2(f10.4, 1x), 3(i5, 1x), F7.2, 1x, F7.2, 4x, A30)

!
!
 close(22)
 close(10) 
 close(17) 
!
!
 end PROGRAM SHPSMSL
!
!
!
!
!
!
