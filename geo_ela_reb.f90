!
! This is program "GEO.F90" 
!
! Created by GS September 2008 for version 2.7 
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! *** Revised GS July 2010 - g95 - Double precision implementation
! ===> Created from GEO.F90 GS December 9, 2010. 
! === Revised GG May 13 2011 - New format for the elareb analysis    
! === Revised GG May 19 2011 - New format for the elareb analysis   
! >>> GS July 2012: 5(f10.5,1x) is now used for the output S, U, N data  
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
!  This program first computes the 4pi-normalized, complex spherical harmonics 
!  at geodetic sites listed in file GEODETIC_DATABASE. Subsequantly, it computes 
!  the 3D displacement (UP, NORTH, EAST), and the scalars S and N at the same 
!  sites. This set of five quantites specifies completely the "geodetic state" 
!  of the site. 
!
!  Input files:
!	-  A list of geodetic sites (e. g., "gps-sites.dat")
!	- 'shs.bin' 
!	- 'shu.bin' 
!       - 'shn.bin'
!       - 'shv.bin'
!  
!  Output files: 
!	- 'geodetic-predictions.dat'	
!
!
!INCLUDE "harmonics.f90"
 PROGRAM GEO
 IMPLICIT NONE 
 INCLUDE "data.inc"
 CHARACTER*22 CJUNK
 CHARACTER*80 NAME(NGEOD) 
 CHARACTER*40 LITHOC
 CHARACTER*20 DATE, TIMC            
 INTEGER I, J, DOM, CODE(NGEOD)  
 REAL*8 AJUNK, LONS(NGEOD), LATS(NGEOD)
 COMPLEX*16 CS(JMAX,0:NN), & 
            CU(JMAX,0:NN), & 
	    CN(JMAX,0:NN), & 
	    CV(JMAX,0:NN)
 COMPLEX*16            Y(JMAX,NGEOD), & 
             GRD_THETA_Y(JMAX,NGEOD), & 
	    GRD_LAMBDA_Y(JMAX,NGEOD)
 REAL*8 RATE_S(NGEOD),  & 
        RATE_N(NGEOD),  & 
	RATE_UP(NGEOD), & 
	RATE_NO(NGEOD), & 
	RATE_EA(NGEOD)        
!
 REAL*8, PARAMETER :: from_m_to_mm=1000.
!
 CHARACTER*3 XCODE(NGEOD)
!
! ---- Declarations for the file format '1' (new as of may 13, 2011) 
  INTEGER, PARAMETER:: large_integer=10000
  INTEGER  IJ, NH,    & 
           STAT_ID,   & 
	   ID(NGEOD), & 
	   JUNK1, JUNK2, JUNK3, & 
	   TMIN(NGEOD), TMAX(NGEOD), NVAL(NGEOD)  
  REAL*8 STAT_LON, STAT_LAT, & 
         XJUNK
  CHARACTER*30  STAT_NAME
  CHARACTER*200 ROW    
   
!
! --------
! --------
! --------
!
!
!
  If       (ELAREB_DATABASE_FORMAT=='0') then   ! ------ format '0'
!
  write(*,*) '    - Reading geodetic sites coordinates'
  OPEN(1,FILE=GEODETIC_DATABASE,STATUS='unknown')
!
  Do 333 i=1, ngeod
       read(1, '(a3)')     xcode(i) 
       read(1,'(a80)')     name(i) 
!
! ---  The first datum is LONGITUDE, the second is LATITUDE   
!
       read(1,*) lons(i), lats(i)  
       if(lons(i).le.0.) lons(i)=360.+lons(i)  
!		
       read(1,'(a22)')cjunk 	
!
333 continue
  Close(1)
!
  ELSEIF   (ELAREB_DATABASE_FORMAT=='1') then    ! ------ format '1'
!  
! --- Counting the data and the header lines  
!
  write(*,*) '    - Reading geodetic sites coordinates'
  OPEN(1,FILE=GEODETIC_DATABASE,STATUS='unknown')
!
  nh=0
  do 3313 j=1, large_integer
           read(1,'(a200)',end=4741) row 
  if(row(1:1)=='#') nh=nh+1 
  3313 continue
  4741 continue 
  Close(1) 
!   
! --- Reading the data and the header lines  
!
  OPEN(1,FILE=GEODETIC_DATABASE,STATUS='unknown')
  Do i=1, nh  
  	Read(1,'(a200)') row 
  Enddo
  Do 303 i=1, ngeod 
        Read(1,119)  IJ,  STAT_ID,          & 
	             STAT_LON, STAT_LAT,    & 
		     JUNK1, JUNK2, JUNK3,   & 
		     XJUNK, XJUNK,          & 
		     STAT_NAME
!
! --- Reading format for type '1'
119  FORMAT (2(i5, 1x), 2(f10.4, 1x), 3(i5, 1x), F7.2, 1x, F7.2, 4x, A30)
!
        lons(i)=STAT_LON
	lats(i)=STAT_LAT
	name(i)=STAT_NAME 
	id(i)  =STAT_ID
	tmin(i)=JUNK1
	tmax(i)=JUNK2
	nval(i)=JUNK3
	if(lons(i).le.0.) lons(i)=360.+lons(i)
!  
303 continue 
! 
  Close(1) 
!
  ENDIF                                          ! ------ format   
!
!
  write(*,*) '    - Computing the harmonics at the geodetic sites'
  do i=1, NGEOD    
       call harmo            (lmax, lons(i), lats(i),            y(:,i)) 
       call grad_theta_harmo (lmax, lons(i), lats(i),  grd_theta_y(:,i))
       call grad_lambda_harmo(lmax, lons(i), lats(i), grd_lambda_y(:,i))
  enddo
!
! 
! Write(*,*) '    - geo.f: Reading the harmonics from binary files'
 open (101,file='shs.bin',status='UNKNOWN',form='UNFORMATTED') 
 open (102,file='shu.bin',status='UNKNOWN',form='UNFORMATTED') 
 open (103,file='shn.bin',status='UNKNOWN',form='UNFORMATTED') 
 open (104,file='shv.bin',status='UNKNOWN',form='UNFORMATTED') 
 read (101) CS
 read (102) CU
 read (103) CN
 read (104) CV
 close(101) ; close(102) ; close(103) ; close(104) 
!
 do 1 i=1, ngeod 
!
       rate_s(i) =0.  ! sea level 
       rate_n(i) =0.  ! geoid undulation 
       rate_up(i)=0.  ! velocity 'up' 
       rate_no(i)=0.  ! velocity 'north' 
       rate_ea(i)=0.  ! velocity 'east'  
! 
      		do 2 j=1, jmax 
	       		rate_s (i) = rate_s (i) + (2.-dom(j))*real(cs(j,nn)*y(j,i))  
 	       		rate_n (i) = rate_n (i) + (2.-dom(j))*real(cn(j,nn)*y(j,i))  
 	       		rate_up(i) = rate_up(i) + (2.-dom(j))*real(cu(j,nn)*y(j,i))  
 	       		rate_no(i) = rate_no(i) + (2.-dom(j))*real(cv(j,nn)*grd_theta_y(j,i))  
 	       		rate_ea(i) = rate_ea(i) + (2.-dom(j))*real(cv(j,nn)*grd_lambda_y(j,i)) 
 2    		continue 
!
      		rate_no(i) = -rate_no(i)
!
 1    continue 
!
!
!
!------- File for scattered geodetic predictions 
 open(82,file='geodetic-predictions.dat',status='unknown')
 call DATE_AND_TIME (date,timc) 
WRITE(82,*)'# SELEN 3.2'     
WRITE(82,*)'# ', date(1:4), '.', date(5:6), '.', date(7:8), & 
           ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6)
WRITE(82,*)"# Ice model: ", trim(adjustl(ice_model))
WRITE(82,*)"# Geodetic database:   ", GEODETIC_DATABASE 
WRITE(82,*)"# Model code (see TABOO User guide): ", CDE
WRITE(82,*)"# SLE mode of solution: ", IMODE
WRITE(82,*)"# SLE iterations: ", SMAX 
WRITE(82,*)"# Maximum harmonic degree: ", LMAX 
WRITE(82,*)"# Tegmark resolution: ", RES 
WRITE(82,*)"#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "
WRITE(82,*)"# 	             /// Geodetic predictions for all the sites in database ///                              "
WRITE(82,*)"#                                                                                                        "
WRITE(82,*)"#  Code   lon       lat      ------ all velocities are in units of mm/yr ------                          "
WRITE(82,*)"#         deg       deg       UP         NORTH      EAST       S-dot      N-dot     Tmin  Tmax  Nval Site"  			     
WRITE(82,*)"#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "
WRITE(82,*)"#	 "
!
!
  do i=1, ngeod 	
!
  IF           (ELAREB_DATABASE_FORMAT=='0')  THEN 
!	
    Write(82,'(3x,a3,2(1x,f8.2),5(f10.5,1x),4x,a80)') & 
		  xcode(i), lons(i), lats(i), rate_up(i) * from_m_to_mm, & 
		                             rate_no(i) * from_m_to_mm, & 
					     rate_ea(i) * from_m_to_mm, & 
					     rate_s (i) * from_m_to_mm, & 
					     rate_n (i) * from_m_to_mm, & 
					     name(i)

  ELSEIF       (ELAREB_DATABASE_FORMAT=='1')  THEN 
!
    Write(82,'(3x,i4,2(1x,f8.2), 5(f10.5,1x),4x, 3(i5, 1x), a80)') & 
		  id(i),   lons(i), lats(i),     & 
		                             rate_up(i) * from_m_to_mm, & 
		                             rate_no(i) * from_m_to_mm, & 
					     rate_ea(i) * from_m_to_mm, & 
					     rate_s (i) * from_m_to_mm, & 
					     rate_n (i) * from_m_to_mm, & 
					     TMIN(i), TMAX(i), NVAL(i),     & 
					     name(i)
!
  ENDIF
!
  enddo
!
  close(82)
!
!
  END PROGRAM GEO
!
!
!
!
!
!
