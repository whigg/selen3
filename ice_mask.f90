! 
!*****************************
! This is program ICE_MASK.F90 
!*****************************
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
!INCLUDE "harmonics.f90"
 PROGRAM ICE_FILTER 
!
! -----------------------------------------------------------
! Program that pixelizes the ice model and creates ICE MASKS
! -----------------------------------------------------------
!
! Input files: - icefile (e. g., ice5g.dat)
!              - pxfile (e. g., pxa.dat)  
!
! Output files: - ice mask files (imask-XX.dat), XX=0, ... NN+2 
!
!
! FC & GS June 18 2008
! FC & GS July 05 2009  
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
! Revised on April 2010 by GS for KL ice and g95 implementation
! *** Revised GS July 2010 - g95 - Double precision implementation
! *** Revised GS December 2010 - implementation of ice model ICE3G
! *** Revised DM 2012 -- Loop parallelization
! Feb 2012: Implementation of the numerical derivative "in the future" 
! Feb 2015 -- Removed the derivative in the future (DM) 
!
! BUT THIS MUST BE RETHINKED 
!
!
 IMPLICIT NONE 
 INCLUDE "data.inc"
 INTEGER I, J, K, IJ, IOC, MUL, NIN, ILAT, NICE 
 INTEGER, PARAMETER :: RANGE=5
 REAL*8, PARAMETER :: ALFA = (180./PI)*SQRT(4./FLOAT(NP))
 REAL*8 COSTETA, TETA 
 REAL*8, ALLOCATABLE :: LONP(:), LATP(:), ICETHICK(:,:)
 REAL*8, ALLOCATABLE :: ICE(:,:) 
 REAL*8 AJ
 CHARACTER*20, PARAMETER :: PXFILE='pxa.dat'
 CHARACTER*20 HEADER
 CHARACTER*2 LABCHAR 
 REAL*8 COSDD, SINDD
!
!
!#---- Declarations for any ice model 
       REAL*8 LONI(NEL), LATI(NEL)
!        
!#---- Declarations for "ICE5"
       INTEGER, PARAMETER :: NH5=24
!        
!#---- Declarations for "ICE3"
       INTEGER, PARAMETER :: NH3=20
!
!#---- Declarations for "ANU05"
       INTEGER, PARAMETER :: NHU=26 
       REAL HAMP(NEL)

!
! Input files:  
!    - A pixel centered ice model (e. g., ice5g.dat) 
!    - An icosahedron pixelization (e. g., pxa.dat) 	   
!
! Output files: 
!    - "imask-XX.dat", mask files for ice/topography 
!    - "XX" is time from the LGM
!
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!   End of Declarations       End of Declarations       End of Declarations  
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!
!
! --- Allocate memory 
!
 ALLOCATE( LONP(NP), LATP(NP), ICETHICK(NP,0:NN+1) )
 ALLOCATE( ICE(NEL,0:NN+1) )
!
 ICE(:,:)     =0.
 ICETHICK(:,:)=0.  
!
!
! IF   (ICE_MODEL(1:4)/='ice5'.and.& 
!       ICE_MODEL(1:4)/='anu0'.and.& 
!       ICE_MODEL(1:4)/='ice3') Then  
!
 IF   (ICE_MODEL(1:4)/='ice5'.and. & 
       ICE_MODEL(1:4)/='anu0') Then  
  write(88,*) "    - Masking  (hence time-evolving shorelines)  is currentlt available "  
  write(88,*) "    - only for ice5*, & anu05*  of models - The program will stop "  
  write(*, *) "    - Masking  (hence time-evolving shorelines)  is currentlt available "  
  write(*, *) "    - only for ice5*, & anu05*  of models - The program will stop " 
 CALL STOP_CONFIG 
 Else
  write(88,*) "    - Masking the ice distribution ..."  
  write(*, *) "    - Masking the ice distribution ..."   
 Endif 
!
 IF   (ICE_MODEL(1:4)=='ice5') THEN 
!
!---- Reading the ICE5G ice file 
!
 open(2,file=ice_model,status='unknown') 
 do j=1, NH5 
	read(2,'(a20)')header 
 enddo
 nice=0
 do j=1, nel 
     Read(2,*) ij, loni(j), ilat, lati(j), (ice(j,k), k=0, nn) 
     nice=nice+1
 enddo
!
 ice(:,nn+1)=ice(:,nn)
!
 close(2) 
!
!
 ELSEIF(ICE_MODEL(1:4)=='anu0') THEN 
!
!---- Reading the ANU ice file 
! 
 open(2,file=ice_model,status='unknown') 
 do j=1, NHU 
	read(2,'(a20)')header 
 enddo
 nice=0
 do j=1, nel  
     Read(2,*) ij, loni(j), lati(j), hamp(j), (ice(j,k), k=0, nn)
 nice=nice+1
 enddo
!
 ice(:,nn+1)=ice(:,nn)
!
 close(2) 
!
 ELSEIF(ICE_MODEL(1:4)=='ice3') THEN 
!
!---- Reading the ICE3G ice file 
! 
 open(2,file=ice_model,status='unknown') 
 do j=1, NH3 
	read(2,'(a20)')header 
 enddo
 nice=0
 do j=1, nel  
     Read(2,111) ij, lati(j), loni(j), hamp(j), aj, (ice(j,k), k=0, 7)
     read(2,112)                                    (ice(j,k), k=8, nn)
     lati(j)=90.-lati(j)
 nice=nice+1
!#---- Formats for ICE3G (formerly provided by JXM) 
111    format(i4,2f9.3,2f6.3,1x,8f5.0)
112    format(15f5.0)
 enddo
!
 ice(:,nn+1)=ice(:,nn)
!
 close(2)
! 
 ENDIF
!
!
!---- A check on the ice elements 
!
 If(nel==nice) then 
 	write(*,*) '    - Read ', nice, '  elements from ', ice_model
 else
 	write(88,*) 'ICE_MASK.F90: Problem encountered reading ice model: ', ice_model
 	write(*, *) 'ICE_MASK.F90: Problem encountered reading ice model: ', ice_model
	CALL STOP_CONFIG
 	write(88,*) 'ICE_MASK.F90: Problem encountered reading ice model: ', ice_model
 Endif 
!
!
!
!---- Reading the pixelization file
!
 open(1,file=pxfile,status='unknown') 
 do i=1, 4 
	read(1,'(a20)')HEADER 
 enddo
 do i=1, np  
	read(1,*) lonp(i), latp(i) 
 enddo
 close(1) 
!
!
! Time loop 
! - - - - -  
!DO 300 K=0, NN+2     ! Updated on February 2012 (it was "nn+1")  
 DO 300 K=0, NN+1     ! Reverted back to nn+1 -- CHECK
!
! open(3,file='junk.dat',status='unknown') 
! if(k<=9) write(3,'(a1,i1)') '0',k  
! if(k> 9) write(3,'(i2)')        k
! close(3)
! open(3,file='junk.dat',status='unknown')
! read(3,'(a2)') labchar
! close(3)
 write(labchar,'(i2.2)') k
!
!
 if(k==0.or.k==nn+2) write(*,*) "    - step ", labchar 
!
 open(18,file='imask-'//labchar//'.dat',status='unknown') 
!
!
! Loop on the pixels 
! - - - - - - - - - -
 DO 400 I=1, NP 
!
 icethick(i,k)=0.
!
 nin=0
!
! Loop on the ice elements 
! - - - - - - - - - - - - 
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP    SHARED(NICE,LONI,LONP,LATI,LATP,I,NIN,ICETHICK,K,ICE) &
!$OMP    PRIVATE(J,COSTETA,TETA) &
!$OMP    SCHEDULE(GUIDED)
 DO j=1, nice 
!
!
! If we are relatively close to the i-th pixel...
!
 If    (abs(loni(j)-lonp(i)).le.range*alfa.and. & 
        abs(lati(j)-latp(i)).le.range*alfa) then 

 	costeta = cosdd(90.-latp(i))*cosdd(90.-lati(j))+ & 
	          sindd(90.-latp(i))*sindd(90.-lati(j))*cosdd(lonp(i)-loni(j))			  	 
 	teta = acos(costeta)*180./PI	
!	
! ... and if an ice element falls within the "radius" of 
! the current Tegmark pixel, update ice thickness.  
! 		
 	if(teta <= alfa) then 
!$OMP CRITICAL	
 	        nin=nin+1 
  		icethick(i,k) = icethick(i,k) + ice(j,k) 
!$OMP END CRITICAL
	endif	
 endif
!	
 END DO
!$OMP END PARALLEL DO
!
! --- Averaged ice thickness within the "pixel radius". Masks bear no information 
!      on the ice thickness. This is left for possible future purposes. 
!
      icethick(i,k) = icethick(i,k)/float(nin) 
!
! Ice is "1", no ice is "0"
!
 if(icethick(i,k)>0)then 
		write(18,*)lonp(i), latp(i), '1', icethick(i,k), icethick(i,k)   	    
	            else
		write(18,*)lonp(i), latp(i), '0', ' 0.0 ',       'NaN'   
 endif	       	
!
 400 continue 
!
 close(18) 
!
 300 continue 
!
!
 DEALLOCATE( LONP, LATP, ICETHICK )
 DEALLOCATE( ICE )
!
!
 Stop
!
END PROGRAM ICE_FILTER 
!
!
!
!
!
!
