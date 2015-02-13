!
! This is program "SH.F90" 
! 
! Last change: GS April 11, 2008 "INTEL port" 
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!)
! *** Reviewed GS & FC November 2009 - Porting under gfortran  
! *** Revised GS July 2010 - g95 - Double precision implementation
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
! WARRANTY; without even the implied warranty of MERCHANTABlILITY or FITNESS 
! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
! details. 
! 
! You should have received a copy of the GNU General Public License along 
! with SELEN.  If not, see <http://www.gnu.org/licenses/>.
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! --------------------------------------------------------------------------- 
! This program computes the ALFs at the latitudes of the "anchor" pixels, and 
! the TRIGs at pixels longitudes for degrees from 0 to LMAX. ALFs and TRIGs 
! are written on the binary, unformatted file 'sh.bin'. ### GS 15-09-2007 ###
! --------------------------------------------------------------------------- 
!
! Input files:
!	- px-lat.dat 
!       - px-table.dat 
!
! Output files:
! 	- sh.bin 
!
!
 INCLUDE "harmonics.f90"
 PROGRAM SH
 IMPLICIT NONE 
 INCLUDE "data.inc"
 INTEGER, PARAMETER :: NANCH=float(NP)/7. 
 INTEGER, PARAMETER :: A_LARGE_INTEGER=10**6
 CHARACTER*12 JUNK 
 INTEGER I, J, K, L, SIZE_ALF, SIZE_LON, ISTAT, NA, IJUNK
 INTEGER SLOTS, LO(100), HI(100)
 REAL*8 LATA, LONG
 REAL*8 COSDD, SINDD
!
 REAL*8,     DIMENSION (:,:), ALLOCATABLE:: ALF
 COMPLEX*16, DIMENSION (:,:), ALLOCATABLE:: LONG_TABLE
!
!
!
  write(*,*)'    - Maximum degree is: ', lmax 
  write(*,*)'    - JMAX is:', (lmax+1)*(lmax+2)/2
  write(*,*)'    - Resolution is: ', res
  write(*,*)'    - The number of pixels is: ', NP 
!
! --- Getting info about the number of anchor pixels 
!
  open(2,file='px-lat.dat',status='unknown') 
  do i=1, 4 
  	read(2,*) JUNK 
  enddo
  na=0
  do 2 i=1, a_large_integer 
   	read(2,*,end=3) JUNK 
	na=na+1 
2 continue 
3 Write(*,*) "    - Found ", na, " anchor pixels in file px-lat.dat"
  close(2) 
!
!
!
!
! --- Output file for ALFs
  open(3,file='sh_alf.bin',status='unknown',form='unformatted') 
!
! --- Computing the ALF at the latitudes of anchor pixels 
  Write(*,*) "    - Pre-computing ALFs at latitudes of anchor pixels"
  open(2,file='px-lat.dat',status='unknown') 
  do i=1, 4 
  	read(2,*) JUNK 
  enddo
!
  size_alf=1000
!
  ALLOCATE(ALF(JMAX,size_alf))
!
  CALL HARMO_BREAKER (NA, size_alf, SLOTS, LO, HI)
!  
  write(*,*) '    - Number of slots: ', slots  
  write(*,*) '    - Size of the slots for ALFs: ', size_alf    
!
  do k=1, SLOTS 	
!
  	write(*,*) '    -', lo(k), hi(k)     
!
  	do i=lo(k), hi(k) 
  	         read(2,*,end=9) ijunk, lata         
       		 call PLMBAR_MOD (lmax, lata, ALF(:,i))        
        enddo
!
9 close(2)
!
  write(3) ALF
!
  Enddo ! On the slots 
! 
  DEALLOCATE(ALF) 
!
  CLOSE(3) 
!
!
!
!
!
!
! --- Output file for the longitudinal component 
  open(3,file='sh_lon.bin',status='unknown',form='unformatted') 
!
! --- Computing the TRIGs at pixels  
  Write(*,*) "    - Pre-computing TRIG functions at pixels"
  open(2,file='px-table.dat',status='unknown') 
  do i=1, 4 
  	read(2,*) JUNK 
  enddo
!
  size_lon=10000
!
  ALLOCATE(LONG_TABLE(0:LMAX,size_lon))
!
  CALL HARMO_BREAKER (NP, size_lon, SLOTS, LO, HI)
!
  write(*,*) '    - Number of slots: ', slots  
  write(*,*) '    - Size of the slots for ALFs: ', size_lon   
!
  do k=1, SLOTS 	
!
  write(*,*) '    -', lo(k), hi(k)     
!
  do i=lo(k), hi(k)   
!
  read(2,*) long           
  	do l=0, lmax
	long_table(l,i)=cmplx(cosdd(l*long),sindd(l*long))
  	enddo	
  enddo
  close(2) 
!  
  write(3) long_table
!
  Enddo  ! On the slots 
!
  close(3)
!
  DEALLOCATE(LONG_TABLE)
!
  write(*,*) '    - ALFs & TRIGs are written ...'
! 
  END PROGRAM SH
!
!
!
