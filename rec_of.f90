! 
! This is program "REC_OF.F90"  
!
! Last modified GS 04-11-2008 "Intel port"
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! *** Reviewed GS & FC November 2009 - Porting under gfortran
! *** Revised GS July 2010 - g95 - Double precision implementation 
! *** Reviewed DM June 2011 - Dynamic memory allocation
! *** Modified DM June 2011 - Parallel execution
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
! ------------------------------------------------------------------
! This program rebuilds the ocean function from its SH coefficients 
! ------------------------------------------------------------------
!
!
! Input files:
!	- px.dat
!	- shof.dat
!       - sh.bin 
!       - anchor.tmp      [NEW]
!
! Output files:
!	- recof.dat
!
!
!INCLUDE "harmonics.f90"
 PROGRAM RECOF
 IMPLICIT NONE
 INCLUDE "data.inc"
 CHARACTER*30 JUNK
 INTEGER I, J, K, MJ, DOM
 INTEGER, ALLOCATABLE :: MM(:), ANC(:)
 INTEGER :: NANCH
 REAL*8, ALLOCATABLE :: ALF(:,:)
 REAL*8 RESH, IMSH  
 REAL*8, ALLOCATABLE :: LON(:), LAT(:)
 REAL*8, ALLOCATABLE :: REC_OF(:)  
 COMPLEX*16, ALLOCATABLE :: OC(:), LONG_TABLE(:,:)
!
!
!
!
!
! --- Getting info about the number of anchor pixels 
!
  open(1,file='anchor.tmp',status='old')
  read(1,*) nanch
  close(1)
!
! --- Allocate memory space
!
 ALLOCATE( MM(JMAX), ANC(NP) )
 ALLOCATE( ALF(JMAX,NANCH) )
 ALLOCATE( LON(NP), LAT(NP) )
 ALLOCATE( OC(JMAX), LONG_TABLE(0:LMAX,NP) )
 ALLOCATE( REC_OF(NP) )
!
! --- Pre-computing the degree 'm' corresponding to 'J'
!
 Write(*,*) '    - Pre-computing the harmonic order'
 do j=1, jmax 
        mm(j)=mj(j) 
 enddo	
!
!
!---- Reading the table of spherical harmonics from <<sh.bin>>
!
 Write(*,*) "    - Reading ALFs and TRIGs from file sh.bin"
 open(3,file='sh.bin',status='unknown',form='unformatted')  
 	read(3)ALF
 	read(3)LONG_TABLE
 Close(3) 
 DO J=1, JMAX 
	ALF(J,:)=ALF(J,:)*(2-DOM(J))
 ENDDO
!
!
!
!---- Reading the SH coefficients from <<sh.bin>>
!
 Write(*,*) "    - Reading the SH OF coefficients from file shof.dat"
 open(3,file='shof.dat',status='unknown')
 do j=1, jmax   
	read(3,*) k, resh, imsh
	oc(j)=cmplx(resh,imsh)
 enddo
 close(3)  
!
!
!
!---- Reading the pixels from file <<px-table.dat>>
!
 Write(*,*) "    - Reading the pixels from file px-table.dat"
 open(2,file='px-table.dat',status='unknown')
 do i=1, 4 
    read(2,'(a30)') JUNK 	
 enddo 
 do i=1, np 
    read(2,*) lon(i), lat(i), anc(i)
 end do
 close(2)
!
! 
! 
!---- Re-building the ocean function
!  
 write(*,*) '    - Re-building the ocean function ...' 	
! 
!
 open(9,file='recof.dat',status='unknown') 
!
 rec_of = 0.
!
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J)      &
!$OMP    SHARED(REC_OF,ALF,OC,LONG_TABLE,MM,ANC)   &
!$OMP    SCHEDULE(GUIDED)
 do i=1, np 
	do j=1, jmax 		 
 		rec_of(i) = rec_of(i) + ALF(j,anc(i))*real(oc(j)*long_table(mm(j),i))   
	enddo
 end do
!$OMP END PARALLEL DO
!
!
! --- Writing the reconstructed OF to rec_of.dat
!
 write(*,*) '    - Writing the reconstructed ocean function to file rec_of.dat ...'
 open(9,file='recof.dat',status='unknown') 
 do i=1,np
     Write(9,*) lon(i), lat(i), rec_of(i)
 end do
 close(9)
!
!
! --- Deallocate memory space
!
 DEALLOCATE( MM, ANC )
 DEALLOCATE( ALF )
 deallocate( lon )
 deallocate( lat )
 deallocate( rec_of )
 DEALLOCATE( OC, LONG_TABLE )
!
!
!
!
 end program RECOF
!
!
!
