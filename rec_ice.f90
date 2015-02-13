!
! REC_ICE.F90
!
! Last modified GS 04-11-2008 "Intel port"
! Touched on July 11, for version 2.6 
! Re-touched July 13, N & S poles 
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
! *** Revised GS July 2010 - g95 - Double precision implementation 
! Modified by GS on June 25 for the multi-stepping for the ice2sea computations 
! *** Revised DM November 2012 - Dynamic memory allocation & loop parallelism
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
! -----------------------------------------------------------  
! This program rebuilds the ice sheets form the coefficients
! -----------------------------------------------------------
!
! Input files:
!	- px-table.dat
!	- sht*.dat
!       - sh.bin 
!
! Output files:
!	- rect*.dat 
!
!
!INCLUDE "harmonics.f90"
 PROGRAM R
 IMPLICIT NONE
 INCLUDE "data.inc"
 CHARACTER*13, PARAMETER :: FMT='(3(f16.8,1X))'
 CHARACTER*22 INPUT_FILENAME, OUTPUT_FILENAME 
 CHARACTER*30 HEADER
 CHARACTER*3 LABCHAR 
 INTEGER I, J, N, K, L, MJ, DOM, JUNK
 INTEGER, ALLOCATABLE :: MM(:)   
 INTEGER :: NANCH
 INTEGER, ALLOCATABLE :: ANC(:) 
!
 REAL*8, ALLOCATABLE :: ALF(:,:)
 REAL*8 RESH, IMSH   
 REAL*8, ALLOCATABLE :: LON(:), LAT(:), REC(:)
 REAL*8 LONX, LATX, ARG
 COMPLEX*16, ALLOCATABLE :: SHI(:), LONG_TABLE(:,:)
 REAL*8, PARAMETER :: EPS = 0.01
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
  allocate( mm(jmax) )
  allocate( anc(np) )
  allocate( alf(jmax,nanch) )
  allocate( lon(np), lat(np), rec(np) )
  allocate( shi(jmax) )
  allocate( long_table(0:lmax,np) )  
!
! --- Pre-computing the degree 'm' corresponding to 'J'
!
! Write(*,*) "    - rec_ice: Pre-computing the harmonic order"
 do j=1, jmax 
        mm(j)=mj(j) 
 enddo	
!
!
!---- Reading the table of spherical harmonics from <<sh.bin>>
!
! Write(*,*) "    - rec_ice: Reading ALFs & TRIGs from <<sh.bin>>"
 open(3,file='sh.bin',status='unknown',form='unformatted')  
 	read(3)ALF
 	read(3)LONG_TABLE
 Close(3) 
 DO J=1, JMAX 
	ALF(J,:)=ALF(J,:)*(2-DOM(J))
 ENDDO 
!
!
!---- Reading the pixels table file... 
!
! Write(*,*) "    - rec_ice: Reading file <<px-table.dat>>"
 open(2,file='px-table.dat',status='unknown')
 do i=1, 4 
    read(2,'(a30)') HEADER 	
 enddo 
 do i=1, np
    read(2,*) lon(i), lat(i), anc(i)
 enddo 
 close(2)  
!
!
write(*,*) '    - Ice reconstruction at all steps'
!  
do k=0,NN+1  
!
        write(labchar,'(i3.3)') k
!          
	if(k==0.or.k==nn+1)&
	write(*,*) '    - rec_ice: Step ', labchar, ' of ',NN+1 
!
	input_filename='sht'//labchar//'.dat'
!
	open(1,file=input_filename,status='unknown') 
!
	do j=1,jmax 
		read(1,*)  junk, resh, imsh
		shi(j)=CMPLX(resh,imsh)
	enddo
	close(1) 
!
!
!
    rec = 0.
!
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J) &
!$OMP   SHARED(REC,ALF,ANC,SHI,LONG_TABLE,MM) &
!$OMP   SCHEDULE(GUIDED)
	do i=1, np 	
		do j=1, jmax 
 			rec(i) = rec(i) + & 
			ALF(j,anc(i))*real(shi(j)*LONG_TABLE(MM(J),I))  
		enddo
    enddo
!$OMP END PARALLEL DO
!
!
!
!
	output_filename='rect'//labchar//'.dat'
!
	open(109,file=output_filename,status='unknown') 
!
!
	do i=1, np 	
!		
!
!
! --- South pole
        If     (lat(i)==-90.) then  
       		latx=lat(i)+eps
       		do l=1, 4 
			lonx=(l-1)*90.   		
            		write(109,fmt) lonx, latx, rec(i)
       		enddo
!
! --- North pole
        elseif (lat(i)==+90.) then   
       		latx=lat(i)-eps
       		do l=1, 4 
			lonx=(l-1)*90.   			
            		write(109,fmt) lonx, latx, rec(i) 
       		enddo
!
!
! --- Elsewhere 
        elseif (lat(i)/=-90.and.lat(i)/=90.) then 
       		latx=lat(i) 
		lonx=lon(i) 
            		write(109,fmt) lonx, latx, rec(i)
        Endif	       
!
    end do
!
	close(109) 
!
  end do
!
write(*,*) '    - Output written on files rect*.dat ' 
!
!
!
 deallocate( mm )
 deallocate( anc )
 deallocate( alf )
 deallocate( lon, lat, rec )
 deallocate( shi )
 deallocate( long_table )
!
!
!
end program R 
!
!
!
