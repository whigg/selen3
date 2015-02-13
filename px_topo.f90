!INCLUDE "harmonics.f90"
PROGRAM PX_TOPO 
!
! 
! - - - - - - - - - - -
! FC & GS June 18 2008 
! - - - - - - - - - - - 
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
! *** Reviewed GS December 2010 - REAL*8 in COSDD, SINDD
! *** Parallelized DM December 2011
! Feb 2012: Retouched for minor issues 
! *** Modified DM April 2012 to handle ETOPO1 and ETOPO2
!     (dynamic allocation of LONE, LATE and TOPO arrays)
!
! -------------------------------------
! Program that pixelizes the ETOPO file
! -------------------------------------
!
! This program produces a icosahedron pixelization of current topography. 
! Output file 'px-topo.dat' gives pixels lon-lat, topo, 0/1 (ocean mask),
! which represents a present-day pixelization of the ocean function - 
!
! Input files: - etopo5.xyz (a topofile) 
!              - pxa.dat (a sorted pxfile)  
!
! Output files: - px-topo.dat (pixelized topography)
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
! SELEN is distributed in the /hope that it will be useful, but WITHOUT ANY 
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
! details. 
! 
! You should have received a copy of the GNU General Public License along 
! with SELEN.  If not, see <http://www.gnu.org/licenses/>.
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!
!
 IMPLICIT NONE 
 INCLUDE "data.inc"
!
 CHARACTER*20, PARAMETER :: PXFILE='pxa.dat'
 CHARACTER*10 HEADER
 INTEGER I, J, K, ANC, NIN, NTOPO 
 INTEGER, PARAMETER :: MAXEL=10**7
 INTEGER, PARAMETER :: RANGE=3
 REAL*8, PARAMETER :: ALFA = (180./PI)*SQRT(4./FLOAT(NP))
 REAL*8 COSTETA, TETA, LONP, LATP, ACCA, JUNK 
 REAL*8, ALLOCATABLE :: LONE(:), LATE(:), TOPO(:) 
 REAL*8 COSDD, SINDD
!
!
!
!
! --- Count the number of lines in TOPO file
!
Write(*,*) "    - px-topo.f90: Analyzing file: ", topo_file 
open(2,file=topo_file,status='old') 
ntopo=0
! 	do j=1, maxel 
	10 read(2,*,end=11) junk
	   ntopo=ntopo+1
!	enddo
    goto 10
11 close(2)
write(*,*) "    - File ", topo_file, " contains ", ntopo, " elements" 
!
! --- Allocate memory
!
allocate( lone(ntopo), late(ntopo), topo(ntopo) )
!
! --- Read TOPO file
!
Write(*,*) "    - Reading file: ", topo_file 
open(2,file=topo_file,status='old') 
 	do j=1, ntopo 
	read(2,*) lone(j), late(j), topo(j) 
	if(lone(j) <= 0.) lone(j) = 360. + lone(j) 
	if(mod(j,1000000)==0)& 
	Write(*,*) "    - ", j
	enddo
 close(2)
!write(*,*) "    - Read", ntopo, " elements from file: ", topo_file 
!
!
!
open(1,file=pxfile,status='unknown') 
!
!Write(*,*) "    - Reading file: ", pxfile 
Write(*,*) "    - Pixelizing the topography"             
do i=1, 4 ; read(1,'(a20)')HEADER ; enddo
!
open (3,file=pxtopo_file,status='unknown') 
Write(3,*) "---------------------------------------------------------"
Write(3,*) "  * Present day pixelized topography & ocean function *  "
Write(3,*) "    lon(deg), lat(deg), anchor, pixel#, wet(1/0), topo   "
Write(3,*) "---------------------------------------------------------"
! 
!
!
!
!---- Loop on Tegmark pixels 
!
do 100 i=1, np 
!
	read(1,*) lonp, latp, anc, k  
!
	if(mod(i,10000)==0) & 
	write(*,*) '    - ', i, 'of ', np 
! 
        nin = 0 
	acca = 0.0  	
!
!
!---- Loop on ETOPO pixels 
!
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(J,COSTETA,TETA) &
!$OMP    SHARED(NTOPO,LONE,LATE,LONP,LATP,ACCA,NIN,TOPO) &
!$OMP    SCHEDULE(GUIDED)
 	do j=1, ntopo 
!	
	if(abs(lone(j)-lonp).le.range*alfa.and.&
	   abs(late(j)-latp).le.range*alfa) then 		
!
	   costeta = cosdd(90.-latp)*cosdd(90.-late(j))+ & 
  		     sindd(90.-latp)*sindd(90.-late(j))*cosdd(lonp-lone(j))	
                     teta = acos(costeta)*180./pi 
!	
! If a topography pixels falls within the "radius" of 
! the current Tegmark pixel, update the topography.  
!   		
 		if(teta <= alfa) then 
!$OMP CRITICAL
	         	nin=nin+1 
	         	acca=acca + topo(j) 
!$OMP END CRITICAL
		endif 
        endif
!	
        enddo
!$OMP END PARALLEL DO
!
! Averaged topography within the "pixel radius" 
!
        acca=acca/float(nin)
!	 
	If(acca.ge.0.) then 
		write(3,*) lonp, latp, anc, k, '0', acca
		      else
		write(3,*) lonp, latp, anc, k, '1', acca
	endif
		         
!
100 continue
! 
    close(1) 
	close(3) 
!
    deallocate( lone, late, topo )
!
!
stop
!
end program PX_TOPO
!
!
!
