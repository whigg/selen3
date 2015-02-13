!
! This is program "SH_OF.F90" 
!
! Last modified GS 04-11-2008  Intel port 
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
! *** Revised GS July 2010 - g95 - Double precision implementation
! *** Modified DM November 2011 - Dynamic memory allocation & SMP parallelization
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
! --------------------------------------------------------------------------- 
! Computes the harmonic coefficients of the ocean function up to degree lmax.
! Input files are the pixels file "px-table.dat" and the set of CSHs stored 
! in the binary file "sh.bin". The harmonic coefficients are reported on file 
! 'shof.dat' in the format j(l,m), real part, imaginary part. 
! ---------------------------------------------------------------------------
!
! Input files:
!	- px-table.dat 
!	- sh.bin
!       - anchor.tmp  [NEW]
!	
! Output files:
!	- shof.dat
!
!
!INCLUDE "harmonics.f90"
 PROGRAM OF 
 IMPLICIT NONE 
 INCLUDE "data.inc"
 CHARACTER*12 JUNK
 INTEGER :: NANCH
 INTEGER I, J, K, MJ, NWET, NDRY, NPIX     
 INTEGER, ALLOCATABLE :: MM(:), WET(:), ANC(:)     ! was MM(JMAX), WET(NP), ANC(NP) 
 REAL*8 LLON, LLAT
 REAL*8, ALLOCATABLE :: LONP(:), LATP(:)           ! was LONP(NP), LATP(NP)
 REAL*8, ALLOCATABLE :: ALF(:,:)                   ! was ALF(JMAX,NANCH)
 REAL*8 OCR, OCI
 COMPLEX*16, ALLOCATABLE :: LONG_TABLE(:,:), OC(:) ! was LONG_TABLE(0:LMAX,NP), OC(JMAX)  
!
!
!
! --- Getting info about the number of anchor pixels 
!
  open(1,file='anchor.tmp',status='old')
  read(1,*) nanch
  close(1)
!
!
!
! --- Allocate memory space
!
 allocate( mm(jmax) )
 allocate( oc(jmax) )
 allocate( wet(np), anc(np) )
 allocate( lonp(np), latp(np) )
 allocate( long_table(0:lmax,np) )
 allocate( alf(jmax,nanch) )
!
!
!
!
! --- Pre-computing the degree 'm' corresponding to 'J'
!
 Write(*,*) '    - Pre-computing the harmonic order'
 do j=1, jmax 
        mm(j)=mj(j) 
 enddo	
!
!
! --- Reading the ALFs table...  
 	write(*,*) '    - Reading ALFs & TRIGs from file sh.bin'
 	open(3,file='sh.bin',status='unknown',form='unformatted') 
		read(3)ALF
		read(3)LONG_TABLE
 	close(3) 
!
!
! --- Examining the pixels table & extracting information
!
	Open(1,file='px-table.dat',status='unknown') 
	Do i=1, 4 
		Read(1,'(a12)')junk
	Enddo
	nwet=0
	ndry=0
	npix=0
	do i=1, 2*np 
		read(1,*,end=10) lonp(i), latp(i), anc(i), k, wet(i) 
		npix=npix+1
		if(wet(i)==1)nwet=nwet+1
		if(wet(i)==0)ndry=ndry+1	
	enddo
10      Close(1) 
	Write(*,*) "    - The number of pixels in the px-table is ", npix 
!If(npix==np)      then 
!Write(*,*) "    - This is consistent with parameters in data.inc"	
!else
!Write(*,*) "    - This is NOT consistent with parameters in data.inc"	
!Endif		     
	Write(*,*) "    - The number of WET pixels in the px-table is ", nwet 
	Write(*,*) "    - The number of DRY pixels in the px-table is ", ndry 
!If(nwet+ndry==np) then 
!Write(*,*) "    - This is consistent with parameters in data.inc"	
!else
!Write(*,*) "    - This is NOT consistent with parameters in data.inc"	
!Endif		     
!
!	
! --- Computing the OF SH coefficients- I first evaluate the "Continent 
!     function" (CF) SH coefficients since this involves a smaller number 
!     of pixels. Then I transform them into the OF coefficients taking 
!     into account that the following relationship holds: CF + OF = 1.   
!	
 write(*,*) '    - Building the ocean function coefficients...'
!
 OC=0.
!
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J) SHARED(OC,WET,ALF,ANC,LONG_TABLE,MM,NPIX)
 do j=1,jmax 
! 
!	If(jmax>=10000) then 
!        	        If(mod(jmax,10000)==0) Write(*,*) "    - wnw.f: ", j, " harmonics of ", jmax	
!        elseif(jmax>=1000.and.jmax<10000) then 
!        	        If(mod(jmax,1000)==0)  Write(*,*) "    - wnw.f: ", j, " harmonics of ", jmax 	
!	                Endif
 	do i=1,np 
	   if(wet(i)==0) OC(J) = OC(J) + ALF(J,ANC(I))*CONJG(LONG_TABLE(MM(J),I))
 	enddo
 enddo 
!$OMP END PARALLEL DO
!
 OC=OC/FLOAT(NP)
!
! --- Writing the result on "sh-of.dat"
!
 open(1,file='shof.dat',status='unknown')  
 do j=1, jmax   
   	if(j==1) write(1,*) j, 1.- real(oc(j)), -aimag(oc(j))
   	if(j/=1) write(1,*) j,   - real(oc(j)), -aimag(oc(j))
 enddo
 close(1) 
!
!
! --- Deallocate memory space
!
 DEALLOCATE( MM, WET, ANC )
 DEALLOCATE( LONP, LATP )
 DEALLOCATE( ALF )
 DEALLOCATE( LONG_TABLE, OC )
!!	
!
!
 END PROGRAM OF 
!
!
!
!
!
!
