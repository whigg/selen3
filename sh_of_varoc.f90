!
! This is program "SH_OF_VAROC.F90"  
!
! Last modified GS 04-11-2008  Intel port 
! Retouched June 2008 [varying coastlines]  
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
! Revised by GS April 25, 2010 - Cleaning up
! *** Revised GS July 2010 - g95 - Double precision implementation
! *** Loop parallelization by DM December 2011
! Feb 2012: Implementation of the numerical derivative "in the future"
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
!	- px-table-XX.dat 
!	- sh.bin
!       - anchor.tmp   [NEW]
!	
! Output files:
!	- shof-XX.dat
!
!
!INCLUDE "harmonics.f90"
 PROGRAM SH_OF_VAROC 
 IMPLICIT NONE 
 INCLUDE "data.inc"
 CHARACTER*12 JUNK
 CHARACTER*2 LABCHAR 
 INTEGER :: NANCH
 INTEGER I, J, K, KT, MJ, NWET, NDRY, NPIX
 INTEGER, ALLOCATABLE :: MM(:), WET(:), ANC(:)
 REAL*8, ALLOCATABLE :: LONP(:), LATP(:) 
 REAL*8 LLON, LLAT
 REAL*8, ALLOCATABLE :: ALF(:,:)
 REAL*8 OCR, OCI
 COMPLEX*16, ALLOCATABLE :: LONG_TABLE(:,:), OC(:,:)  
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
! Updated Feb. 2012 (NB: nn+1 as before, but of(nn+1) is NOT imposed to of(nn))
!
  ALLOCATE( MM(JMAX), WET(NP), ANC(NP) )
  ALLOCATE( LONP(NP), LATP(NP) )
  ALLOCATE( ALF(JMAX,NANCH) )
  ALLOCATE( LONG_TABLE(0:LMAX,NP), OC(JMAX,0:NN+1) ) 
!
! --- Pre-computing the degree 'm' corresponding to 'J'
!
 	do j=1, jmax 
        	mm(j)=mj(j) 
 	enddo	
!
!
! --- Reading the ALFs table...  
	write(*,*) '    - sh_of_varoc.f: Reading ALFs & TRIGs from file <<sh.bin>>'
 	open(3,file='sh.bin',status='unknown',form='unformatted') 
		read(3)ALF
		read(3)LONG_TABLE
 	close(3) 
!
!
!+++++++++++++++++++++++++++++++++++++++++++++++++
	do 300 kt = 0, nn+1   ! Updated Feb. 2012    
!+++++++++++++++++++++++++++++++++++++++++++++++++
!
!	open(3,file='junk.dat',status='unknown') 
!	if(kt<=9) write(3,'(a1,i1)') '0',kt  
!	if(kt> 9) write(3,'(i2)')        kt
!	close(3)
!       open(3,file='junk.dat',status='unknown')
!	read(3,'(a2)') labchar
!	close(3)
        write(labchar,'(i2.2)') kt
!
	if(kt==0.or.kt==nn+1) write(*,*) "    - OF SH harmonics for step ", labchar 
!
!
! --- Examining the pixels table & extracting information
!
	Open(17,file='px-table-'//labchar//'.dat',status='unknown') 
	Do i=1, 4 
		Read(17,'(a12)')junk
	Enddo
	nwet=0
	ndry=0
	npix=0
	do i=1, 2*np 
		read(17,*,end=10) lonp(i), latp(i), anc(i), k, wet(i) 
		npix=npix+1
		if(wet(i)==1)nwet=nwet+1
		if(wet(i)==0)ndry=ndry+1	
	enddo
10      Close(17) 
!
!	
! --- Computing the OF SH coefficients- I first evaluate the "Continent 
!     function" (CF) SH coefficients since this involves a smaller number 
!     of pixels. Then I transform them into the OF coefficients taking 
!     into account that the following relationship holds: CF + OF = 1.   
!
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP      PRIVATE(I,J) SHARED(KT,OC,WET,ALF,ANC,LONG_TABLE,MM)
 do j=1,jmax 
! 
!	If(jmax>=10000) then 
!        	        If(mod(jmax,10000)==0) & 
!			Write(*,*) "    - sh_of_varoc.f: ", j, " harmonics of ", jmax	
!        elseif(jmax>=1000.and.jmax<10000) then 
!        	        If(mod(jmax,1000)==0)  & 
!			Write(*,*) "    - sh_of_varoc.f: ", j, " harmonics of ", jmax 	
!	                Endif
	OC(J,KT)=0.
 	do i=1,np 
	   if(wet(i)==0) OC(J,KT) = OC(J,KT) + & 
	   ALF(J,ANC(I))*CONJG(LONG_TABLE(MM(J),I))
 	enddo
	OC(J,KT)=OC(J,KT)/FLOAT(NP)
 enddo 
!$OMP END PARALLEL DO
!
300 Continue
!
! ----- Updated Feb. 2012: This constraint is removed
!   oc(:,nn+1) = oc(:,nn)
!
    do 500 k=0, nn+1      
!
	write( labchar, '(i2.2)' ) k
!	open(3,file='junk.dat',status='unknown') 
!	if(k<=9) write(3,'(a1,i1)') '0',k  
!	if(k> 9) write(3,'(i2)')        k
!	close(3)
!       open(3,file='junk.dat',status='unknown')
!	read(3,'(a2)') labchar
!	close(3)    
     
! --- Writing the result on "sh-of-XX.dat"
!
 open(1,file='shof-'//labchar//'.dat',status='unknown')  

 do j=1, jmax   
   	if(j==1) write(1,*) j, 1.- real(oc(j,k)), -aimag(oc(j,k))
   	if(j/=1) write(1,*) j,   - real(oc(j,k)), -aimag(oc(j,k))
 enddo
 close(1) 
!	   
 500 continue 
!
!
!
! --- Free up memory space
!
 DEALLOCATE( MM, WET, ANC )
 DEALLOCATE( LONP, LATP )
 DEALLOCATE( ALF )
 DEALLOCATE( LONG_TABLE, OC ) 
!
 END PROGRAM SH_OF_VAROC 
!
!
!
!
!
!
