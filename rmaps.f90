!
! This is program "RMAPS.F90" 
!
! Last modified GS 04-11-2008 "Intel port"
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
! *** Revised GS July 2010 - g95 - Double precision implementation 
! === Revised GS & FC May 21 2011 - DELTA parameter in ice history     
! === Parallelized DM December 2011
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
! 
! Produces maps of the rate of present-day sealevel change (dot S), of 
! vertical uplift (dot U), and of geoid change (dot N). Units are mm/yr
!
! Input files: 
! 	- shs.bin
!	- shu.bin
!	- shn.bin
!	- shz.bin 
!	- shice.dat 
! 	- sh.bin
! 	- px-table.dat
!       - anchor.tmp      [NEW]
!
! Output files:
! 	- sdotmap.dat 
!	- udotmap.dat
!	- ndotmap.dat
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!INCLUDE "harmonics.f90"
 PROGRAM RMAPS
 IMPLICIT NONE 
 INCLUDE "data.inc"
 CHARACTER*30 HEADER
 CHARACTER*11, PARAMETER :: UF='UNFORMATTED', F='FORMATTED', UN='UNKNOWN' 
 INTEGER :: NANCH
 INTEGER I, J, K, DOM, MJ
 INTEGER, ALLOCATABLE :: MM(:), DM(:), ANC(:)
 REAL*8, PARAMETER :: EPS = 0.01
 REAL*8, ALLOCATABLE :: RATE_S(:), RATE_U(:), RATE_N(:)
 REAL*8, ALLOCATABLE :: LON(:), LAT(:), ALF(:,:)
 REAL*8 LONX, LATX
 COMPLEX*16, ALLOCATABLE :: CS(:,:), CU(:,:), CN(:,:), LONG_TABLE(:,:)
 CHARACTER*13, PARAMETER :: FMT='(3(f10.6,1X))'
!
 REAL*8, PARAMETER :: DDELTA=DELTA 
 INTEGER :: IT1, IT2
 REAL*8  :: DT
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
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
! --- Allocating memory
! Updated February 2012 with "nn+1"
!
  ALLOCATE( MM(JMAX), DM(JMAX), ANC(NP) )
  ALLOCATE( RATE_S(NP), RATE_U(NP), RATE_N(NP) )
  ALLOCATE( LON(NP), LAT(NP) )
  ALLOCATE( ALF(JMAX,NANCH) )
  ALLOCATE( CS(JMAX,0:NN+1), CU(JMAX,0:NN+1), CN(JMAX,0:NN+1) )
  ALLOCATE( LONG_TABLE(0:LMAX,NP) )
!
!
! --- Pre-computing the degree 'm' corresponding to 'J'
!
!Write(*,*) "    - rmaps.f: Pre-computing the harmonic order"
 do j=1, jmax 
        mm(j)=mj(j)
	dm(j)=2-dom(j) 
 enddo	
!
!
!---- Reading the pixels table file... 
!
!Write(*,*) "    - rmaps.f: Reading data from the pixels table..."
 open(2,file='px-table.dat',status='unknown')
 do i=1, 4 
    read(2,'(a30)') HEADER 	
 enddo 
!Write(*,*) "    - rmaps.f: reading file <<px-table.dat>>"
 do i=1, np
    read(2,*) lon(i), lat(i), anc(i)
 enddo 
 close(2)  
!
!
!---- Reading the table of spherical harmonics from <<sh.bin>>
!
 Write(*,*) "    - Reading ALFs & TRIGs from <<sh.bin>>"
 open(3,file='sh.bin',status='unknown',form='unformatted')  
 	read(3)ALF
 	read(3)LONG_TABLE
 Close(3) 
 DO J=1, JMAX 
	ALF(J,:)=ALF(J,:)*DM(J)
 ENDDO
!
!
!
!---- Computing S-dot, U-dot, and N-dot at pixels 
!  
!Write(*,*) "    - rmaps.f: Reading the harmonics from binary files"
!
 open (101,file='shs.bin',status='unknown',form=uf) 
 open (102,file='shu.bin',status='unknown',form=uf) 
 open (103,file='shn.bin',status='unknown',form=uf) 
 read (101) CS
 read (102) CU
 read (103) CN
 close(101)
 close(102)
 close(103)  
!
!
!
 IF( IDER==1 ) THEN
    IT1 = NN-1
    IT2 = NN
    DT  = DDELTA
 ENDIF
 IF( IDER==2 ) THEN
    IT1 = NN-1
    IT2 = NN+1
    DT  = 2.D0*DDELTA
 ENDIF      
!
!
!
 Write(*,*) "    - Computing S, U and N-dot at", np, " pixels "
!
!
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J) &
!$OMP    SHARED(RATE_S,RATE_U,RATE_N,ALF,ANC,CS,CU,CN,LONG_TABLE,MM,IT1,IT2,DT) &
!$OMP    SCHEDULE(GUIDED)
 do i=1, np 
!
       rate_s(i)=0. 
       rate_u(i)=0. 
       rate_n(i)=0. 
!
       do j=1, jmax 
!
           rate_s(i) = rate_s(i) + & 
	          ALF(j,anc(i))*real(((cs(j,it2)-cs(j,it1))/DT)*LONG_TABLE(MM(J),I))  
! 
           rate_u(i) = rate_u(i) + & 
	          ALF(j,anc(i))*real(((cu(j,it2)-cu(j,it1))/DT)*LONG_TABLE(MM(J),I)) 
!
           rate_n(i) = rate_n(i) + & 
	          ALF(j,anc(i))*real(((cn(j,it2)-cn(j,it1))/DT)*LONG_TABLE(MM(J),I)) 
!
       enddo
 end do
!$OMP END PARALLEL DO      
!
!
!
!
!---- Writing S-dot, U-dot, and N-dot at pixels 
!  
 Write(*,*) "    - Writing S, U and N-dot on [s/u/n]dotmap.dat"
 open(101,file='sdotmap.dat',status='unknown') 
 open(102,file='udotmap.dat',status='unknown') 
 open(103,file='ndotmap.dat',status='unknown') 
!
!
 do i=1,np
!
! --- South pole
       If     (lat(i)==-90.) then  
       		latx=lat(i)+eps
       		do k=1, 4 
			lonx=(k-1)*90.   		
            		write(101,FMT) lonx, latx, rate_s(i) 
            		write(102,FMT) lonx, latx, rate_u(i) 
            		write(103,FMT) lonx, latx, rate_n(i) 
       		enddo
!
! --- North pole
       elseif (lat(i)==+90.) then   ! "North pole"
       		latx=lat(i)-eps
       		do k=1, 4 
			lonx=(k-1)*90.   			
            		write(101,FMT) lonx, latx, rate_s(i) 
            		write(102,FMT) lonx, latx, rate_u(i) 
            		write(103,FMT) lonx, latx, rate_n(i) 
       		enddo
!
! --- Elsewhere 
       elseif (lat(i)/=-90.and.lat(i)/=90.) then 
       		latx=lat(i) 
		lonx=lon(i) 
            		write(101,FMT) lonx, latx, rate_s(i) 
            		write(102,FMT) lonx, latx, rate_u(i) 
            		write(103,FMT) lonx, latx, rate_n(i) 
       Endif	       
!
 end do
!
   close(100)
   close(101)
   close(102)  
!
!
!
! --- Deallocate memory space
!
 DEALLOCATE( RATE_S, RATE_U, RATE_N )
 DEALLOCATE( MM, ANC, DM )
 DEALLOCATE( LON, LAT, ALF )
 DEALLOCATE( CS, CU, CN, LONG_TABLE )
   !
!
 END PROGRAM RMAPS
!
!
!
