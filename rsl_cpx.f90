!
! This is program "RSL_CPX.F90"  
!
! Last modified GS & FC 06-11-2008 "Varying coastlines"
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
! *** Revised GS December 2010 - & cleanup & Double precision implementation
! Feb 2012: Implementation of the numerical derivative "in the future" 
! *** Revised DM July 2012 - OpenMP parallelization
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
! ----------------------------------------------------------------
!  Computes the RSL curves at a number of <<virtual RSL sites>> 
! ----------------------------------------------------------------
! 
! Input files:
!	- pxa.dat
! 	- shs.bin
!       - px-topo.dat 
!       - ice-mask-XX.dat 
!       - anchor.tmp [NEW]
!
!
! Output files:
!       - px-table-XX.dat 
!	- oc-XX.dat 
!       - ptopo-XX.dat 
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!INCLUDE "harmonics.f90"
 PROGRAM SL
 IMPLICIT NONE
 INCLUDE "data.inc"
 INTEGER :: NANCH
 CHARACTER*20, PARAMETER :: PXFILE='pxa.dat'
 CHARACTER*20, PARAMETER :: FILEOUT1 = 'rslpx.dat'
 CHARACTER*30 HEADER
 CHARACTER*20 DATE, TIMC                
 CHARACTER*22 FILENAME, FILEPX, FILEMASK, FILEOC, FILEPTOPO 
 CHARACTER*42 JUNK
 CHARACTER*2 LABCHAR 
 CHARACTER*11, PARAMETER :: UF='UNFORMATTED', F='FORMATTED', UN='UNKNOWN' 
 CHARACTER*13, PARAMETER :: FMT='(3(f10.6,1X))'
 REAL*8, ALLOCATABLE :: LON(:), LAT(:), TOPO(:), ALF(:,:)
 REAL*8, ALLOCATABLE :: SLC(:,:), RSL(:,:)
 REAL*8 X, Y
 REAL*8, ALLOCATABLE :: ICETHICK(:,:), ICET(:)
 REAL*8 NEWTOPO, PALEOTOPO 
 COMPLEX*16, ALLOCATABLE :: CS(:,:), LONG_TABLE(:,:)
 INTEGER I, J, K, JJ, MJ, DOM, IJUNK, ICEMASK
 INTEGER, ALLOCATABLE :: ANC(:), MM(:), DM(:)  
!! REAL*8 LON(NP), LAT(NP), TOPO(NP), ALF(JMAX,NANCH)
!! REAL*8 SLC(NP,0:NN), RSL(0:NN,NP), X, Y
!! REAL*8 ICETHICK(NP,0:NN+1), NEWTOPO, PALEOTOPO 
!! COMPLEX*16 CS(JMAX,0:NN), LONG_TABLE(0:LMAX,NP)
!! INTEGER I, J, K, JJ, MJ, DOM, IJUNK, ANC(NP), MM(JMAX), DM(JMAX), ICEMASK  
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
! --- Allocate memory space   ! Updated variables as of February 2012 
!
  allocate( lon(np), lat(np), topo(np) )
  allocate( alf(jmax,nanch) )
  allocate( slc(np,0:NN+1), rsl(-1:NN,np) )
  allocate( icethick(np,0:nn+1), icet(np) )
  allocate( cs(jmax,0:nn+1) )
  allocate( long_table(0:lmax,np) )
  allocate( anc(np), mm(jmax), dm(jmax) )
!
! --- Reading the present-day pixelized TOPO 
!
 open(3,file=pxtopo_file,status='unknown') 
 do i=1, 4 
    read(3,'(a30)') HEADER 	
 enddo 
 do i=1, np 
      read(3,*) x, y, ijunk, ijunk, ijunk, topo(i)     	        
 enddo
 close(3)
!
!
! --- Pre-computing the degree 'm' corresponding to 'J'
 do j=1, jmax 
        mm(j)=mj(j)
	dm(j)=2-dom(j) 
 enddo	
!
!
!---- Reading the pixels table file... 
 open(2,file=pxfile,status='unknown')
 do i=1, 4 
    read(2,'(a30)') HEADER 	
 enddo 
 do i=1, np
    read(2,*) lon(i), lat(i), anc(i)
 enddo 
 close(2)  
!
!
!---- Reading the table of spherical harmonics from <<sh.bin>>
 open(3,file='sh.bin',status='unknown',form='unformatted')  
 	read(3)ALF
 	read(3)LONG_TABLE
 Close(3) 
 DO J=1, JMAX 
	ALF(J,:)=ALF(J,:)*DM(J)
 ENDDO
!
!
!---- Reading the SH coefficients for "S" 
 open (101,file='shs.bin',status='unknown',form='UNFORMATTED') 
 read (101) CS
 close(101)  
!
!
!---- Computing "S" 
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP    SHARED(SLC,ALF,CS,LONG_TABLE,MM,ANC) &
!$OMP    PRIVATE(I,J) SCHEDULE(GUIDED)
  do i=1, np
	slc(i,:) = 0. 
	do j=1, jmax
	slc(i,:) = slc(i,:) + ALF(j,anc(i))*real((cs(j,:))*LONG_TABLE(MM(J),I))
	enddo
  enddo
!$OMP END PARALLEL DO
!
	Write(*,*) "    - Computing RSL change at all pixels at time"
!
! Here "k" is time BP
! Updated Feb. 2012 with a step in the future
!	
	do k=-1, nn
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(RSL,SLC,K) PRIVATE(I) SCHEDULE(GUIDED)	
 			do i=1, np
		        rsl(k,i) = -(slc(i,nn)-slc(i,nn-k))
			enddo
!$OMP END PARALLEL DO			
	enddo 
        Write(*,*) "    - Paleotopography by backward steps"	
!
!
	do 300 k=-1, nn 
!
! Here "k" is time BP
!
! jj=0 gives the OF before the beginning of melting 
! jj=N is for the last step.   
        jj=nn-k     ! with k above, this provides jj=nn+1 (it was 'nn')
!
	!open(3,file='junk.dat',status='unknown') 
	!if(jj<=9) write(3,'(a1,i1)') '0',jj  
	!if(jj> 9) write(3,'(i2)')        jj
 	!close(3)
	!open(3,file='junk.dat',status='unknown'); read(3,'(a2)') labchar; close(3) 
    write(labchar,'(i2.2)') jj
!          
	If(k==0.or.k==nn) write(*,'(a12,a2,a4,i2)') '     - Step ', labchar, ' of ',NN
!
	filepx    = 'px-table-'//labchar//'.dat'
	fileoc    = 'oc-'//labchar//'.dat'
	filemask  = 'imask-'//labchar//'.dat'
	fileptopo = 'ptopo-'//labchar//'.dat'
!
	Open(14,file=filepx,    status='unknown')
	Open(15,file=fileoc,    status='unknown')
	Open(16,file=filemask,  status='unknown')
	Open(17,file=fileptopo, status='unknown')  
!
!
	write(14,*) "-------------------------------------------"
	write(14,*) " Time-varying OCEAN function for time ", jj  
	write(14,*) "       **** GS & FC - June 2008 ****"
	write(14,*) "-------------------------------------------"
!
        do 200 i=1, NP 
! 
! Reading the ice thickness corresponding to time "jj=nn-k"
!
		read(16,*) x, y, icemask, icet(i) 
!
        	newtopo   = topo(i) - rsl(k,i)
!	
		paleotopo = newtopo + icet(i) 
!
!New topo is NEGATIVE and there is NO ice ===> OCEAN   
      		if(newtopo<=0.and.icemask==0) then 
              		write(14,*) lon(i), lat(i), anc(i), i, '1' 
              		write(15,*) lon(i), lat(i) 
      		endif 

!New topo is NEGATIVE and there is ice ===> CONTINENT		      
      		if(newtopo<=0.and.icemask==1) then	      
              		write(14,*) lon(i), lat(i), anc(i), i, '0' 
      		endif   

!New topo is POSITIVE ===> CONTINENT		      
      		if(newtopo>0) then  
              		write(14,*) lon(i), lat(i), anc(i), i, '0'      
      		endif 
!
     	        Write (17,*) lon(i), lat(i), newtopo, paleotopo
!
!
! OLD PIECE OF CODE...OLD PIECE OF CODE...OLD PIECE OF CODE...OLD PIECE OF CODE...
!
! New topo is NEGATIVE and there is NO ice ===> OCEAN 	
!	if(paleotopo<=0) then 
!		write(14,*) lon(i), lat(i), anc(i), i, '1' 
!		write(15,*) lon(i), lat(i) 
!	endif 
!
! New topo is NEGATIVE and there is ice ===> CONTINENT			
!	if(newtopo<=0.and.icemask==1) then        	
!		write(14,*) lon(i), lat(i), anc(i), i, '0' 
!	endif 	
!
! New topo is POSITIVE ===> CONTINENT			
!	if(paleotopo>0) then  
!		write(14,*) lon(i), lat(i), anc(i), i, '0' 	
!	endif 
!
!	Write (17,*) lon(i), lat(i), newtopo, paleotopo
!
! OLD PIECE OF CODE...OLD PIECE OF CODE...OLD PIECE OF CODE...OLD PIECE OF CODE...
!
!	
 200    continue
!
        close(14)
	close(15)  
	close(16) 
	close(17) 
!
 300    CONTINUE
!
!
!
! --- Free memory space
!
  deallocate( lon, lat, topo )
  deallocate( alf )
  deallocate( slc, rsl )
  deallocate( icethick )
  deallocate( cs )
  deallocate( long_table )
  deallocate( anc, mm, dm )
  deallocate( icet )
!
!
!
 END PROGRAM SL
!
!
