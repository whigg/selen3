!
! Program "GMAPS_ELA_MS_DOT.F90"
! 
! This is largely BASED on program "GMAPS.F90" 
! Created by GS on February 22, 2012 
! Based on previous "GMAPS_ELA_MS_DOT.F90"
! Modified in the subsequent weeks a number of times 
! GS March 08 2012 - REAL -> DBLE and expressions rearranged to preserve precision in
!                    the computation of time derivatives. 
! Including day March 11, 2012 for the U maps
! Option GIA on April 19 (Urbino Hospital...) ---> now working 
! Modified by GS on June 25 for the multi-stepping for the ice2sea computations 
! Modified by DM on Nov 09 2012 - dynamic memory allocation & external NANCH
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
! For a multi-step elastic rebound, this program produces maps of: 
! RATE of Sea level change (S-dot)  
! RATE of Sea surface (N-dot)  

!
! Input files: 
! 	- shs.bin, shn.bin. shu.bin
! 	- sh.bin
! 	- px-table.dat
!
! Output files:
! 	- sdotmap.dat 
! 	- ndotmap.dat 
! 	- udotmap.dat 
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!INCLUDE "harmonics.f90"
 PROGRAM G 
 IMPLICIT NONE 
 INCLUDE "data.inc"
!
 CHARACTER*30 HEADER
 CHARACTER*13, PARAMETER :: FMT ='(3(f10.6,1X))'
 CHARACTER*13, PARAMETER :: FMT2='(3(f14.6,1X))'

 CHARACTER*11, PARAMETER :: UF='UNFORMATTED', F='FORMATTED', UN='UNKNOWN' 
 REAL*8, PARAMETER :: EPS = 0.01
 CHARACTER*3 LABCHAR 
 CHARACTER*22 INPUT_FILENAME, OUTPUT_FILENAME 
 INTEGER :: NANCH
 INTEGER I, J, K, KT, DOM, MJ
 INTEGER, ALLOCATABLE :: MM(:), ANC(:), DM(:)  
 REAL*8 LONX, LATX
 REAL*8, ALLOCATABLE :: LON(:), LAT(:), ALF(:,:) 
 REAL*8, PARAMETER :: DDELTA=DELTA 
 REAL*8 DOTS, DOTN, DOTU 
!
! Updated February 2012 with "nn+1"
!
 COMPLEX*16 LONG_TABLE(0:LMAX,NP)
!
 COMPLEX*16 CS(JMAX,0:NN+1), & 
            CN(JMAX,0:NN+1), & 
            CU(JMAX,0:NN+1) 	    
!
!
! PATH to the GIA corrections and misc variables 
!
 CHARACTER*30, PARAMETER :: GIA_PATH="./DATA/GIA-corrections/" 
 CHARACTER*40, PARAMETER :: file_dots_gia_corr = trim(adjustl(GIA_PATH))//"sdotmap.dat"
 CHARACTER*40, PARAMETER :: file_dotn_gia_corr = trim(adjustl(GIA_PATH))//"ndotmap.dat"
 CHARACTER*40, PARAMETER :: file_dotu_gia_corr = trim(adjustl(GIA_PATH))//"udotmap.dat"
 REAL*8 XLON, XLAT
 REAL*8 SDOT_GIA_CORR, UDOT_GIA_CORR, NDOT_GIA_CORR  
 REAL*8 AJUNK, BJUNK, CJUNK  
 REAL*8, PARAMETER :: ZERO=0D0, MNONA=-90D0, PNONA=+90D0
!
!
!	   
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
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
!
  allocate( lon(np), lat(np), anc(np) )
  allocate( mm(jmax), dm(jmax) )
  allocate( alf(jmax,nanch) )!
!
!
! --- Pre-computing the degree 'm' corresponding to 'J'
!
 do j=1, jmax 
        mm(j)=mj(j)
	dm(j)=2-dom(j) 
 enddo	
!
!
!---- Reading the pixels table file... 
!
 open(2,file='px-table.dat',status='unknown')
! 
 if(GLOB_ELA_S_DOT==1)  Write(*,*) "    - Computing S-dot at ", np, " pixels "
 if(GLOB_ELA_N_DOT==1)  Write(*,*) "    - Computing N-dot at ", np, " pixels "
 if(GLOB_ELA_U_DOT==1)  Write(*,*) "    - Computing U-dot at ", np, " pixels "
!
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
!
 Write(*,*) "    - Reading ALFs & TRIGs from sh.bin"
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
!---- Reading the tables of spherical harmonics  
!
 if(GLOB_ELA_S_DOT==1) then 
 	    open (1,file='shs.bin', status='unknown',form=uf) 
 	    read (1) CS ; close(1) 
 endif
  if(GLOB_ELA_N_DOT==1) then 
 	    open (1,file='shn.bin', status='unknown',form=uf) 
 	    read (1) CN ; close(1) 
 endif
  if(GLOB_ELA_U_DOT==1) then 
 	    open (1,file='shu.bin', status='unknown',form=uf) 
 	    read (1) CU ; close(1) 
 endif
!
!
!

!
!---- Syncronizing the GIA corrections files 
!
      IF(GIA_CORR==1) THEN  

      DO K=1, 3 
      
      IF    (K==1) then 
      		OPEN(600,file=file_dots_gia_corr,status='unknown')
      		OPEN(601,file="dots.tmp",status='unknown')      
      ELSEIF(K==2) then       
      		OPEN(600,file=file_dotu_gia_corr,status='unknown')
      		OPEN(601,file="dotn.tmp",status='unknown')      
      ELSEIF(K==3) then       
      		OPEN(600,file=file_dotn_gia_corr,status='unknown')
      		OPEN(601,file="dotu.tmp",status='unknown')      
      ENDIF
!
!
      DO I=1, NP+6 
!      
          READ (600,FMT)  AJUNK, BJUNK, CJUNK   
!		   
          IF(I==1)                  THEN 
                   WRITE(601,FMT)               ZERO, MNONA, CJUNK 
!		   
		   ELSEIF(I==NP+6)  THEN 
		               WRITE(601,FMT)   ZERO, PNONA, CJUNK 
!		   
		   ELSEIF(I.GE.5.AND.I.LE.NP+3)   THEN 
		               WRITE(601,FMT)  AJUNK, BJUNK, CJUNK  		                 
		   ENDIF		   
!      
      ENDDO
!
	CLOSE(600)
	CLOSE(601) 
!
      ENDDO
!      
      ENDIF
!
!
!
!
! ***********************************************************************
! 
! ***********************************************************************
!
!---  Loop on the time steps ----
!
! ***********************************************************************
! 
! ***********************************************************************
!
!
 DO 100 KT=1,NN-1    ! Note: from "1" to "N-1" to avoid border effects   
!
!open(3,file='junk.dat',status='unknown') 
!if(kt<=9) write(3,'(a1,i1)') '0',kt  
!if(kt> 9) write(3,'(i2)')        kt ; close(3)
!open(3,file='junk.dat',status='unknown'); read(3,'(a2)') labchar; close(3)         
!
! Modified by GS on June 25 for the multi-stepping for the ice2sea computations 
!
open(3,file='junk.dat',status='unknown')
if(kt<=9)              write(3,'(a2,i1)') '00',kt  
if(kt>=10.and.kt <=99) write(3,'(a1,i2)')  '0',kt 
if(kt>=100)            write(3,   '(i3)')      kt ; close(3)
open(3,file='junk.dat',status='unknown'); read(3,'(a3)') labchar; close(3) 
!
		if(kt==0.or.kt==nn+1) write(*,*) '    - Step ', labchar, ' of ',NN+1 
!
                IF(GLOB_ELA_S_DOT==1) then 
		        input_filename='sdotmap-'//labchar//'.dat'
		        open(101,file=input_filename,status='unknown')
		Endif
                IF(GLOB_ELA_N_DOT==1) then 
		        input_filename='ndotmap-'//labchar//'.dat'
		        open(102,file=input_filename,status='unknown')
		Endif
                IF(GLOB_ELA_U_DOT==1) then 
		        input_filename='udotmap-'//labchar//'.dat'
		        open(103,file=input_filename,status='unknown')
		Endif
!
!
      IF(GIA_CORR==1) THEN  
      		OPEN(401,file="dots.tmp",status='unknown')
      		OPEN(402,file="dotn.tmp",status='unknown')
      		OPEN(403,file="dotu.tmp",status='unknown')
      ENDIF
!
!---  Loop on pixels ----
!
  DO 1 I=1, NP 
!
     IF(GLOB_ELA_S_DOT==1) DOTS   = 0D0 
     IF(GLOB_ELA_N_DOT==1) DOTN   = 0D0
     IF(GLOB_ELA_U_DOT==1) DOTU   = 0D0
!
     SDOT_GIA_CORR=0D0
     NDOT_GIA_CORR=0D0
     UDOT_GIA_CORR=0D0
!
      IF    (GIA_CORR==1) THEN  
      		READ(401,FMT) XLON, XLAT, SDOT_GIA_CORR 
      		READ(402,FMT) XLON, XLAT, NDOT_GIA_CORR 
      		READ(403,FMT) XLON, XLAT, UDOT_GIA_CORR 
      ENDIF
!
!
!--- Loop on the SHs ----
!
     DO 2 J=1, JMAX 
!
       IF(GLOB_ELA_S_DOT==1) DOTS = DOTS + & 
                             ALF(J,ANC(I))*DBLE(((CS(J,KT+1)-CS(J,KT-1)))*LONG_TABLE(MM(J),I))
!                                           
       IF(GLOB_ELA_N_DOT==1) DOTN = DOTN + & 
                             ALF(J,ANC(I))*DBLE(((CN(J,KT+1)-CN(J,KT-1)))*LONG_TABLE(MM(J),I))
!
       IF(GLOB_ELA_U_DOT==1) DOTU = DOTU + & 
                             ALF(J,ANC(I))*DBLE(((CU(J,KT+1)-CU(J,KT-1)))*LONG_TABLE(MM(J),I))

!
2    CONTINUE
!
! --- South pole
       If     (lat(i)==-90.) then  
       		latx=lat(i)+eps
       		do k=1, 4 
			lonx=(k-1)*90.   		
            		if(GLOB_ELA_S_DOT==1) write(101,FMT) lonx, latx, dots/2D0/DDELTA + SDOT_GIA_CORR
            		if(GLOB_ELA_N_DOT==1) write(102,FMT) lonx, latx, dotn/2D0/DDELTA + NDOT_GIA_CORR
            		if(GLOB_ELA_U_DOT==1) write(103,FMT) lonx, latx, dotu/2D0/DDELTA + UDOT_GIA_CORR
       		enddo
!
! --- North pole
       elseif (lat(i)==+90.) then   
       		latx=lat(i)-eps
       		do k=1, 4 
			lonx=(k-1)*90.   		
            		if(GLOB_ELA_S_DOT==1) write(101,FMT) lonx, latx, dots/2D0/DDELTA + SDOT_GIA_CORR
            		if(GLOB_ELA_N_DOT==1) write(102,FMT) lonx, latx, dotn/2D0/DDELTA + NDOT_GIA_CORR
            		if(GLOB_ELA_U_DOT==1) write(103,FMT) lonx, latx, dotu/2D0/DDELTA + UDOT_GIA_CORR
       		enddo
!
! --- Elsewhere 
       elseif (lat(i)/=-90.and.lat(i)/=90.) then 
       		latx=lat(i) 
		lonx=lon(i) 
            	        if(GLOB_ELA_S_DOT==1) write(101,FMT) lonx, latx, dots/2D0/DDELTA + SDOT_GIA_CORR
            		if(GLOB_ELA_N_DOT==1) write(102,FMT) lonx, latx, dotn/2D0/DDELTA + NDOT_GIA_CORR
            		if(GLOB_ELA_U_DOT==1) write(103,FMT) lonx, latx, dotu/2D0/DDELTA + UDOT_GIA_CORR
!
       Endif	       
!
1 CONTINUE
!
  CLOSE(101) 
  CLOSE(102)
  CLOSE(103) 
!
  CLOSE(401)
  CLOSE(402)
  CLOSE(403)
! 
100 CONTINUE
!
! --- Release memory
!
  deallocate( lon, lat, anc )
  deallocate( mm, dm )
  deallocate( alf )
!
 END PROGRAM G 
!
!
!
!
!
!
