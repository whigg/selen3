!
! This is program "SHICE.F90" 
!
! Created by GS July 8 2008  -Intel version 2.6- 
! Updated GS July 24 for implementation of the "disk load"
! Modified August 2008 for implementation of ANU ice on v. 2.7
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
! Revised GS May 2010 for the ice breaker mission 
! Reviewed GS on July 2010 for the implementation of the Greenland ice sheet 
! On July 19 found a BAD MISTAKE -  Now Fixed
! Reviewed GS on July 2010 for the implementation of the small glaciers 
! ***  Revised GS July 2010 - g95 - Double precision implementation
! ***  Revised GS & GR Nov 2010 - Implementation of MED1
! ***  Revised GS & GG May 2011 - Implementation of the Antarctic model "ANTA"
! ***  Revised GS & GG June 2011 - Implementation of the "TRANS" ice model 
! ***  Revised GS July 2011 - Implementation of the "RECT" ice model 
! ---  Added some lines for a backup of ice models (ICE-3G, ... ) by GS
! Feb 2012: Implementation of the numerical derivative "in the future"
! Feb 2012: GS is adding a step "nn+2" into the "h" array 
! === February 2012: Introduced the "GRD" ice model type by GS
! +++ July 1, 2012: A scaling factor is introduced for the elastic rebound files
! Revised DM Nov 2012 for performance improvement
! Jul 1, 2013; The "FLO" ice for coupling with her models 
! April 11, 2014. The BENOIT (BENO) files for the collaboration with LEGOS 
! 
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
!  This program computes the SH coefficients of the ice models in the library 
!  ./ICE-MODELS, based on pre-computed shape factors. This routine collects 
!  in a single script various program units that were written "ad-hoc" for 
!  individual ice models until version 2.5 of SELEN (i.e. SH1.F90, etc..)
! ---------------------------------------------------------------------------
!
!
!
!      INCLUDE "harmonics.f90"
       PROGRAM SHICE 
       IMPLICIT NONE        
       INCLUDE "data.inc"
!
!#---- General declarations 
!      COMPLEX*16 TTTT(JMAX,0:NN+1), IIII(JMAX,0:NN) 
!  
!      COMPLEX*16 TTTT(JMAX,0:NN+2), IIII(JMAX,0:NN+1)   ! New Feb 2012
       COMPLEX*16, ALLOCATABLE :: TTTT(:,:), IIII(:,:)   
!
       INTEGER, PARAMETER :: NH1=20, NH3=20, NH5=24 
       INTEGER, PARAMETER :: NHA=4,  NHI=18, NHD=20 
       INTEGER, PARAMETER :: NHU=26, NHC=20, NHG=10  
       INTEGER, PARAMETER :: NHL=10, NHM=20, NHR=10 
       INTEGER, PARAMETER :: NHT=20, NHY=20, NHZ=10 
       INTEGER, PARAMETER :: NHF=4,  NHB=4      
!            
       INTEGER I, J, K, L, N, KK, IJ, JJ, NH, I1, I2, ISLOT, IDEG, IEL, ICR(0:NN) 
       REAL*8 AJ, AJ1, AJ2, AJ3, AJ4  
       REAL*8, ALLOCATABLE :: CR(:), CCRR(:) ! CR(0:NN), CCRR(0:NN)
!
!      REAL*8 H(NEL,0:NN+1)
!      REAL*8 H(NEL,0:NN+2) ! New Feb 2012             
       REAL*8, ALLOCATABLE :: H(:,:)
!
       CHARACTER*2 CJ
!
!
       INTEGER, PARAMETER :: SLOT_SIZE=1000 ! Consistent with ICE_BREAKER.F90
       INTEGER II, KS, JDEG, NELL, NSLOTS, LO(151), HI(151) 
       CHARACTER*2 LABCHAR 
       CHARACTER*20 broken_ice_file(151)
!       COMPLEX*16 ARRAY(JMAX,1:SLOT_SIZE)
       COMPLEX*16, ALLOCATABLE :: ARRAY(:,:)


!
!
!#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
! - Output files:  
        CHARACTER*20, PARAMETER :: FILEOUT='shice.dat'
!       shtXX.dat, SH files at time XX (XX=00, 01, ..,N, N+1) 
!#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!
!
!#---- Allocate memory
!     
       allocate( array(jmax,1:slot_size) )
       allocate( TTTT(JMAX,0:NN+2), IIII(JMAX,0:NN+1) )
       ALLOCATE( H(NEL,0:NN+2), CR(0:NN), CCRR(0:NN) )     
!
!
!
!#---- Reading the header lines of the ice table
!	    	
            OPEN(10,FILE=ICE_MODEL,STATUS='UNKNOWN') 
!
	    	  if(ice_model(1:4)=='ice1')nh=nh1 
	    	  if(ice_model(1:4)=='ice3')nh=nh3
	    	  if(ice_model(1:4)=='imed')nh=nhm
	    	  if(ice_model(1:4)=='ice5')nh=nh5
	    	  if(ice_model(1:4)=='alps')nh=nha
	    	  if(ice_model(1:4)=='flor')nh=nhf
	    	  if(ice_model(1:4)=='ij05')nh=nhi
	    	  if(ice_model(1:4)=='disk')nh=nhd
		  if(ice_model(1:4)=='anu0')nh=nhu
	    	  if(ice_model(1:4)=='icap')nh=nhc
		  if(ice_model(1:4)=='gree')nh=nhg
		  if(ice_model(1:4)=='anta')nh=nhr		 
		  if(ice_model(1:4)=='glac')nh=nhl
		  if(ice_model(1:4)=='tran')nh=nht		  
		  if(ice_model(1:4)=='rett')nh=nhy
                  if(ice_model(1:4)=='beno')nh=nhb					
                  if(ice_model(1:3)=='grd') nh=nhz		
!		  
            	do j=1, nh 
	        	read(10,'(a20)') cj
!write(*,'(a20)') cj
!read(*,*) 
	    	enddo	     
!
		IF     (ICE_MODEL(1:4)=='ice3'.or.ICE_MODEL(1:4)=='ice5') then 
		Open(44,file=ice_model//".bck",status='unknown')
	 	do ii=1, nh 
	 		write(44,*)"A backup of ice model ", ice_model
	 	enddo
		ENDIF 

!	    
!#---- Reading the time-dependent thickness of the ice elements 
!	    	
  DO 10 I=1, NEL   
!
!
!
!#---- "RECT" ---- -----
!
 IF(ICE_MODEL(1:4)=='rett') THEN 
  Read (10,*) ij, aj, aj, aj, h(i,0) 
  Read (10,*) (h(i,k),k=1 ,10)
  Read (10,*) (h(i,k),k=11,20)
  Read (10,*) (h(i,k),k=21,30)
  Read (10,*) (h(i,k),k=31,40)
  Read (10,*) (h(i,k),k=41,50)
  h(i,nn+1)=h(i,nn)
  h(i,nn+2)=h(i,nn)    ! New Feb 2012
!
!
!#---- "TRANS" ---- -----
!
 ELSEIF(ICE_MODEL(1:4)=='tran') THEN 
  Read (10,*) ij, aj, aj, aj, h(i,0) 
  Read (10,*) (h(i,k),k=1 ,10)
  Read (10,*) (h(i,k),k=11,20)
  Read (10,*) (h(i,k),k=21,30)
  Read (10,*) (h(i,k),k=31,40)
  Read (10,*) (h(i,k),k=41,50)
  h(i,nn+1)=h(i,nn)
  h(i,nn+2)=h(i,nn)    ! New Feb 2012
!
!
!#---- "GRD" ---- -----          ! New Feb 2012
!
 ELSEIF(ICE_MODEL(1:3)=='grd') THEN 
  Read (10,*) ij, aj, aj, aj, (h(i,k),k=0,nn+1)
  h(i,nn+2)=h(i,nn+1)  
!
 ELSEIF(ICE_MODEL(1:4)=='beno') THEN 
  Read (10,*) ij, aj, aj, aj, (h(i,k),k=0,nn)
  h(i,nn+1)=h(i,nn)   
  h(i,nn+2)=h(i,nn+1)   

!
!#---- ------ ---- -----
!	    
!#---- "ICE1" ----             	    
	    ELSEIF(ICE_MODEL(1:4)=='ice1') THEN 
        	   Read (10,*) ij, ij,                (icr(k),k=nn,0,-1) 
		   cr(:)=icr(:)			
!
!#---- "ICE3" ----   
	    ELSEIF(ICE_MODEL(1:4)=='ice3') THEN 
       		   Read (10,111) ij, aj1, aj2, aj3, aj4,  (cr(k),k=nn,11,-1)
       		   Read (10,112)                          (cr(k),k=10, 0,-1)
!
!#---- "IMED" ----   
	    ELSEIF(ICE_MODEL(1:4)=='imed') THEN 
       		   Read (10,111) ij, aj, aj, aj, aj,  (cr(k),k=nn,11,-1)
       		   Read (10,112)                      (cr(k),k=10, 0,-1)
!
!#---- "ANU05" ----   
	    ELSEIF(ICE_MODEL(1:4)=='anu0') THEN 
       		   Read (10,115) ij, aj, aj, aj,      (h(i,k),k=0,nn+1)
!
!#---- "DISK" ----   
	    ELSEIF(ICE_MODEL(1:4)=='disk') THEN 
       		   Read (10,*)   ij, aj, aj, aj,      (cr(k),k=nn,11,-1)
       		   Read (10,*)                        (cr(k),k=10, 0,-1)
!
!#---- "ICAP" ----   
	    ELSEIF(ICE_MODEL(1:4)=='icap') THEN 
       		   Read (10,*)   ij, aj, aj, aj,      (cr(k),k=nn,11,-1)
       		   Read (10,*)                        (cr(k),k=10, 0,-1)
!
!#---- "ICE5" ----
	    ELSEIF(ICE_MODEL(1:4)=='ice5') THEN 
!	    
	        If (ice_model/='ice5g26.dat')then 
       		    Read (10,113) ij, aj1, ij, aj2,     (cr(k),k=nn,0,-1) 			
		else
		    Read (10,114) ij, aj1, ij, aj2,     (cr(k),k=nn,0,-1) 
		Endif
!		        		
!#---- "ANTARCTICA" ----
	    ELSEIF(ICE_MODEL(1:4)=='anta') THEN 
       		   Read (10,*) ij, aj, aj, aj,        (cr(k),k=nn,0,-1) 
!do k=nn, 0, -1
!cr(k)=cr(k)*SCALING_ICE
!enddo
!		   
!		        		
!#---- "GREENLAND" ----
	    ELSEIF(ICE_MODEL(1:4)=='gree') THEN 
       		   Read (10,*) ij, aj, aj, aj,        (cr(k),k=nn,0,-1) 
!do k=nn, 0, -1
!cr(k)=cr(k)*SCALING_ICE
!enddo
!		        		
!#---- "Small Glaciers..." ----
	    ELSEIF(ICE_MODEL(1:4)=='glac') THEN 
       		   Read (10,*) ij, aj, aj, aj,        (cr(k),k=nn,0,-1) 
!do k=nn, 0, -1
!cr(k)=cr(k)*SCALING_ICE
!enddo
!
!		        		
!#---- "ALPS" ----
	    ELSEIF(ICE_MODEL(1:4)=='alps') THEN 
       		   Read (10,*) ij, aj, aj, aj,        (cr(k),k=nn,0,-1) 

!		        		
!#---- "FLOR" ----   !!!! NEW HERE 
	    ELSEIF(ICE_MODEL(1:4)=='flor') THEN 
!           write(*,*) ICE_MODEL
            Read (10,*) ij, aj, aj, aj,        (cr(k),k=nn,0,-1) 
!write (*,*) ij, aj, aj, aj,        (cr(k),k=nn,0,-1) 
!read(*,*) 

!
!        		
!#---- "GRD" ----
!ELSEIF(ICE_MODEL(1:3)=='grd')  THEN 
!Read (10,*) ij, aj, aj, aj,        (cr(k),k=nn,0,-1) 

!
!#---- "IJ05" ----
	    ELSEIF(ICE_MODEL(1:4)=='ij05') THEN 
       		   Read (10,*)  ij, aj, aj, aj, aj,   (cr(k),k=nn,9,-1) 
        	   Read (10,*) 	                      (cr(k),k= 8,0,-1)  ! 2nd row
!
!#---- "End of available ice models..." 
            ENDIF
!
!
!#---- Conversion to our default for ice thickness
!
         If(ICE_MODEL(1:4)/='anu0'.and.& 
	    ICE_MODEL(1:4)/='tran'.and.&
	    ICE_MODEL(1:4)/='rett'.and.& 
	    ICE_MODEL(1:4)/='beno'.and.&     	    
	    ICE_MODEL(1:3)/='grd')      then 
         	h(i,0)=cr(nn)              
       		do k=1,nn                   
                     h(i,k)=cr(nn-k)   
       		enddo 
                h(i,nn+1) = h(i,nn) 
		h(i,nn+2) = h(i,nn)    ! New Feb 2012
!write(*,*) "PIPPO", (h(i,k),k=0,nn+2) OK OK OK OK OK OK OK 
	Endif
!
!#---- Backup of the ice file  --------------------------------------------------- new new new 
!
	IF    (ICE_MODEL(1:4)=='ice3') THEN 	 
!
         do k=0,nn                   
                     ccrr(k)=cr(k)-cr(0)  
         enddo 
!	 
       		   Write (44,111) ij, aj1, aj2, aj3, aj4,  (ccrr(k),k=nn,11,-1)
       		   Write (44,112)                          (ccrr(k),k=10, 0,-1)		      
!
	ELSEIF(ICE_MODEL(1:4)=='ice5') THEN 	 
!
         do k=0,nn                   
                     ccrr(k)=cr(k)-cr(0)  
         enddo 
!	 
       		   Write (44,113) ij, aj1, ij, aj2,        (ccrr(k),k=nn,0,-1) 			
!
	Endif
!
!#---- Backup of the ice file  --------------------------------------------------- new new new
!
!
10          CONTINUE
!
	    close(10) 
!
            Write(*,*)'    - Read ',nel, '  elements from ', ice_model
            Write(*,*)'    - SH coefficients at time... '
!
            IF(ICE_MODEL(1:4)=='ice3') Open(44,file=ice_model//".bck",status='unknown')
!
!#---- Formats for ICE3G
111    format(i4,2f9.3,2f6.3,1x,8f5.0)
112    format(15f5.0)
!
!#---- Formats for ICE5G
113    format(i4,f10.4,1x,i4,f10.4,1x,22(f10.4))
114    format(i4,f10.4,1x,i4,f10.4,1x,27(f10.4)) 	
!
!#---- Formats for ANU05
115    format(i3,1x,3(f10.4,1x),32(f10.4,1x))
!
!
! ----------------------------------------------------------------------- 
!
! The ice spherical harmonics file are found broken in to small parts in separate
! files. Each has a size (1:JMAX,1:SLOT_SIZE), where SLOT_SIZE is fixed to 1000. 
!
! The ice breaking is performed by routine SHAPE_FACTORS_MOD.F90 
!
! Bounds on the broken files
!
 call ICE_BREAKER(NEL, NSLOTS, LO, HI)
!
 Write(*,*) '    - The ice model is broken in to ', NSLOTS, ' small piece(s) '
!
! The first index runs of the harmonics (J), 
! the second on the ice elements... 
 ARRAY(:,:)= DCMPLX(0.,0.) 
!
! This array is to store the harmonics of ice
! all time steps. Forst index is for harmonics
! (J), the second is for time (kk=0, 1, N, N+1)
! This statemente was BADLY wrong. Corrected on July 
! 19, 2010-. 
! 
 tttt(:,:) = DCMPLX(0.,0.)
!
!
!
!
! This big LOOP is over the broken ice files 
!
 DO 1000 KS=1, NSLOTS 
!
        call INTEGER_2_CHAR2(ks,labchar)
!
  	broken_ice_file(ks)='shice_broken_'//labchar//'.bin'
!
  	Open(100+ks,file=broken_ice_file(ks),status='unknown',form='unformatted')
!
  	write(*,'(i4,1x,i6,1x,i6,1x,a20)') ks, lo(ks), hi(ks), broken_ice_file(ks)
!
! - Reading each block 
!	do 12 islot =  1, slot_size 	
!  	do 12 ideg  =  1, jmax  
!		read(100+ks,*,end=122) ii, kk, array(ii,islot)
! 12 continue 
! 
! - Exists at EOF  
!122 continue
!
! - Reading each block
!
    read(100+ks) array
!
!
! - Reading the huge file "shice.bin" is now avoided... it was: 
!
!   open(7,file=filein,status='unknown',form='unformatted')
!   read(7) pppp
!   close(7) 
!
!
!#----- Computing the harmonic coefficients
!
 do k=0,nn+2     ! New Feb 2012 (it was "nn+1")
!
        call INTEGER_2_CHAR2(k,cj)
!
	if(k==0.or.k==nn+2)& 
	write(*,'(a7,1x,a2,a3,i4)') '     - ', cj, ' of', nn+1     
!
!
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP    PRIVATE(J,I) SHARED(TTTT,K,H,ARRAY,LO,HI,KS) SCHEDULE(GUIDED)
		do j=1,jmax 
        do i=lo(ks), hi(ks)
			 tttt(j,k)=tttt(j,k)+h(i,k) * array(j,i-(ks-1)*slot_size)  
        enddo
		enddo
!$OMP END PARALLEL DO
!
!				     
	open(1,file='sht'//cj//'.dat',status='unknown') 
		do j=1,jmax 
		write(1,*) j, real(tttt(j,k)), aimag(tttt(j,k)) 
		enddo
	close(1) 
!
 end do
!
! Endo of LOOP on the broken ice files 
1000 continue
!
!
!#----- Computing array "I"  
!
! This is the array of ice thickness variation. The first index 
! denotes the harmonic degree, the second is the time step.      
!
        write(*,*) '    - Computing array "I"' 
!
        iiii(:,:)=(0.,0.) 
!
        do k=0,nn +1   ! New Feb 2012 (it was 'nn') 
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(J,N) &
!$OMP     SHARED(IIII,K,TTTT) SCHEDULE(GUIDED)        
	      do j=1,jmax 
	        !iiii(j,k)=(0.,0.) 
		    do n=0, k
        	   iiii(j,k)=iiii(j,k)+(tttt(j,n+1)-tttt(j,n))
            end do
          end do
!$OMP END PARALLEL DO
        end do
!
!#----- Storing array "I" 
!
    write(*,*) '    - Array "I" is stored in the formatted file ', FILEOUT 
    open(7,file=fileout,status='unknown')  
 	write(7,*) iiii
 	close(7) 
!
!#------ Release memory
!
    deallocate( array )
    deallocate( TTTT, IIII )
    DEALLOCATE( H, CR, CCRR )     
!
 END PROGRAM SHICE 
!
!
!
!
!
!
