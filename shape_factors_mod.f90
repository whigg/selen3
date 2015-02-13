!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
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
!  This program computes the shape factors for the ice models in the library 
!  ./ICE-MODELS. Various methods are used - as explained below - according to
!  the geometrical features of each load. This routines collects in a single 
!  script various program units that where written "ad-hoc" for the various 
!  ice  models until version 2.5 of SELEN (i.e. those ending with "_C.F90"). 
!
!      Created by GS July 7 2008  "Intel port" Version 2.6 
!      Written by GS July 09 2008 ===
!      Modified GS July 24 for 'Zonal Disk' 
!      Also modified August 08 for ANU implementation on v. 2.7 
!      Touched March 2008 for the 'cap' ice load... 
!      Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
!      Reviewed GS & FC November 2009 - Porting under gfortran 
!      Modified May 2010 for the "ice breaking" method 
!      Reviewed GS on July 2010 for the implementation of the Greenland ice sheet
!      Reviewed GS on July 2010 for the implementation of the "small glaciers"
!      Retouched GS on July 2010 for the ANU ice model wit ice breaker 
! ***  Revised GS July 2010 - g95 - Double precision implementation
!      [Mainly for the call to ROMINT - see also HARMTOOLS.F90]
! ***  Revised GS & GR Nov 2010 - Implementation of MED1
! >>>  GS has fixed a big BUG in the ICE1 implementation... on February 24, 2011. 
! ***  Revised GS & GG May 2011 - Implementation of the Antarctic model "ANTA"
! ***  Revised GS & GG June 2011 - Implementation of the "TRANS" ice model 
! ***  Revised GS July 2011 - Implementation of the "RETT" ice model 
!      GS, July 2011: PROBLEMS with the LOAD SIZE in ROMINT  ???????????????????????????
!      GS Aug 011: After discussion with Daniele, MAXE is input /=0 in ROMINT 
!      GS SEP 011: Minor issues for the SLE benchmark for the RETT ice model 
!      FUNCLM === Modified to REAL*8 by GS september 30, 2011
!      February 2012: Introduced the "GRD" ice model type 
!      Last touch: February 22, 2012; the anniversary of Flo 
!      ANTA+GREE+GLAC models test for negative longitude === GS 01-07-2012
!      Revised DM Nov 2012 - performance improvements
!      Fixed a BUG by GS on MArch 2013, in the ICE1 script!
!      Jul 1, 2013; The "FLOR" ice for coupling with her models 
!      April 11, 2014. The BENOIT (BEN) files for the collaboration with LEGOS 
! 
!
! ---------------------------------------------------------------------------
!
!
!      INCLUDE "harmonics.f90"
       PROGRAM SHAPE_FACTORS 
       IMPLICIT NONE        
       INCLUDE "data.inc"
!
!
!#---- General declarations 
       INTEGER, PARAMETER :: LLMAX=256, JJMAX=(LLMAX+1)*(LLMAX+2)/2  
       INTEGER, PARAMETER :: SLOT_SIZE=1000 ! Consistent with ICE_BREAKER.F90
       INTEGER II, NELL, NSLOTS, LO(151), HI(151) 
       CHARACTER*2 LABCHAR 
       CHARACTER*20 broken_ice_file(151)
       REAL*8, ALLOCATABLE :: LEG(:), T(:,:)       
       INTEGER I, J, K, L, JJ, KS, IJ
       REAL*8 AJ
       REAL*8, ALLOCATABLE :: TETAC(:), LONGC(:), LATIC(:), ALFA(:)   
       CHARACTER*20 CJ
       REAL*8 COSDD, SINDD
       COMPLEX*16, ALLOCATABLE :: PPPP(:,:)              
!
! Old declaration for a huge array: COMPLEX*8 PPPP(JMAX,1:NEL)   
!     
       INTEGER NHEADER
!
!
!#---- Declarations for "RETT"
       INTEGER, PARAMETER :: NHY=20 
!REAL*4 AMPT
!REAL*4 longc4(nel), latic4(nel), tetac4(nel)
       REAL*8 AMPT
       REAL*8, allocatable :: longc4(:), latic4(:), tetac4(:)

!
!#---- Declarations for "TRANSIENT ICE"
       INTEGER, PARAMETER :: NHT=20 
       REAL*8 JUNK
!
!#---- Declarations for "SMALL GLACIERS"
       INTEGER, PARAMETER :: NHL=10 
!
!#---- Declarations for "GREENLAND"
       INTEGER, PARAMETER :: NHG=10 
!
!#---- Declarations for "ANTARCTICA"
       INTEGER, PARAMETER :: NHR=10
!
!#---- Declarations for "ALPS"
       INTEGER, PARAMETER :: NHA=4 
!
!#---- Declarations for "FLOR"
       INTEGER, PARAMETER :: NHF=4 
!
!#---- Declarations for "IJ05"
       INTEGER, PARAMETER :: NHI=18  
!
!#---- Declarations for "ICE5"
       INTEGER, PARAMETER :: NH5=24
       INTEGER, PARAMETER :: NLON=360, NLAT=180   ! it was: NLON=512, NLAT=256  
       INTEGER ILON, ILAT
       INTEGER, ALLOCATABLE :: ACTIVE(:,:)
       REAL*8, PARAMETER :: AMP5 = 1.0            ! it was: 0.70313 
       REAL*8, ALLOCATABLE :: PPLM(:)  
       REAL*8 ARG
!
!#---- Declarations for "IMED"
       INTEGER, PARAMETER :: NHM=20 
!
!#---- Declarations for "ICE3"
       INTEGER, PARAMETER :: NH3=20 
       COMPLEX*16, ALLOCATABLE :: ARMOY(:)          
!
!#---- Declarations for "ANU05"
       INTEGER, PARAMETER :: NHU=26 
!       
!#---- Declarations for "ICE1"
       INTEGER, PARAMETER :: NH1=20
       COMMON/AREA_DEG/LL,MM
       EXTERNAL FUNCLM      
!      REAL*4, PARAMETER :: AMP1 = 5. 
!      REAL*4  FUNCLM, AP, AM, S, T1, T2, ERR, CONV 
       REAL*8, PARAMETER :: AMP1 = 5. 
       REAL*8  FUNCLM, AP, AM, S, T1, T2, ERR, CONV, EPS

       REAL*8  DAP, DAM        
       INTEGER LJ, MJ, LL, MM, I1, I2, NNRO, MAXE
       INTEGER, PARAMETER :: MAXE_IN=100
       COMPLEX*16 C
!
!#---- Declarations for "DISK"
       INTEGER, PARAMETER :: NHD=20 
!
!#---- Declarations for "ICAP"
       INTEGER, PARAMETER :: NHC=20 
       REAL*8 CSI, COSA, COSLP0A, COSLP1A, COSLP2A, COSLM1A  
!
!#---- Declarations for "GRD"
       INTEGER, PARAMETER :: NHZ=10 
!
!#---- Declarations for "BEN"
       INTEGER, PARAMETER :: NHB=4 

!
!#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!      - Input file: ice_model (from 'data.inc') 
!      
!      - Output file:
       CHARACTER*20, PARAMETER :: FILEOUT='shicec.bin'
!#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!
! --- Allocate memory space
!
 ALLOCATE( LEG(0:LMAX+1), T(0:LMAX,NEL) )
 ALLOCATE( TETAC(nel), LONGC(nel), LATIC(nel), ALFA(nel) )
 ALLOCATE( PPPP(JMAX,SLOT_SIZE) )   
 ALLOCATE( longc4(nel), latic4(nel), tetac4(nel) )
 ALLOCATE( ACTIVE(NLON, NLAT) ) 
 ALLOCATE( PPLM(JMAX), ARMOY(JMAX) )  !
! 
!
! ----------------------------------------------------------------------- 
!
! The ice breaker routine breaks the ice spherical harmonics file into 
! small parts. Each has a size (1:JMAX,1:SLOT_SIZE), where SLOT_SIZE is
! fixed to 1000. 
!
!# Breaking the NEL-elements ice model in NSLOTS parts of length SLOT_SIZE
!
 call ICE_BREAKER(NEL, NSLOTS, LO, HI)
!
!# Creating a suitable number of files for hosting the parts... 	
!
 do 2 KS=1, NSLOTS 
!
       call INTEGER_2_CHAR2(KS,labchar)
!       
       broken_ice_file(KS)='shice_broken_'//labchar//'.bin'
!
       Open(100+KS,file=broken_ice_file(KS),status='unknown',form='unformatted')	
!
2 CONTINUE
!
 Write(*,*) '    - The ice model harmonics are broken in to ', & 
                   NSLOTS, ' small piece(s) '
!
! ----------------------------------------------------------------------- 
!
!
!#     -----------------------
!#====> An ice sheet like TRANS
!#     -----------------------
! Updated on May 25, 2010
!
	IF(ICE_MODEL(1:4)=='tran') THEN 
!
!#---- Reading the header lines  	    	
            open(10,file=ice_model,status='unknown') 
            do j=1, nht 
	         read(10,'(a20)') cj
	    enddo
!
!#----- Reads 'ice_model' 
            do i=1, nel
       			Read (10,*) ij, longc(i), latic(i), alfa(i) 
			do jj=1, 5
			Read (10,*) JUNK
			enddo
            enddo
       close(10)
       write(*,*)'    - Read ', nel, ' elements from ', ice_model       
!
!#---- Shape factors  	
        DO 71014 KS=1, NSLOTS 
!
        write(*,*)'    - Writing on file ', broken_ice_file(KS) 
!
!#--- Shape factors 
	    do i=lo(ks), hi(ks) 	! MODIFIED BY DM NOV-2012
           if(mod(i,100)==0)& 
	          write(*,*)'    - Shape factors for element ', i, ' of ', nel 
		   call pleg(lmax+1, 90.-alfa(i), leg) 
           call plmbar_mod(lmax, latic(i), pplm)		 
	  	   t(0,i)=(1.-cosdd(alfa(i)))/2.
	  	   do l=1, lmax
	     		t(l,i)= (-leg(l+1)+leg(l-1))/(2.*l+1.)/2.
	  	   enddo 
!
!#--- Shape factors by rotation 		 
        	do j=1, jmax
  		      arg=mj(j)*longc(i)
              pppp(j,i-(ks-1)*slot_size)=t(lj(j),i)*pplm(j)*cmplx(cosdd(arg),-sindd(arg))
            end do
         end do
!
         write(100+ks) pppp
!		    if(i.ge.lo(KS).and.i.le.hi(KS)) & 
!		    write(100+KS,*) j, i, t(lj(j),i)*pplm(j)*cmplx(cosdd(arg),-sindd(arg))
!		
1180 	continue 
71014   CONTINUE
!
!
!#      ----------------------  ---------------------- ---------------------- NEW NEW NEW 
!#====> An ice sheet like RECT  An ice sheet like RECT An ice sheet like RECT 
!#      ----------------------	---------------------- ---------------------- NEW NEW NEW 
!
! Updated on July 16, 2011
! Also updated on July 26 (in touch with VB)
!
	ELSEIF(ICE_MODEL(1:4)=='rett') THEN 
! 
!#---- Reading the header lines  	    	
            open(10,file=ice_model,status='unknown') 
            do j=1, nhy 
	         read(10,'(a20)') cj
	    enddo
!	    
!#---- Reading co-latitudes and longitudes of the ice elements  	    	
       	    do i=1, nel  
       		read (10,*) i1, longc4(i), latic4(i), AMPT
!
		if(longc4(i).le.0.0) longc4(i)=360.+longc4(i)
		tetac4(i)=90.0-latic4(i) 
!
! write(*,*) i, longc4(i), tetac4(i), AMPT
!
            do jj=1, 5
            Read (10,*) JUNK
            enddo
! 
            enddo
!
            write(*,*) '    - Read ', nel, ' elements from ', ice_model
	    close(10) 	
!
!#---- Shape factors  	    		
!
           conv=pi/180. 
!
	   DO 11018 KS=1, NSLOTS 
!
	   do i=lo(ks), hi(ks)
!
!write(*,*) i 
!
	       t1= (tetac4(i)-ampt/2.)*conv
	       t2= (tetac4(i)+ampt/2.)*conv 
!write(*,*) t1, t2
!
	   do j=1, jmax
	   
!if(mod(j,100)==0) write(*,*) i, j 
	   
	    	ll=lj(j)
	   	    mm=mj(j)  		   
!
	       if(mm==0) c = cmplx(ampt*conv,0.)  
	       if(mm/=0) then 
		       dap=float(mm)*(longc4(i)+ampt/2.)
		       dam=float(mm)*(longc4(i)-ampt/2.)
		       c = dcmplx(sindd(dap)-sindd(dam),cosdd(dap)-cosdd(dam))/float(mm)
		       endif
	       c=c/4.0/pi  
	       maxe=maxe_in
!	       call romint (s, err, 1.e-4, t1, t2, nnro, maxe, funclm)
!
	       call romint (s, err, 1.d-4, t1, t2, nnro, maxe, funclm)
!
	       If(maxe.gt.maxe_in) then 
			Write(* ,*) "Warning: the MAXE parameter from ROMINT exceeds the input value"
			Write(88,*) "Warning: the MAXE parameter from ROMINT exceeds the input value"			
	       Endif
!	       
! it was: pppp(j,i)= c*s
!
!		if(i.ge.lo(KS).and.i.le.hi(KS)) write(100+KS,*) j, i, c*s
        pppp(j,i-(ks-1)*slot_size)= c*s

       end do
       end do
       
       write(100+ks) pppp
       
11018	    CONTINUE
!
!
!
!
!
!#      ----------------------  ---------------------- ----------------------
!#====> An ice sheet like ICE1  An ice sheet like ICE1 An ice sheet like ICE1 
!#      ----------------------	---------------------- ----------------------
! Updated on May 25, 2010
!
	ELSEIF(ICE_MODEL(1:4)=='ice1') THEN 
! 
!#---- Reading the header lines  	    	
            open(10,file=ice_model,status='unknown') 
            do j=1, nh1 
	         read(10,'(a20)') cj
	    enddo
!	    
!#---- Reading co-latitudes and longitudes of the ice elements  	    	
       	    do i=1, nel   
       		read (10,*) i1, i2
       		tetac(i)=float(i1)
		longc(i)=float(i2)
            enddo
10          write(*,*) '    - Read ', nel, ' elements from ', ice_model
	    close(10) 	
!
!#---- Shape factors  	    		
!
           conv=pi/180. 
!
! FIXED by DM August 18, 2011 - MAXE initialization for ROMINT  
	   maxe=1000000
	   eps=1.d-6
!
	   DO 11011 KS=1, NSLOTS 
!
! FIXED by GS March 16, 2013 - it was "1, lo(ks), hi(ks)" !!!  
	   do  i=lo(ks), hi(ks) 
	       t1= (tetac(i)-amp1/2.)*conv
	       t2= (tetac(i)+amp1/2.)*conv 
!
	   do  j=1, jmax
	    	ll=lj(j)
	      	mm=mj(j)  
!	    	if(mod(j,20)==0)& 
!	   		Write(*,*)"    - Shape factors for degree ", j, " of ", jmax 		   
!
	       if(mm==0) c = cmplx(amp1*conv,0.)  
	       if(mm/=0) then 
		       dap=float(mm)*(longc(i)+amp1/2.)
		       dam=float(mm)*(longc(i)-amp1/2.)
		       c = dcmplx(sindd(dap)-sindd(dam),cosdd(dap)-cosdd(dam))/float(mm)
		       endif
	       c=c/4.0/pi  
!
	       call romint (s, err, EPS, t1, t2, nnro, maxe, funclm)
!
!write(*,*) ks, i, lo(ks), hi(ks), j 
!write(*,*) c
!write(*,*) s
!write(*,*) err
!write(*,*) eps
!write(*,*) t1, t2 
!write(*,*) nnro
!write(*,*) maxe
!
        pppp(j,i-(ks-1)*slot_size)=c*s
!
        end do
        end do
!        
        write(100+ks) pppp
!
11011	    CONTINUE
!
! #### Old piece of Code (before February 24, 2011) is below this line ------ 
!
!	    DO 11011 KS=1, NSLOTS 
!
!	    do 20 j=1, jmax
!            ll=lj(j)
!	    mm=mj(j)  
!            if(mod(j,20)==0)& 
!	    Write(*,*)"    - Shape factors for degree ", j, " of ", jmax    		    
!	    conv=pi/180. 
!
!	    do 20 i=1, nel 
!		if(mm==0) c = cmplx(amp1*conv,0.)  
!		if(mm/=0) then 
!		  	dap=float(mm)*(longc(i)+amp1/2.)
!		  	dam=float(mm)*(longc(i)-amp1/2.)
!		  	c = dcmplx(sindd(dap)-sindd(dam),cosdd(dap)-cosdd(dam))/float(mm)
!		  	endif
!		c=c/4.0/pi  
!		t1= (tetac(i)-amp1/2.)*conv
!		t2= (tetac(i)+amp1/2.)*conv 
!        	call romint (s, err, 1.e-4, t1, t2, nnro, maxe, funclm)
!		
! it was: pppp(j,i)= c*s
!
!                if(i.ge.lo(KS).and.i.le.hi(KS)) write(100+KS,*) j, i, c*s
!
!20          continue 		
!11011       CONTINUE
!
!
! #### Old piece of Code (before February 24, 2011) is below this line ------ 
!
!
!#     -----------------------
!#====> An ice sheet like ICE3
!#     -----------------------
! Updated on May 25, 2010
!
	ELSEIF(ICE_MODEL(1:4)=='ice3'.or.ICE_MODEL(1:4)=='imed') THEN 
!
!#---- Reading the header lines  	    	
            open(10,file=ice_model,status='unknown') 
 	    IF     (ICE_MODEL(1:4)=='ice3') then 
            do j=1, nh3 
	         read(10,'(a20)') cj
	    enddo
	    elseif (ICE_MODEL(1:4)=='imed') then
            do j=1, nhm 
	         read(10,'(a20)') cj
	    enddo	    
	    endif
!	    
!#---- Reading co-latitudes and longitudes of the ice elements  
       	    do i=1, nel  
       		Read (10,111) ij, tetac(i), longc(i), alfa(i)
       		Read (10,112) aj
       		latic(i)=90.-tetac(i)
       	    enddo
            write(*,*)'    - Read ',nel, '  elements from ', ice_model
	    close(10)
!
!#---- Degree-dependent factors
	    do 40 i=1, nel	  
	        call pleg(lmax+1, 90.-alfa(i), leg) 	  
	        t(0,i)=(1.-cosdd(alfa(i)))/2.
	        do 40 l=1, lmax
	     	        t(l,i)= (-leg(l+1)+leg(l-1))/(2.*l+1.)/2.
40          continue
!
!#---- Shape factors 
!
	    DO 11012 KS=1, NSLOTS 
! 
	    	write(*,*)'    - Writing on file ', broken_ice_file(KS) 
!	
	    	do i=lo(ks), hi(ks) 
!
            	if(mod(i,250)==0)& 
		write(*,*)'    - Shape factors for element ', i, ' of ', nel
!
	    	call harmo(lmax, longc(i), latic(i), armoy) 
!
            	do  j=1, jmax 
!
!                if(i.ge.lo(KS).and.i.le.hi(KS)) & 
!		write(100+KS,*) j, i, t(lj(j),i)*conjg(armoy(j))
        pppp(j,i-(ks-1)*slot_size)=t(lj(j),i)*conjg(armoy(j))
!
! it was: pppp(j,i) = t(lj(j),i)*conjg(armoy(j))  
!
        end do
        end do

        write(100+ks) pppp

11012       CONTINUE	
!
!
!#     ------------------------
!#====> An ice sheet like ANU05
!#     ------------------------
! Updated on May 25, 2010
!
	ELSEIF(ICE_MODEL(1:4)=='anu0') THEN 
!
!#---- Reading the header lines  	    	
            open(10,file=ice_model,status='unknown') 
            do j=1, nhu 
	         read(10,'(a20)') cj
	    enddo
!	    
!#---- Reading co-latitudes and longitudes of the ice elements  
       	    do i=1, nel  
       		Read (10,115) ij, longc(i), latic(i), alfa(i)
       		tetac(i)=90.-latic(i)
       	    enddo
            write(*,*)'    - Read ',nel, '  elements from ', ice_model
	    close(10)
!
!#---- Degree-dependent factors
	    do 41 i=1, nel	  
	        call pleg(lmax+1, 90.-alfa(i), leg) 	  
	        t(0,i)=(1.-cosdd(alfa(i)))/2.
	        do 41 l=1, lmax
	     	        t(l,i)= (-leg(l+1)+leg(l-1))/(2.*l+1.)/2.
41          continue
!
!#---- Shape factors  	
	    DO 11013 KS=1, NSLOTS 
!
	    	write(*,*)'    - Writing on file ', broken_ice_file(KS) 
!
	    	do  i = lo(ks), hi(ks) 
!
	    	if(mod(i,10000)==0)& 
		write(*,*)'    - Shape factors for element ', i, ' of ', nel 
!
            	call harmo(lmax, longc(i), latic(i), armoy) 
!
            	do j=1, jmax 
!
!                if(i.ge.lo(KS).and.i.le.hi(KS)) & 
!		write(100+KS,*) j, i, t(lj(j),i)*conjg(armoy(j))
		pppp(j,i-(ks-1)*slot_size) = t(lj(j),i)*conjg(armoy(j))
!
! it was: pppp(j,i) = t(lj(j),i)*conjg(armoy(j))  
!
        end do
        end do
        
        write(100+ks) pppp
!
11013       CONTINUE	
!
!
!#     ----------------------
!#====> An ice sheet like ICE5
!#     ----------------------
!
	ELSEIF(ICE_MODEL(1:4)=='ice5') THEN 
!
!#---- Reading the header lines  	    	
            open(10,file=ice_model,status='unknown') 
            do j=1, nh5 
	         read(10,'(a20)') cj
	    enddo
!	    
!#---- Reading longitudes and latitudes of the ice elements            
	    active(:,:)=0 	
       	    do 60 i=1, nel
       		Read (10,113) ilon, longc(i), ilat, latic(i)
       		active(ilon,ilat)=1
60          continue
            write(*,*)'    - Read ', nel, '  elements from ', ice_model
	    close(10) 
!
!#---- Shape factors 
	    i=0 
	    ks=1
!
	    do 70 ilat=1, nlat  
	    do 70 ilon=1, nlon	
!
	    IF(ACTIVE(ILON,ILAT)==1) THEN 
! 	
	    i=i+1 
!
	    if(mod(i,5000)==0)& 
	    write(*,*)'    - Shape factors for element ', i, ' of ', nel 
!
! --- Half-amplitude of the disc load with equal area of the quadrilateral 	
	    alfa(i)=acos(1.-sindd(90.-latic(i))*sindd(amp5/2.)*amp5/180.)*180./pi
!
! --- New Plm's are computed only when latitude changes		 
	    if(i==1.or.(i>=2.and.latic(i)/=latic(i-1))) then 
	       call pleg(lmax+1, 90.-alfa(i), leg) 
               call plmbar_mod(lmax, latic(i), pplm)
            endif  
!
! --- Degree-dependent factors	 
	    t(0,i)=(1.-cosdd(alfa(i)))/2.
	    do l = 1, lmax
	     	 t(l,i)= (-leg(l+1)+leg(l-1))/(2.*l+1.)/2.
	    enddo 
!
!#--- Shape factors by rotation 		 
!
!#--- Filling the broken ice files to avoid memory issues when the 
!     maximum harmonic degree is "large"  ------- May 2010 ------- 
!
! If ice element "i" is found in some slot, the slot file is updated 
! with the shape factor of degree "j" for the ice element "i"  - The 
! large array "PPPP" is not necessary anymore.  
!  
! ---- OLD CODE ----   MODIFIED BY DM Nov 2012
!            do k=1, nslots 	    
!
!                do j=1, jmax
!                arg=mj(j)*longc(i)
!
! it was:       pppp(j,i) = t(lj(j),i)*pplm(j)*cmplx(cosdd(arg),-sindd(arg))
!
!		    if(i.ge.lo(k).and.i.le.hi(k)) then 
!		    write(100+k,*) j, i, t(lj(j),i)*pplm(j)*cmplx(cosdd(arg),-sindd(arg))
!		    Endif
!                enddo		   		   
!	    enddo		 
!		 
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP    PRIVATE(J,ARG) SHARED(LONGC,PPPP,T,I,PPLM,KS) &
!$OMP        SCHEDULE(GUIDED)
     do j=1, jmax
         arg=mj(j)*longc(i)
         pppp(j,i-(ks-1)*slot_size)=t(lj(j),i)*pplm(j)*cmplx(cosdd(arg),-sindd(arg))
     end do
!$OMP END PARALLEL DO
!
!
!
! --- Detect a 'slot change'
!
     if( i.eq.hi(ks) ) then
            write(100+ks) pppp      ! Dump the current SLOT
            ks = ks + 1             ! Increment slot counter
     end if
!       
!                    
!
!
!
     ENDIF ! active elements 
!
70 continue   
!	
!###############################################################
!
!	
!
!#     -----------------------
!#====> An ice sheet like FLO 
!#     -----------------------
! Updated on Jul 1, 2013
!
	ELSEIF(ICE_MODEL(1:4)=='flor') THEN 
!
!#---- Reading the header lines  	    	
            open(10,file=ice_model,status='unknown') 
            do j=1, nhf 
	         read(10,'(a20)') cj
	    enddo
!
!#----- Reads 'ice_model' 
            do i=1, nel
       		Read (10,*) ij, longc(i), latic(i), alfa(i) 
            enddo
       close(10)
       write(*,*)'    - Read ', nel, ' elements from ', ice_model       
!
!#---- Shape factors  	
        DO 11814 KS=1, NSLOTS 
!
        write(*,*)'    - Writing on file ', broken_ice_file(KS) 
!
!#--- Shape factors 
	do  i=lo(ks), hi(ks)
        if(mod(i,100)==0)& 
	write(*,*)'    - Shape factors for element ', i, ' of ', nel 
		call pleg(lmax+1, 90.-alfa(i), leg) 
        	call plmbar_mod(lmax, latic(i), pplm)		 
	  	t(0,i)=(1.-cosdd(alfa(i)))/2.
	  	do l=1, lmax
	     		t(l,i)= (-leg(l+1)+leg(l-1))/(2.*l+1.)/2.
	  	enddo 
!
!#--- Shape factors by rotation 		 
        	do j=1, jmax
		    arg=mj(j)*longc(i)
!
!		    if(i.ge.lo(KS).and.i.le.hi(KS)) & 
!		    write(100+KS,*) j, i, t(lj(j),i)*pplm(j)*cmplx(cosdd(arg),-sindd(arg))
                  pppp(j,i-(ks-1)*slot_size)= t(lj(j),i)*pplm(j)*cmplx(cosdd(arg),-sindd(arg))
            end do
!
! it was: pppp(j,i) = t(lj(j),i)*pplm(j)*cmplx(cosdd(arg),-sindd(arg))
!		
    end do
    
    write(100+ks) pppp
    
11814   CONTINUE
!
!
!








!
!###############################################################
!
!
!#     -----------------------
!#====> An ice sheet like ALPS
!#     -----------------------
! Updated on May 25, 2010
!
	ELSEIF(ICE_MODEL(1:4)=='alps') THEN 
!
!#---- Reading the header lines  	    	
            open(10,file=ice_model,status='unknown') 
            do j=1, nha 
	         read(10,'(a20)') cj
	    enddo
!
!#----- Reads 'ice_model' 
            do i=1, nel
       		Read (10,*) ij, longc(i), latic(i), alfa(i) 
!               WARNING - After the update of April 22 2008, ALFA(I) = 0.09 deg,  
!               consistent with the number of pixels used to build the model -  
!               alfa(i)=0.09
            enddo
       close(10)
       write(*,*)'    - Read ', nel, ' elements from ', ice_model       
!
!#---- Shape factors  	
        DO 11014 KS=1, NSLOTS 
!
        write(*,*)'    - Writing on file ', broken_ice_file(KS) 
!
!#--- Shape factors 
	do  i=lo(ks), hi(ks)
        if(mod(i,100)==0)& 
	write(*,*)'    - Shape factors for element ', i, ' of ', nel 
		call pleg(lmax+1, 90.-alfa(i), leg) 
        	call plmbar_mod(lmax, latic(i), pplm)		 
	  	t(0,i)=(1.-cosdd(alfa(i)))/2.
	  	do l=1, lmax
	     		t(l,i)= (-leg(l+1)+leg(l-1))/(2.*l+1.)/2.
	  	enddo 
!
!#--- Shape factors by rotation 		 
        	do j=1, jmax
		    arg=mj(j)*longc(i)
!
!		    if(i.ge.lo(KS).and.i.le.hi(KS)) & 
!		    write(100+KS,*) j, i, t(lj(j),i)*pplm(j)*cmplx(cosdd(arg),-sindd(arg))
                  pppp(j,i-(ks-1)*slot_size)= t(lj(j),i)*pplm(j)*cmplx(cosdd(arg),-sindd(arg))
            end do
!
! it was: pppp(j,i) = t(lj(j),i)*pplm(j)*cmplx(cosdd(arg),-sindd(arg))
!		
    end do
    
    write(100+ks) pppp
    
11014   CONTINUE
!
!
!
!#     -----------------------
!#====> An ice sheet like "GRD"
!#     -----------------------
! February 2012 from Giorgio Spada ("ice2sea")
!
	ELSEIF(ICE_MODEL(1:3)=='grd') THEN 
!
!#---- Reading the header lines  	    	
            open(10,file=ice_model,status='unknown') 
            do j=1, nhz 
	         read(10,'(a20)') cj
	    enddo
!
!#----- Reads 'ice_model' 
            do i=1, nel
       		Read (10,*) ij, longc(i), latic(i), alfa(i) 
            enddo
       close(10)
       write(*,*)'    - Read ', nel, ' elements from ', ice_model       
!
!#---- Shape factors  	
        DO 11714 KS=1, NSLOTS 
!
        write(*,*)'    - Writing on file ', broken_ice_file(KS) 
!
!#--- Shape factors 
!       do 807 i=1, nel    
	do i = lo(ks), hi(ks)	
        if(mod(i,100)==0)& 
	write(*,*)'    - Shape factors for element ', i, ' of ', nel 
		call pleg(lmax+1, 90.-alfa(i), leg) 
        	call plmbar_mod(lmax, latic(i), pplm)		 
	  	t(0,i)=(1.-cosdd(alfa(i)))/2.
	  	do l=1, lmax
	     		t(l,i)= (-leg(l+1)+leg(l-1))/(2.*l+1.)/2.
	  	enddo 
!
!#--- Shape factors by rotation 		 
        	do j=1, jmax
		    arg=mj(j)*longc(i)
!
!		    if(i.ge.lo(KS).and.i.le.hi(KS)) & 
!		    write(100+KS,*) j, i, t(lj(j),i)*pplm(j)*cmplx(cosdd(arg),-sindd(arg))
                  pppp(j,i-(ks-1)*slot_size)= t(lj(j),i)*pplm(j)*cmplx(cosdd(arg),-sindd(arg))
            end do
    end do
    
    write(100+ks) pppp

!
! it was: pppp(j,i) = t(lj(j),i)*pplm(j)*cmplx(cosdd(arg),-sindd(arg))
!		
11714   CONTINUE
!
!
!
!-----------------
!
!#     -----------------------
!#====> An ice sheet like "BENO"
!#     -----------------------
!
! April 2012 from Giorgio Spada ("collaboration with LEGOS")
!
	ELSEIF(ICE_MODEL(1:4)=='beno') THEN 
!
!#---- Reading the header lines  	    	
            open(10,file=ice_model,status='unknown') 
            do j=1, nhb 
	         read(10,'(a20)') cj
	    enddo
!
!#----- Reads 'ice_model' 
            do i=1, nel
       		Read (10,*) ij, longc(i), latic(i), alfa(i) 
!               write (*,*) ij, longc(i), latic(i), alfa(i) 
            enddo
       close(10)
       write(*,*)'    - Read ', nel, ' elements from ', ice_model       
!
!#---- Shape factors  	
        DO 1114 KS=1, NSLOTS 
!
        write(*,*)'    - Writing on file ', broken_ice_file(KS) 
!
!#--- Shape factors 
!       do 807 i=1, nel    
	do i = lo(ks), hi(ks)	
        if(mod(i,100)==0)& 
	write(*,*)'    - Shape factors for element ', i, ' of ', nel 
		call pleg(lmax+1, 90.-alfa(i), leg) 
        	call plmbar_mod(lmax, latic(i), pplm)		 
	  	t(0,i)=(1.-cosdd(alfa(i)))/2.
	  	do l=1, lmax
	     		t(l,i)= (-leg(l+1)+leg(l-1))/(2.*l+1.)/2.
	  	enddo 
!
!#--- Shape factors by rotation 		 
        	do j=1, jmax
		    arg=mj(j)*longc(i)
!
!		    if(i.ge.lo(KS).and.i.le.hi(KS)) & 
!		    write(100+KS,*) j, i, t(lj(j),i)*pplm(j)*cmplx(cosdd(arg),-sindd(arg))
                  pppp(j,i-(ks-1)*slot_size)= t(lj(j),i)*pplm(j)*cmplx(cosdd(arg),-sindd(arg))
            end do
    end do
    
    write(100+ks) pppp

!
! it was: pppp(j,i) = t(lj(j),i)*pplm(j)*cmplx(cosdd(arg),-sindd(arg))
!		
1114   CONTINUE
!
!




!-----------------
!
!#     ----------------------- ------------------------------------------
!#====> An ice sheet like GREENLAND or "Antarctica" or "SMALL ice sheets"
!#     ----------------------- ------------------------------------------
!
! Updated on May 25, 2010
! Refreshed on JULY 2010- many versiona tested. The behavior of the
! ice breaker appeared WRONG because of a subtle mistake in SHICE_MOD.F90. 
! Now fixed. 
! Implemented on July 24, 2010 also for the small glaciers
! *** Revised GS & GG May 2011 - Implementation of the Antarctic model "ANTA" 
!
	ELSEIF(ICE_MODEL(1:4)=='gree'.or.ICE_MODEL(1:4)=='glac'.or.ICE_MODEL(1:4)=='anta') THEN 
!
        if(ICE_MODEL(1:4)=='gree')nheader=nhg 
        if(ICE_MODEL(1:4)=='glac')nheader=nhl 
        if(ICE_MODEL(1:4)=='anta')nheader=nhr 
!
!#---- Reading the header lines  	    	
            open(10,file=ice_model,status='unknown') 
            do j=1, nheader 
	         read(10,'(a20)') cj
	    enddo
!
!#----- Reads 'ice_model' 
            do i=1, nel
       		Read (10,*) ij, longc(i), latic(i), alfa(i) 
		if(longc(i).le.0.0) longc(i)=longc(i)+360.0
            enddo
       close(10)
       write(*,*)'    - Read ', nel, ' elements from ', ice_model       
!
!#---- Degree-dependent factors
! Assuming a priori different sizes for the ice emelents...
!
	    do 141 i=1, nel	  
	        call pleg(lmax+1, 90.-alfa(i), leg) 	  
	        t(0,i)=(1.-cosdd(alfa(i)))/2.
	        do 141 l=1, lmax
	     	        t(l,i)= (-leg(l+1)+leg(l-1))/(2.*l+1.)/2.
141          continue
!
!
!#---- Shape factors  	
	    DO 11713 KS=1, NSLOTS 
!
	    	write(*,*)'    - Writing on file ', broken_ice_file(KS) 
!
		    do i = lo(ks), hi(ks) 
!
            	call harmo(lmax, longc(i), latic(i), armoy) 
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(PPPP,T,ARMOY,I,KS) PRIVATE(J) &
!$OMP       SCHEDULE(GUIDED)
            	do  j=1, jmax 
		!write(100+KS,*) j, i, t(lj(j),i)*conjg(armoy(j))
		            pppp(j,i-(ks-1)*slot_size) = t(lj(j),i)*conjg(armoy(j))
   	            end do
!$OMP END PARALLEL DO	    
            end do
            
            write(100+ks) pppp
           


11713       CONTINUE
!
!
!#      ----------------------
!#====> An ice sheet like IJ05
!#      ----------------------
! Updated on May 25, 2010 ***********
!
	ELSEIF(ICE_MODEL(1:4)=='ij05') THEN 
!
!#---- Reading the header lines  	    	
            open(10,file=ice_model,status='unknown') 
            do j=1, nhi 
	         read(10,'(a20)') cj
	    enddo
!
! --- Reads 'ice_model'  
       do i=1, nel  
       		Read (10,*) ij, longc(i), latic(i), alfa(i)
       		Read (10,*) aj 
       enddo
       close(10)
       WRITE(*,*)'    - Read ',nel, ' ice elements from ', ice_model
!
!#--- Degree-dependent factors	  
	do i=1,nel
	  call pleg(lmax+1, 90.-alfa(i), leg) 	  
	  t(0,i)=(1.-cosdd(alfa(i)))/2.
	  	do l=1, lmax
	     		t(l,i)= (-leg(l+1)+leg(l-1))/(2.*l+1.)/2.
	  	enddo 	    
        enddo
!
!#--- Shape factors 
        DO 11015 KS=1, NSLOTS 
! 
	    	write(*,*)'    - Writing on file ', broken_ice_file(KS) 
!
	do  i=lo(ks), hi(ks)
		if(mod(i,100)==0)& 
		write(*,*)'    - Shape factors for element ', i, ' of ', nel 
	        call harmo(lmax, longc(i), latic(i), armoy) 
			do j=1, jmax 
!
            !    	if(i.ge.lo(KS).and.i.le.hi(KS)) & 
			!write(100+KS,*) j, i, t(lj(j),i)*conjg(armoy(j))
			     pppp(j,i-(ks-1)*slot_size) = t(lj(j),i)*conjg(armoy(j))  
!
! it was: pppp(j,i) = t(lj(j),i)*conjg(armoy(j))  
!
            end do
    end do

    write(100+ks) pppp

11015   CONTINUE
!
!
!#     -----------------------
!#====> An ice sheet like DISK
!#     -----------------------
! Updated on May 25, 2010 ***********
!
	ELSEIF(ICE_MODEL(1:4)=='disk') THEN 
!
!#---- Reading the header lines  	    	
            open(10,file=ice_model,status='unknown') 
            do j=1, nhd 
	         read(10,'(a20)') cj
	    enddo
!	    
!#---- Reading co-latitudes and longitudes of the ice elements  
       	    do i=1, nel  
       		Read (10,*) ij, tetac(i), longc(i), alfa(i)
       		Read (10,*) aj
       		latic(i)=90.-tetac(i)
       	    enddo
            write(*,*)'    - Read ',nel,'  elements from ', ice_model
	    close(10)
!
!#---- Degree-dependent factors
	    do 90 i=1, nel	  
	        call pleg(lmax+1, 90.-alfa(i), leg) 	  
	        t(0,i)=(1.-cosdd(alfa(i)))/2.
	        do 90 l=1, lmax
	     	        t(l,i)= (-leg(l+1)+leg(l-1))/(2.*l+1.)/2.
90          continue
!
!#---- Shape factors  	
        DO 11016 KS=1, NSLOTS 
!
	    	write(*,*)'    - Writing on file ', broken_ice_file(KS) 
!
        do  i=lo(ks),hi(ks) 
		if(mod(i,100)==0)& 
		write(*,*)'    - Shape factors for element ', i, ' of ', nel 
	        call harmo(lmax, longc(i), latic(i), armoy) 
		do j=1, jmax 
!
!                	if(i.ge.lo(KS).and.i.le.hi(KS)) & 
!			write(100+KS,*) j, i, t(lj(j),i)*conjg(armoy(j))
            pppp(j,i-(ks-1)*slot_size) = t(lj(j),i)*conjg(armoy(j))
        end do
        end do
!
        write(100+ks) pppp
!
!
!
! it was: pppp(j,i) = t(lj(j),i)*conjg(armoy(j))  
!
11016   CONTINUE
!
!
!#     ----------------------
!#====> An ice sheet like CAP
!#     ----------------------
! Updated on May 25, 2010 ***********

	ELSEIF(ICE_MODEL(1:4)=='icap') THEN 
!
!#---- Reading the header lines  	    	
            open(10,file=ice_model,status='unknown') 
            do j=1, nhd 
	         read(10,'(a20)') cj
	    enddo
!	    
!#---- Reading co-latitudes and longitudes of the ice elements  
       	    do i=1, nel  
       		Read (10,*) ij, tetac(i), longc(i), alfa(i)
       		Read (10,*) aj
       		latic(i)=90.-tetac(i)
       	    enddo
            write(*,*)'    - Read ',nel, '  elements from ', ice_model
	    close(10)
!
!#---- Degree-dependent factors
	    do 94 i=1, nel	  
!			  
	        t(0,i)=(1.-cosdd(alfa(i)))/3.
!				
		do 94 l = 1, lmax
!		
		       cosa    = cosdd(alfa(i))
		       coslp0a = cosdd((l+0.)*alfa(i))
		       coslp1a = cosdd((l+1.)*alfa(i))
		       coslp2a = cosdd((l+2.)*alfa(i))
		       coslm1a = cosdd((l-1.)*alfa(i))		    
!               
	        csi    = -(0.75/(1.-cosa)/(1.-cosa))*( (coslp1a-coslp2a)/(1.*l+1.5) - & 
		                                       (coslm1a-coslp0a)/(1.*l-0.5) )		 		
!		
		t(l,i) = (1.-cosa)*csi/3./(2.*l+1.)
!
94          continue
!
!
!#---- Shape factors  	
        DO 11017 KS=1, NSLOTS 
!
	    	write(*,*)'    - Writing on file ', broken_ice_file(KS) 
!
!#---- Shape factors  	
	    do  i=1, lo(ks), hi(ks) 
		if(mod(i,100)==0)& 
	        write(*,*)'    - Shape factors for element ', i, ' of ', nel
	        call harmo(lmax, longc(i), latic(i), armoy) 
		do j=1, jmax 
!
!                	if(i.ge.lo(KS).and.i.le.hi(KS)) & 
!			write(100+KS,*) j, i, t(lj(j),i)*conjg(armoy(j))
		    pppp(j,i-(ks-1)*slot_size) = t(lj(j),i)*conjg(armoy(j))
        end do
        end do
!
! it was: pppp(j,i) = t(lj(j),i)*conjg(armoy(j))  
!
 
        write(100+ks) pppp

11017    CONTINUE 
!
	ENDIF
!
!#      --------------------------
!#====> End of available models...
!#      --------------------------
!
	Write(*,*)'    - Computed ', nel, ' shape factors'
!
!#---- Format for ICE3 elements...   
111    FORMAT(I4,2F9.3,2F6.3,1X,8F5.0)
112    FORMAT(15F5.0) 
!
!#---- Format for ICE5G elements...   
113    FORMAT(I4,F10.4,1X,I4,F10.4,1X)
!
!#---- Formats for ANU05
115    format(i3,1x,3(f10.4,1x),32(f10.4,1x))
!
!
       do 222 KS=1, NSLOTS 
             CLOSE(KS)
 222   Continue
!
! -----------------------------------------------------------------------------
!
! In this last part of the code, the array PPPP was written on file "shice.bin"
! It is clear now that this is not any more necessary since the ice harmonics
! are broken into a suite of files with names: 
!
! 'shice_broken_**.dat', where '**' = 01, 02, .... 
!
! Write(*,*)'    - Writing the shape factors on "shicec.bin"'
! open (7,file=fileout,status='unknown')
! write(7,*) pppp
! close(7) 
! Write(*,*)'    - done!'
!
! -----------------------------------------------------------------------------
!
!
 DEALLOCATE( LEG, T )
 DEALLOCATE( TETAC, LONGC, LATIC, ALFA )
 DEALLOCATE( PPPP )   
 DEALLOCATE( longc4, latic4, tetac4 )
 DEALLOCATE( ACTIVE ) 
 DEALLOCATE( PPLM, ARMOY )  
! 
 END PROGRAM SHAPE_FACTORS 
!
!
!
!
!
!
