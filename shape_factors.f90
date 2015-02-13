!
! Created by GS July 7 2008  "Intel port" Version 2.6 
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
! *** Revised GS July 2010 - g95 - Double precision implementation
!     [Mainly for the call to ROMINT - see also HARMTOOLS.F90]
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
!  This program computes the shape factors for the ice models in the library 
!  ./ICE-MODELS. Various methods are used - as explained below - according to
!  the geometrical features of the load. This routines collects in a single 
!  script various program units that where written "ad-hoc" for the various 
!  ice  models until version 2.5 of SELEN (i.e. those ending with "_C.F90"). 
!
!  === Written by GS July 09 2008 ===
!      Modified GS July 24 for 'Zonal Disk' 
!      Also modified August 08 for ANU implementation on v. 2.7 
!      Touched March 2008 for the 'cap' ice load... 
! ---------------------------------------------------------------------------
!
!
!      INCLUDE "harmonics.f90"
       PROGRAM SHAPE_FACTORS 
       IMPLICIT NONE        
       INCLUDE "data.inc"
!
!#---- General declarations 
       INTEGER, PARAMETER :: LLMAX=256, JJMAX=(LLMAX+1)*(LLMAX+2)/2  
!
       COMPLEX*16 PPPP(JMAX,1:NEL)   
!
       REAL*8 LEG(0:LMAX+1), T(0:LMAX,NEL) 
       INTEGER I, J, K, L, IJ
       REAL*8 AJ, TETAC(NEL), LONGC(NEL), LATIC(NEL), ALFA(NEL)   
       CHARACTER*20 CJ
       REAL*8 COSDD, SINDD
!
!#---- Declarations for "ALPS"
       INTEGER, PARAMETER :: NHA=4 
!
!#---- Declarations for "IJ05"
       INTEGER, PARAMETER :: NHI=18  
!
!#---- Declarations for "ICE5"
       INTEGER, PARAMETER :: NH5=24
       INTEGER, PARAMETER :: NLON=360, NLAT=180   ! it was: NLON=512, NLAT=256  
       INTEGER ILON, ILAT, ACTIVE(NLON, NLAT)    
       REAL*8, PARAMETER :: AMP5 = 1.0            ! it was: 0.70313 
       REAL*8 PPLM(JJMAX)  
       REAL*8 ARG
!
!#---- Declarations for "ICE3"
       INTEGER, PARAMETER :: NH3=20 
       COMPLEX*16 ARMOY(JJMAX)          
!
!#---- Declarations for "ANU05"
       INTEGER, PARAMETER :: NHU=26 
!
!#---- Declarations for "ICE1"
       INTEGER, PARAMETER :: NH1=20
       COMMON/AREA_DEG/LL,MM
       EXTERNAL FUNCLM      
       REAL*4, PARAMETER :: AMP1 = 5. 
       REAL*4  FUNCLM, AP, AM, S, T1, T2, ERR, CONV 
       REAL*8  DAP, DAM 
       INTEGER LJ, MJ, LL, MM, I1, I2, NNRO, MAXE
       COMPLEX*16 C
!
!#---- Declarations for "DISK"
       INTEGER, PARAMETER :: NHD=20 
!
!#---- Declarations for "ICAP"
       INTEGER, PARAMETER :: NHC=20 
       REAL*8 CSI, COSA, COSLP0A, COSLP1A, COSLP2A, COSLM1A  
!
!
!
!#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!      - Input file: ice_model (from 'data.inc') 
!      
!      - Output file:
       CHARACTER*20, PARAMETER :: FILEOUT='shicec.bin'
!#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!
!
!
!#      ----------------------
!#====> An ice sheet like ICE1 
!#      ----------------------
!
	IF(ICE_MODEL(1:4)=='ice1') THEN 
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
		longc(i) =float(i2)
            enddo
10          write(*,*) & 
	    '    - Read ', nel, & 
	    ' elements from ', ice_model
	    close(10) 	
!
!#---- Shape factors  	    	
	    do 20 j=1, jmax
            ll=lj(j)
	    mm=mj(j)  
            if(mod(j,20)==0) Write(*,*)& 
	    "    - Degree ", j, " of ", jmax    		    
	    conv=pi/180. 
!
	    do 20 i=1, nel 
		if(mm==0) c = cmplx(amp1*conv,0.)  
		if(mm/=0) then 
		  	dap=float(mm)*(longc(i)+amp1/2.)
		  	dam=float(mm)*(longc(i)-amp1/2.)
		  	c = cmplx(sindd(dap)-sindd(dam),cosdd(dap)-cosdd(dam))/float(mm)
		  	endif
		c=c/4.0/pi  
		t1= (tetac(i)-amp1/2.)*conv
		t2= (tetac(i)+amp1/2.)*conv 
        	call romint (s, err, 1.e-4, t1, t2, nnro, maxe, funclm)
		pppp(j,i)= c*s
20          continue 		
!
!
!
!#     -----------------------
!#====> An ice sheet like ICE3
!#     -----------------------
!
	ELSEIF(ICE_MODEL(1:4)=='ice3') THEN 
!
!
!#---- Reading the header lines  	    	
            open(10,file=ice_model,status='unknown') 
            do j=1, nh3 
	         read(10,'(a20)') cj
	    enddo
!	    
!#---- Reading co-latitudes and longitudes of the ice elements  
       	    do i=1, nel  
       		Read (10,111) ij, tetac(i), longc(i), alfa(i)
       		Read (10,112) aj
       		latic(i)=90.-tetac(i)
       	    enddo
            write(*,*)&
	    '    - Read ',nel, &
	    '  elements from ', ice_model
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
	    do 50 i=1, nel 
		if(mod(i,250)==0)write(*,*)&
            '    - ', i, ' of ', nel 
	        call harmo(lmax, longc(i), latic(i), armoy) 
		do 50 j=1, jmax 
	            pppp(j,i) = t(lj(j),i)*conjg(armoy(j))  
50	    continue
!	
!
!#     ------------------------
!#====> An ice sheet like ANU05
!#     ------------------------
!
	ELSEIF(ICE_MODEL(1:4)=='anu0') THEN 
!
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
            write(*,*)&
	    '    - Read ',nel, &
	    '  elements from ', ice_model
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
	    do 51 i=1, nel 
		if(mod(i,1000)==0)write(*,*)&
            '    - ', i, ' of ', nel 
	        call harmo(lmax, longc(i), latic(i), armoy) 
		do 51 j=1, jmax 
	            pppp(j,i) = t(lj(j),i)*conjg(armoy(j))  
51	    continue
!	
!
!
!#     ----------------------
!#====> An ice sheet like ICE5
!#     ----------------------
!
	ELSEIF(ICE_MODEL(1:4)=='ice5') THEN 
!
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
            write(*,*) &
	    '    - Read ', nel, &
	    '  elements from ', ice_model
	    close(10) 
!
!#---- Shape factors 
	    i=0 
!
	    do 70 ilat=1, nlat  
	    do 70 ilon=1, nlon	
	    if(active(ilon,ilat)==1) then 
! 	
	    i=i+1 
!
	    if(mod(i,5000)==0)write(*,*)&
	    '    -  ', & 
	    i, ' of ', nel 
!
! --- Half-amplitude of the disc load with equal area of the quadrilateral 	
	    alfa(i)=acos(1.-sindd(90.-latic(i))*sindd(amp5/2.)*amp5/180.)*180./pi
!write(*,*) i, alfa(i)
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
            do j=1, jmax
		    arg=mj(j)*longc(i)
	            pppp(j,i) = t(lj(j),i)*pplm(j)*cmplx(cosdd(arg),-sindd(arg))
	    enddo		 
!		 
	    endif	
!
70          continue 
!  
!	
!
!#     -----------------------
!#====> An ice sheet like ALPS
!#     -----------------------
!
	ELSEIF(ICE_MODEL(1:4)=='alps') THEN 
!
!#---- Reading the header lines  	    	
            open(10,file=ice_model,status='unknown') 
            do j=1, nha 
	         read(10,'(a20)') cj
	    enddo
!
!
!#----- Reads 'ice_model' 
            do i=1, nel
       		Read (10,*) ij, longc(i), latic(i), alfa(i) 
!               WARNING - After the update of April 22 2008, ALFA(I) = 0.09 deg,  
!               consistent with the number of pixels used to build the model -  
                alfa(i)=0.09
            enddo
       close(10)
       write(*,*) &
       '    - Read ', nel, &
       ' elements from ', ice_model       
!
!#--- Shape factors 
	do 80 i=1, nel 	
        if(mod(i,100)==0)write(*,*)&
	'     -  ', i, &
	' of ', nel 
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
	            pppp(j,i) = t(lj(j),i)*pplm(j)*cmplx(cosdd(arg),-sindd(arg))
		enddo		 
!		
80 	continue 
!
!
!#      ----------------------
!#====> An ice sheet like IJ05
!#      ----------------------
!
	ELSEIF(ICE_MODEL(1:4)=='ij05') THEN 

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
       WRITE(*,*) & 
       '    - Read ',nel, & 
       ' ice elements from ', ice_model
!
!#--- Degree-dependent factors	  
	do i=1,nel
	  call pleg(lmax+1, 90.-alfa(i), leg) 	  
!
	  t(0,i)=(1.-cosdd(alfa(i)))/2.
	  	do l=1, lmax
	     		t(l,i)= (-leg(l+1)+leg(l-1))/(2.*l+1.)/2.
	  	enddo 	    
        enddo
!
!#--- Shape factors 
	do i=1, nel 
		if(mod(i,100)==0)write(*,*)&
		'    -  ', & 
		i, ' of ', nel 
	        call harmo(lmax, longc(i), latic(i), armoy) 
			do j=1, jmax 
	            		pppp(j,i) = t(lj(j),i)*conjg(armoy(j))  
			enddo
	enddo
!
!
!#     -----------------------
!#====> An ice sheet like DISK
!#     -----------------------
!
	ELSEIF(ICE_MODEL(1:4)=='disk') THEN 
!
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
            write(*,*)&
	    '    - Read ',nel, &
	    '  elements from ', ice_model
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
	    do 95 i=1, nel 
		if(mod(i,100)==0)write(*,*)&
            '    - ', i, ' of ', nel 
	        call harmo(lmax, longc(i), latic(i), armoy) 
		do 95 j=1, jmax 
	            pppp(j,i) = t(lj(j),i)*conjg(armoy(j))  
95	    continue
!
!
!#     ----------------------
!#====> An ice sheet like CAP
!#     ----------------------
!
	ELSEIF(ICE_MODEL(1:4)=='icap') THEN 
!
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
            write(*,*)&
	    '    - Read ',nel, &
	    '  elements from ', ice_model
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
!#---- Shape factors  	
	    do 96 i=1, nel 
		if(mod(i,100)==0)write(*,*)&
            '    - ', i, ' of ', nel 
	        call harmo(lmax, longc(i), latic(i), armoy) 
		do 96 j=1, jmax 
	            pppp(j,i) = t(lj(j),i)*conjg(armoy(j))  
96	    continue	
!
	ENDIF
!
!#      --------------------------
!#====> End of available models...
!#      --------------------------
!
!
!
	Write(*,*)'    - Computed ', nel, ' shape factors'
!
	Write(*,*)'    - Writing the shape factors on "shicec.bin"'

 open (7,file=fileout,status='unknown',form='unformatted')
 write(7) pppp
 close(7) 

! open (7,file=fileout,status='unknown',form='unformatted')
! write(7) pppp
! close(7) 
! Write(*,*)'    - done!'
!
!
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
       END PROGRAM SHAPE_FACTORS 
!
!
!
