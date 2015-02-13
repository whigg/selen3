!
! This is program "ESL.F90"  
!
! Created by GS July 8 2008  -Intel version 2.6- 
! Updated GS July 24 for implementation of the "disk load"
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
! Reviewed GS on July 2010 for the implementation of the Greenland ice sheet 
! Reviewed GS July 2010 for the implementation of "small glaciers..."
! GS July 2010: Volumes and ESl are computed in DOUBLE PRECISION
! *** Revised GS July 2010 - g95 - Double precision implementation
! *** Revised GS & GR Nov 2010 - Implementation of MED1
! *** Revised GS & GG May 2011 - Implementation of the Antarctic model "ANTA"
! *** Revised GS & GG June 2011 - Implementation of the "TRANS" ice model 
! *** Revised GS July 2011 - Implementation of the "RECT" ice model 
!      February 2012: Introduced the "GRD" ice model type by GS
!      ANTA+GREE+GLAC models test for negative longitude === GS 01-07-2012
! Jul 1, 2013; The "FLOR" ice for coupling with her models 
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
!  ESL.F90 computes the Equivalent Sea Level (ESL) of the ice models in the 
!  library ./ICE-MODELS. This routine substitutes all the SH*_ESL.F90 program 
!  units that were written "ad-hoc" for individual ice models until v. 2.5 
! ---------------------------------------------------------------------------
!
!
!
!      INCLUDE "harmonics.f90"
       PROGRAM EQUIVALENT_SEALEVEL 
       IMPLICIT NONE        
       INCLUDE "data.inc"
!
!#---- General declarations 
!
!--- file headers 
       INTEGER, PARAMETER :: NH1=20, NH3=20, NH5=24 
       INTEGER, PARAMETER :: NHA=4,  NHI=18, NHD=20 
       INTEGER, PARAMETER :: NHU=26, NHC=20, NHG=10 
       INTEGER, PARAMETER :: NHM=20, NHR=10, NHT=20 
       INTEGER, PARAMETER :: NHY=20, NHZ=10, NHF=4     
!--- ... 
       REAL*8, PARAMETER :: AMP5 = 1.0, AMP1 = 5. !,  AMP5 = 0.70313
       INTEGER, PARAMETER :: ILARGE=1000     
       INTEGER I, J, K, L, N, IJ, NH, I1, I2, ICR(0:NN)
       REAL*8 TETAC(NEL), LONGC(NEL), LATIC(NEL), ALFA(NEL) 
       REAL*8 AJ, H(NEL,0:NN+1), CR(0:NN)      
       REAL*8 RATIO, VMIN, VMAX, VESL        
       REAL*8 VOL, ESL(0:NN+1) 
       REAL*8 COSDD, SINDD
       REAL*8 AMPT
       CHARACTER*2 CJ
!
!
!#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!      - Input files:  
!        Ice_model (from 'data.inc') 
         CHARACTER*20, PARAMETER :: FILEIN='shof.dat'
!
!      - Output files: 
!        esl.dat 
!	 esl-thin.dat
!        esl-tot.dat     
!#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!
!
!
!#---- Reading the header lines of the ice table
!	    	
  OPEN(10,FILE=ICE_MODEL,STATUS='UNKNOWN') 
!
   if(ice_model(1:3)=='grd') nh=nhz
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
   if(ice_model(1:4)=='gree'.or. & 
      ice_model(1:4)=='glac'.or. & 
      ice_model(1:4)=='anta')nh=nhg		  		  
   if(ice_model(1:4)=='tran')nh=nht
   if(ice_model(1:4)=='rett')nh=nhy
!		  
   do j=1, nh 
       read(10,'(a20)') cj
   enddo	     
!
!	    
!#---- Reading the time-dependent thickness of the ice elements 
!
!	    	
 DO 10 I=1, NEL   
!
!
!
!
!#---- "TRANS" ---- 
!
 IF(ICE_MODEL(1:4)=='tran') THEN 
  Read (10,*) ij, longc(i), latic(i), alfa(i), h(i,0) 
  Read (10,*) (h(i,k),k=1 ,10)
  Read (10,*) (h(i,k),k=11,20)
  Read (10,*) (h(i,k),k=21,30)
  Read (10,*) (h(i,k),k=31,40)
  Read (10,*) (h(i,k),k=41,50)
  h(i,nn+1)=h(i,nn)
!
!#---- ------ ---- 
!
!
!#---- "RECT" ---- NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW -----
!
 ELSEIF (ICE_MODEL(1:4)=='rett') THEN 
  Read (10,*) ij, longc(i), latic(i), AMPT, h(i,0) 
!
  if(longc(i).le.0.0) longc(i) = 360.+longc(i)
  tetac(i) = 90.-latic(i)
!
  Read (10,*) (h(i,k),k=1 ,10)
  Read (10,*) (h(i,k),k=11,20)
  Read (10,*) (h(i,k),k=21,30)
  Read (10,*) (h(i,k),k=31,40)
  Read (10,*) (h(i,k),k=41,50)
  h(i,nn+1)=h(i,nn)
!
!#---- ------ ---- NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW -----
!	    
!#---- "ICE1" ----
!             	    
 ELSEIF    (ICE_MODEL(1:4)=='ice1') THEN 
  Read (10,*) i1, i2, (icr(k),k=nn,0,-1) 
  cr(:)=icr(:)	
  tetac(i)=float(i1) 
  longc(i)=float(i2) 		
!
!
!#---- "ICE3" ----   
!
 ELSEIF(ICE_MODEL(1:4)=='ice3') THEN 
  Read (10,111) ij, tetac(i), longc(i), alfa(i), aj, (cr(k),k=nn,11,-1)
  Read (10,112)                                      (cr(k),k=10, 0,-1)		   		
  latic(i)=90.-tetac(i) 
!
!
!#---- "IMED" ----   
!
 ELSEIF(ICE_MODEL(1:4)=='imed') THEN 
  Read (10,111) ij, tetac(i), longc(i), alfa(i), aj, (cr(k),k=nn,11,-1)
  Read (10,112)                                      (cr(k),k=10, 0,-1)		   		
  latic(i)=90.-tetac(i) 
!
!
!#---- "DISK" ----   
!
 ELSEIF(ICE_MODEL(1:4)=='disk') THEN 
  Read (10,*) ij, tetac(i), longc(i), alfa(i), (cr(k),k=nn,11,-1)
  Read (10,*)                                  (cr(k),k=10, 0,-1)		   		
  latic(i)=90.-tetac(i) 
!
!
!#---- "ICE5" ----
!
 ELSEIF(ICE_MODEL(1:4)=='ice5') THEN 
  If (ice_model/='ice5g26.dat') &  
    Read (10,113) ij, longc(i), ij, latic(i), (cr(k),k=nn,0,-1) 			
  If (ice_model=='ice5g26.dat') &  
    Read (10,114) ij, longc(i), ij, latic(i), (cr(k),k=nn,0,-1) 
  tetac(i)=90.-latic(i)
!
!		        		
!#---- "GRD" ----
!
 ELSEIF(ICE_MODEL(1:3)=='grd') THEN 
    Read (10,*) ij, longc(i), latic(i), alfa(i), (cr(k),k=nn,0,-1) 
!   alfa(i)=0.09
!
!
!		        		
!#---- "ALPS" ----
!
 ELSEIF(ICE_MODEL(1:4)=='alps') THEN 
    Read (10,*) ij, longc(i), latic(i), alfa(i), (cr(k),k=nn,0,-1) 
!alfa(i)=0.09
!
!		        		
!#---- "FLOR" ----
!
 ELSEIF(ICE_MODEL(1:4)=='flor') THEN 
    Read (10,*) ij, longc(i), latic(i), alfa(i), (cr(k),k=nn,0,-1) 
!
!		        		
!#---- "ANTA" ----
!
 ELSEIF(ICE_MODEL(1:4)=='anta') THEN 
    Read (10,*) ij, longc(i), latic(i), alfa(i), (cr(k),k=nn,0,-1)
    if(longc(i).le.0.0) longc(i)=longc(i)+360.0
!do k=nn, 0, -1
!cr(k)=cr(k)*SCALING_ICE
!enddo
!
!		        		
!#---- "GREEN" ----
!
 ELSEIF(ICE_MODEL(1:4)=='gree') THEN 
    Read (10,*) ij, longc(i), latic(i), alfa(i), (cr(k),k=nn,0,-1) 
    if(longc(i).le.0.0) longc(i)=longc(i)+360.0
!do k=nn, 0, -1
!cr(k)=cr(k)*SCALING_ICE
!enddo
!
!		        		
!#---- "GLAC" ----
!
 ELSEIF(ICE_MODEL(1:4)=='glac') THEN 
    Read (10,*) ij, longc(i), latic(i), alfa(i), (cr(k),k=nn,0,-1) 
    if(longc(i).le.0.0) longc(i)=longc(i)+360.0
!do k=nn, 0, -1
!cr(k)=cr(k)*SCALING_ICE
!enddo
!
!		        		
!#---- "ICAP" ----
!
 ELSEIF(ICE_MODEL(1:4)=='icap') THEN 
  Read (10,*) ij, tetac(i), longc(i), alfa(i), (cr(k),k=nn,11,-1)
  Read (10,*)                                  (cr(k),k=10, 0,-1)		   		
  latic(i)=90.-tetac(i) 
!
!
!#---- "IJ05" ----
!
 ELSEIF(ICE_MODEL(1:4)=='ij05') THEN 
  Read (10,*) ij, longc(i), latic(i), alfa(i), aj, (cr(k),k=nn,9,-1) 
  Read (10,*) 	                                   (cr(k),k= 8,0,-1)   		    	                      
!
!
!#---- "ANU05" ----
!
 ELSEIF(ICE_MODEL(1:4)=='anu0') THEN 
!
    Read (10,115) ij, longc(i), latic(i), alfa(i), (h(i,k),k=0,nn+1) 
    tetac(i)=90.-latic(i) 

!
!#---- "End of available ice models..." 
!
 ENDIF
!
!
!#---- Conversion to our default for ice thickness
!
 If(ICE_MODEL(1:4)/='anu0'.and.&
    ICE_MODEL(1:4)/='tran'.and.&
    ICE_MODEL(1:4)/='rett')       then 
!
 h(i,0)=cr(nn)              
      do k=1,nn                   
         h(i,k)=cr(nn-k)   
      enddo 
 h(i,nn+1) = h(i,nn)
 Endif 
!
!
10 CONTINUE
!
   close(10) 
!
!
!  Write(*,*)'    - esl: read ',nel, '  elements from ', ice_model
!
!
!#---- Ratio: Area(oceans)/Area(surface of the Earth) 
!      [i. e., the degree 0 SH coefficient of the OF]
!
   open(7,file=filein,status='old')
   read(7,*) ij, ratio
   close(7)
!
!
!
!#---- Equivalent Sea level vs. time
!
	DO 20 K=0,NN+1
!
	   vol=0.
!	     
	   DO 30 I=1, NEL 
!	    
!#---- "RECT" ----
!
! **************************************************** NEW NEW NEW 
!
           IF    (ICE_MODEL(1:4)=='rett') THEN 
!		  
		  vol = vol + & 
		  h(i,k)*(ampt*pi/180.)*&
                  (cosdd(tetac(i)-ampt/2.)-cosdd(tetac(i)+ampt/2.))
!		   
! **************************************************** NEW NEW NEW 
!	    
!#---- "ICE1" ----
!
           ELSEIF(ICE_MODEL(1:4)=='ice1') THEN 
!		  
		  vol = vol + & 
		  h(i,k)*(amp1*pi/180.)*&
                  (cosdd(tetac(i)-amp1/2.)-cosdd(tetac(i)+amp1/2.)) 			
!
!#---- "TRANS" ----   
!
           ELSEIF(ICE_MODEL(1:4)=='tran') THEN 
!	   
		  vol = vol + & 
		  2.*pi*h(i,k)*(1.-cosdd(alfa(i))) 			
!
!
!#---- "ICE3" ----   
!
           ELSEIF(ICE_MODEL(1:4)=='ice3') THEN 
!	   
		  vol = vol + & 
		  2.*pi*h(i,k)*(1.-cosdd(alfa(i))) 			
!
!#---- "IMED" ----   
!
           ELSEIF(ICE_MODEL(1:4)=='imed') THEN 
!	   
		  vol = vol + & 
		  2.*pi*h(i,k)*(1.-cosdd(alfa(i))) 
!
!#---- "DISK" ----   
!
           ELSEIF(ICE_MODEL(1:4)=='disk') THEN 
!	   
		  vol = vol + & 
		  2.*pi*h(i,k)*(1.-cosdd(alfa(i))) 
!
!
!#---- "ICAP" ----   
!
           ELSEIF(ICE_MODEL(1:4)=='icap') THEN 
!	   
		  vol = vol + & 
		  (4./3.)*pi*h(i,k)*(1.-cosdd(alfa(i)))		 
!
!#---- "ICE5" ----   
! 
           ELSEIF(ICE_MODEL(1:4)=='ice5') THEN 
!
		  vol = vol + & 
		  h(i,k)*(amp5*pi/180.)*&
	          (cosdd(tetac(i)-amp5/2.)-cosdd(tetac(i)+amp5/2.)) 
!		        		
!#---- "GRD" ----
!
           ELSEIF(ICE_MODEL(1:3)=='grd') THEN 
!
		  vol = vol + & 
		  2.*pi*h(i,k)*(1.-cosdd(alfa(i))) 			
!
!		        		
!#---- "ALPS" ----
!
           ELSEIF(ICE_MODEL(1:4)=='alps') THEN 
!
		  vol = vol + & 
		  2.*pi*h(i,k)*(1.-cosdd(alfa(i))) 			
!		        		
!#---- "FLOR" ----
!
           ELSEIF(ICE_MODEL(1:4)=='flor') THEN 
!
		  vol = vol + & 
		  2.*pi*h(i,k)*(1.-cosdd(alfa(i))) 			

!		        		
!#---- "ANTA" ----
!
           ELSEIF(ICE_MODEL(1:4)=='anta') THEN 
!
		  vol = vol + & 
		  2.*pi*h(i,k)*(1.-cosdd(alfa(i))) 
!		        		
!		        		
!#---- "GREEN" ----
!
           ELSEIF(ICE_MODEL(1:4)=='gree') THEN 
!
		  vol = vol + & 
		  2.*pi*h(i,k)*(1.-cosdd(alfa(i))) 
!		        		
!#---- "GLAC" ----
!
           ELSEIF(ICE_MODEL(1:4)=='glac') THEN 
!
		  vol = vol + & 
		  2.*pi*h(i,k)*(1.-cosdd(alfa(i))) 
!
!#---- "IJ05" ----
!
           ELSEIF(ICE_MODEL(1:4)=='ij05') THEN 
!
		  vol = vol + & 
		  2.*pi*h(i,k)*(1.-cosdd(alfa(i))) 			
!
!#---- "ANU05" ----
!
           ELSEIF(ICE_MODEL(1:4)=='anu0') THEN 
!
		  vol = vol + & 
		  2.*pi*h(i,k)*(1.-cosdd(alfa(i))) 
!
!#---- "End of available ice models..." 
!
           ENDIF
!
 30        CONTINUE
!
        	ESL(K)=  (RHOI/RHOW)*VOL/RATIO/4./PI		
!
 20     CONTINUE
!
!
!
!
! -- data for a GMT plot
!    -thick line 
	open(13,file='esl.dat',     status='unknown')       
!    -thin line 
	open(14,file='esl-thin.dat',status='unknown')       
!
! -- Before melting begins ice thickness is constant 	
	write(14,*) "-5 ", esl(0) 
!	
! -- Step wise time history betwen time '0' and 'NN'
	do k=0, NN 
 	    write(14,*) k, esl(k) 
	    write(14,*) k, esl(k+1)	    
 	    write(13,*) k, esl(k) 
	    write(13,*) k, esl(k+1)	    
	enddo
!
! -- After the end of melting ice thickness is constant 
	write(14,*) NN+2, esl(NN+1)  
!
	close(13)	
	close(14)        
!
!
! --- Total ESL = max(ESL)-MIN(ESL)  
	vmax=-999999.
	vmin=+999999. 
	open(4,file='esl.dat',status='old') 
	do i=1, ilarge
		read(4,*,end=1) k, vesl 
		if(vesl<=vmin)vmin=vesl 
		if(vesl>=vmax)vmax=vesl 
	enddo
 1      close(4)
!
! --- Total ESL = max(ESL)-MIN(ESL)  
	if(ICE_MODEL(1:4)=='gree'.or. &
           ICE_MODEL(1:4)=='glac'.or. &
	   ICE_MODEL(1:4)=='anta') THEN 	
           Write(*,*) "    - TOTAL ESL (mm): ", (vmax-vmin)*1000.          
	Else
           Write(*,*) "    - TOTAL ESL (m) : ",  vmax-vmin          
        Endif
! 
! --- This string is later used by "esplot.gmt"
	open(4,file='esl-tot.dat',status='unknown')  
	IF(ICE_MODEL(1:4)=='alps') &
	write(4,'(a29,f8.2,a2)') "-1 1.8 16 0 1 BL Total ESL = ", vmax-vmin, " m"
	IF(ICE_MODEL(1:4)/='alps') &
	write(4,'(a29,f8.2,a2)') "-1 200 16 0 1 BL Total ESL = ", vmax-vmin, " m"
	close(4)  
!
!
!
!
!#---- Formats for ICE3G (formerly provided by JXM) 
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
!#----- 
!
 END PROGRAM EQUIVALENT_SEALEVEL 
!
!
!
