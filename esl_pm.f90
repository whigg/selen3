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
! CREATED ON MERCH 3 2012 by GS by a copy of "ESL.F90"
!      April 11, 2014. The BENOIT (BEN) files for the collaboration with LEGOS 
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
       INTEGER, PARAMETER :: NHZ=10     
       INTEGER, PARAMETER :: NHB=4       
!--- ... 
       INTEGER, PARAMETER :: ILARGE=1000     
       INTEGER I, J, K, L, N, IJ, NH, I1, I2, ICR(0:NN)
       REAL*8 TETAC(NEL), LONGC(NEL), LATIC(NEL), ALFA(NEL) 
       REAL*8 AJ, H(NEL,0:NN+2), CR(0:NN)      
       REAL*8 RATIO, VMIN, VMAX, VESL        
       REAL*8 VOL, ESL(0:NN+2) 
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
   if(ice_model(1:4)=='beno') nh=nhb
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
	        		
!#---- "GRD" ----
!
     IF(ICE_MODEL(1:3)=='grd') then 
     Read (10,*) ij, longc(i), latic(i), alfa(i), (h(i,k),k=0,nn+1) 
     h(i,nn+2)=h(i,nn+1) 
     endif    

     IF(ICE_MODEL(1:4)=='beno') then 
     Read (10,*) ij, longc(i), latic(i), alfa(i), (h(i,k),k=0,nn) 
     h(i,nn+1)=h(i,nn) 
     h(i,nn+2)=h(i,nn+1) 
     endif    
!
10 CONTINUE
!
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
	DO 20 K=0,NN+2
!
	   vol=0.
!	     
	   DO 30 I=1, NEL 
!		        		
!#---- "GRD" ----
!
        IF(ICE_MODEL(1:3)=='grd'.or.ICE_MODEL(1:4)=='beno') THEN 
!
		  vol = vol + & 
		  2.*pi*h(i,k)*(1.-cosdd(alfa(i))) 			
!
!#---- "End of available ice models..." 
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
	open(13,file='esl-pm.dat',     status='unknown')       
!    -thin line 
	open(14,file='esl-pm-thin.dat',status='unknown')       
!
! -- Before melting begins ice thickness is constant 	
	write(14,*) -DELTA*1000., esl(0) 
!	
! -- Step wise time history betwen time '0' and 'NN'
!
! Time is printed in units of YEARS!
!
	do k=0, NN+1 
 	    write(14,*) float(k)*DELTA*1000., esl(k) 
	    write(14,*) float(k)*DELTA*1000., esl(k+1)	    
 	    write(13,*) float(k)*DELTA*1000., esl(k) 
	    write(13,*) float(k)*DELTA*1000., esl(k+1)	    
	enddo
!
! -- After the end of melting ice thickness is constant 
	write(14,*) float(nn+2)*DELTA*1000., esl(NN+1)  
!
	close(13)	
	close(14)        
!

 END PROGRAM EQUIVALENT_SEALEVEL 
!
!
!
