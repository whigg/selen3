!
! This is program "SHICE.F90" 
!
! Created by GS July 8 2008  -Intel version 2.6- 
! Updated GS July 24 for implementation of the "disk load"
! Modified August 2008 for implementation of ANU ice on v. 2.7
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
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
       COMPLEX*8 TTTT(JMAX,0:NN+1), IIII(JMAX,0:NN), PPPP(JMAX,NEL)    
       INTEGER, PARAMETER :: NH1=20, NH3=20, NH5=24 
       INTEGER, PARAMETER :: NHA=4,  NHI=18, NHD=20 
       INTEGER, PARAMETER :: NHU=26, NHC=20       
       INTEGER I, J, K, L, N, IJ, NH, I1, I2, ICR(0:NN) 
       REAL*4 AJ, H(NEL,0:NN+1), CR(0:NN)      
       CHARACTER*2 CJ
!
!
!#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!      - Input files:  
!      ice_model (from 'data.inc'), and:  
       CHARACTER*20, PARAMETER :: FILEIN='shicec.bin'      

!      - Output files:  
       CHARACTER*20, PARAMETER :: FILEOUT='shice.dat'
!      shtXX.dat, SH files at time XX (XX=00, 01, ..,N, NN+1) 
!#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!
!
!
!#---- Reading the header lines of the ice table
!	    	
            OPEN(10,FILE=ICE_MODEL,STATUS='UNKNOWN') 
!
	    	  if(ice_model(1:4)=='ice1')nh=nh1 
	    	  if(ice_model(1:4)=='ice3')nh=nh3
	    	  if(ice_model(1:4)=='ice5')nh=nh5
	    	  if(ice_model(1:4)=='alps')nh=nha
	    	  if(ice_model(1:4)=='ij05')nh=nhi
	    	  if(ice_model(1:4)=='disk')nh=nhd
		  if(ice_model(1:4)=='anu0')nh=nhu
	    	  if(ice_model(1:4)=='icap')nh=nhc
!		  
            	do j=1, nh 
	        	read(10,'(a20)') cj
	    	enddo	     
!	    
!#---- Reading the time-dependent thickness of the ice elements 
!	    	
       	    DO 10 I=1, NEL   
!	    
!#---- "ICE1" ----             	    
	    IF    (ICE_MODEL(1:4)=='ice1') THEN 
        	   Read (10,*) ij, ij,                (icr(k),k=nn,0,-1) 
		   cr(:)=icr(:)			
!
!#---- "ICE3" ----   
	    ELSEIF(ICE_MODEL(1:4)=='ice3') THEN 
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
       		    Read (10,113) ij, aj, ij, aj,     (cr(k),k=nn,0,-1) 			
		else
		    Read (10,114) ij, aj, ij, aj,     (cr(k),k=nn,0,-1) 
		Endif
!		        		
!#---- "ALPS" ----
	    ELSEIF(ICE_MODEL(1:4)=='alps') THEN 
       		   Read (10,*) ij, aj, aj, aj,        (cr(k),k=nn,0,-1) 
!
!#---- "IJ05" ----
	    ELSEIF(ICE_MODEL(1:4)=='ij05') THEN 
       		   Read (10,*)  ij, aj, aj, aj, aj,   (cr(k),k=nn,9,-1) 
        	   Read (10,*) 	                      (cr(k),k= 8,0,-1)   		    	                            ! 2nd row
!
!#---- "End of available ice models..." 
            ENDIF
!
!
!#---- Conversion to our default for ice thickness
!
	If(ice_model(1:4)/='anu0') then 
         	h(i,0)=cr(nn)              
       		do k=1,nn                   
                     h(i,k)=cr(nn-k)   
       		enddo 
                h(i,nn+1) = h(i,nn) 
	Endif
!
10          CONTINUE
!
	    close(10) 
!
            Write(*,*)&
	    '    - Read ',nel, '  elements from ', ice_model
            Write(*,*)&
	    '    - SH coefficients at time... '
!
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
!#----- Reading the shape factors 
!
 open(7,file=filein,status='unknown',form='unformatted') !,form='unformatted')
 read(7) pppp
 close(7) 
!
!
!#----- Computing the harmonic coefficients
!
 do 20 k=0,nn+1  
	open(3,file='junk.dat',status='unknown') 
	  if(k<=9) write(3,'(a1,i1)') '0',k  
	  if(k> 9) write(3,'(i2)')        k
 	close(3)
	open(3,file='junk.dat',status='unknown') 
	  read(3,'(a2)') cj 
 	close(3)  
	if(k==0.or.k==nn+1) write(*,'(a7,1x,a2,a3,i4)') '     - ', cj, & 
	           ' of', nn+1     
!				     
	open(1,file='sht'//cj//'.dat',status='unknown') 
		do j=1,jmax 
		tttt(j,k)=(0.,0.)
			do i=1,nel
			tttt(j,k)=tttt(j,k)+h(i,k)*pppp(j,i)
			enddo
		write(1,*) j, real(tttt(j,k)),aimag(tttt(j,k)) 
		enddo
!
	close(1) 
20 continue
!
!
!#----- Computing array "I"     
!
        open(7,file=fileout,status='unknown',form='formatted') 
        write(*,*) '    - Computing array "I"' 
        do 30 k=0,nn 
	      do 30 j=1,jmax 
	        iiii(j,k)=(0.,0.) 
		do 30 n=0, k
        	iiii(j,k)=iiii(j,k)+(tttt(j,n+1)-tttt(j,n))
30      continue
!
!#----- Storing array "I" 
!
 	write(7,*) iiii
 	close(7) 
!
 END PROGRAM SHICE 
!
!
!
