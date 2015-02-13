!
! This is program "AM.F90" 
!
! Created GS & FC July 2009 -  "ALMA coupling" 
! Updated on April 13, 2010 by GS and FC - ALMA with degree 1
!             - Green Functions updated and re-organized
!             - A Geoid Green Function now explicitly given
!             - Clearing of all parts 
! Re-edited on May 2010 by GS (last time it was May 30)
! === Modified GS June 2010 - Free air gravity anomaly implementation 
! === Modified GS July 2010 - FA & SS anomalies revised 
! Feb 2012: Implementation of the numerical derivative "in the future" 
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
! ------------------------------------------------------------------------------
!  AM.F90 builds the 'BETA' arrays from the Love numbers computed by ALMA. It
!  has thus similar purposes of TB.F90, which computes the BETA arrays using
!  the TABOO (Normal Modes) Love numbers. - GS & FC July 2009 "ALMA coupling" -
! ------------------------------------------------------------------------------
!
! =============================================================================
!
! Computes the arrays Beta(l,k) and E(l), and Beta^u(l,k) and E^u(l), the only 
! rheology--dependent quantities involved in the solution of the SLE by the PS 
! method.  Here the program TABOO is used as a subroutine to compute the LDCs. 
!
! Input files:
!	- ALMA 'h', 'l', and 'k' tables of Love numbers  
!
! Output files:  
!	- ebs.dat    ! Sea level green function 
!	- ebu.dat    ! Vertical displacement Green function
!       - ebn.dat    ! Geoid elevation Green function
!       - ebv.dat    ! Hhorizontal displacement Green function
!       - ebfa.dat    ! Free air gravity anomaly Green function   <======
!       - ebss.dat    ! Solid surface gravity anomaly Green function   <======
!
! =============================================================================
!
!********** 
 PROGRAM AM
!**********
!
 IMPLICIT NONE        
 INCLUDE "data.inc"
!
 INTEGER, PARAMETER :: NH=18 
 INTEGER, PARAMETER :: LMIN=1, LMP1=LMIN+1, LSTP=1 
! ------------------------------
 REAL*4 ES(0:LMAX), & 
        EU(0:LMAX), &
        EN(0:LMAX), &  
        EV(0:LMAX), & 
	EFA(0:LMAX),&
	ESS(0:LMAX)		 
 REAL*4 BETAS(0:LMAX, 0:NN+1), &    ! Feb 2012: updated 
        BETAU(0:LMAX, 0:NN+1), & 
        BETAN(0:LMAX, 0:NN+1), & 
        BETAV(0:LMAX, 0:NN+1), & 
        BETAFA(0:LMAX,0:NN+1),&    
        BETASS(0:LMAX,0:NN+1) 
! ------------------------------
 REAL*4 HPW(LMIN:LMAX,0:NN+1), &   ! Feb 2012: updated 
        LPW(lMIN:LMAX,0:NN+1), & 
        KPW(LMIN:LMAX,0:NN+1)
 REAL*4  TIME(0:NN+1)              ! Feb 2012: updated 
 REAL*4  DEN
 INTEGER I, J, K, L, ND
 CHARACTER*10 CJUNK
!
!
!
! --- Opening the PW Heaviside Love numbers files 
 open(101,file='hpw.dat',status='unknown') 
 open(102,file='lpw.dat',status='unknown') 
 open(103,file='kpw.dat',status='unknown') 
!
! --- Reading the nh header lines 
 do i=1,3 
   do j=1,nh
         read(100+i,'(a10)')CJUNK
   enddo
 enddo
!
! --- Reading the PW Love Numbers ! Feb 2012: updated 
 do k=0,nn+1
	read(101,'(1024(e19.8))') TIME(K), (HPW(L,K),L= LMIN, LMAX, LSTP)
	read(102,'(1024(e19.8))') TIME(K), (LPW(L,K),L= LMIN, LMAX, LSTP)
	read(103,'(1024(e19.8))') TIME(K), (KPW(L,K),L= LMIN, LMAX, LSTP)
 enddo
!
 close(101)
 close(102)
 close(103) 
!
!
!###################################################
! Building the "E" arrays (Elastic Green functions)
!###################################################
!
   ES(:)=0.   
   EU(:)=0.
   EN(:)=0.
   EV(:)=0.
   EFA(:)=0. 
   ESS(:)=0.  
!
! --------------------------------------------------------------
! ====== Setting the ELASTIC Green Functions for Degree 1 ======
! --------------------------------------------------------------
!
 if(deg1==1) then 				
!
! ... if the degree 1 needs to be included in SLE modeling, the Green 
! functions at degree nd=1 are set by the corresponding Love numbers 
!
   nd=1
!
! ---- GSC or Elastic solution
	If(imode==1.or.imode==2.or.imode==5)  then 
   		EU(nd) =     hpw(nd,0)  
		EN(nd) =            0. 
   		EV(nd) =     lpw(nd,0) 	
		ES(nd) = EN(nd)-EU(nd)
!
	EFA(nd) = float(nd+2)	   
        ESS(nd) = float(nd)+2*hpw(nd,0)-float(nd+1)*kpw(nd,0) 
        ESS(nd) = -ESS(nd)
	Endif
!
! ---- Eustatic or solution 
	If(imode==3)  then 
		EU(nd) = 0.
		EN(nd) = 0. 
		EV(nd) = 0. 	
	        ES(nd) = 0.
!               EFA(nd) = ??????????????????????????????
!               ESS(nd) = ??????????????????????????????
	Endif 
!
! ---- Woodward solution 
	If(imode==4)  then 
		EU(nd) = 0.
		EN(nd) = 1. 
		EV(nd) = 0. 	
	        ES(nd) = 0.
!               EFA(nd) = ??????????????????????????????
!               ESS(nd) = ??????????????????????????????
	Endif 	
!
   den=(2.*float(nd)+1.) 	
!
   EU(nd) = EU(nd)/den 
   EN(nd) = EN(nd)/den 
   EV(nd) = EV(nd)/den 
   ES(nd) = ES(nd)/den 	
   EFA(nd) = EFA(nd)/den 	 	
   ESS(nd) = ESS(nd)/den 
!   
   Elseif(deg1==0) then 
!
! ... if the degree 1 needs NOT to be included, I set the Green functions at 
! degree nd=1 to zero, as it was the case in  previous versions of SELEN.
!
   nd=1 
!
   EU(nd) = 0.
   EN(nd) = 0. 
   EV(nd) = 0. 	
   ES(nd) = 0.
   EFA(nd) = 0.	 	
   ESS(nd) = 0.
!
 ENDIF 
!
!
! ----------------------------------------------------------------
! ====== Setting the ELASTIC Green Functions for Degree > 1 ======
! ----------------------------------------------------------------
!
 Do 10 l=LMP1, LMAX
!
	EU(l) =        hpw(l,0)  
	EN(l) =     1.+kpw(l,0) 
	EV(l) =        lpw(l,0)   
	ES(l) =   EN(l) - EU(l)	
!
      	EFA(l) =  float(l+2)              - float(l-1)*kpw(l,0)   	  
      	ESS(l) =  float(l)   + 2*hpw(l,0) - float(l+1)*kpw(l,0)  	  
        ESS(l) = -ESS(l)
!	
!         
! ---- Eustatic solution 
	if    (imode==3) then  
			EU(l)= 0.
			EN(l)= 0. 
			EV(l)= 0.
			ES(l)= EN(l)-EU(l)
!                       EFA(nd) = ??????????????????????????????
!                       ESS(nd) = ??????????????????????????????
!
! ---- Woodward solution	
	elseif(imode==4) then    
			EU(l)= 0.  
			EN(l)= 1.      
			EV(l)= 0.         
			ES(l)= EN(l)-EU(l)
!                       EFA(nd) = ??????????????????????????????
!                       ESS(nd) = ??????????????????????????????
	endif		
!
        den=	(2.*float(l)+1.) 
!
        EU(l) = EU(l)/den 
	EN(l) = EN(l)/den 
	EV(l) = EV(l)/den 
	ES(l) = ES(l)/den 
	EFA(l) = EFA(l)/den 
	ESS(l) = ESS(l)/den 
!
 10 continue
!
!
!############################################################
! Building the "BETA" arrays (Visco-Elastic Green functions)
!############################################################
!
!
 BETAS(:,:)=0. 
 BETAU(:,:)=0. 
 BETAN(:,:)=0.
 BETAV(:,:)=0.
 BETAFA(:,:)=0. 
 BETASS(:,:)=0. 

!
!
! --------------------------------------------------------------------
! ====== Setting the Visco-Elastic Green Functions for Degree 1 ======
! --------------------------------------------------------------------
!
! ... if the degree 1 needs to be included in SLE modeling, the Green 
! functions at degree nd=1 are set by the corresponding Love numbers 
!
 if(deg1==1) then 
	
	do k=0, NN+1  ! Feb 2012: updated 
		nd=1 
   		den=	(2.*float(nd)+1.) 
!
                BETAU(nd,k)= (hpw(nd,0)-hpw(nd,k)) /den
                BETAN(nd,k)=  0.0                         !(kpw(nd,0)-kpw(nd,k)) /den
                BETAV(nd,k)= (lpw(nd,0)-lpw(nd,k)) /den				    
		BETAS(nd,k)=  BETAN(nd,k)-BETAU(nd,k)              
!
                BETAFA(nd,k) = 0.0		
                BETASS(nd,k) = 2.*BETAU(nd,k) + FLOAT(nd+1)/DEN 		
	        BETASS(nd,k) = -BETASS(nd,k)
!
        enddo
!	
	elseif(deg1==0) then 
!
! ... if the degree 1 needs NOT to be included, I set the Green functions at 
! degree nd=1 to zero, as it was the case in  previous versions of SELEN.
!	
 do k=0, NN+1  ! Feb 2012: updated 	
     nd=1 
     BETAS(nd,k)=0.0     
     BETAU(nd,k)=0.0			     
     BETAN(nd,k)=0.0     
     BETAV(nd,k)=0.0			      
     BETAFA(nd,k)=0. 
     BETASS(nd,k)=0. 
 enddo
!	
 Endif
!
!
! ----------------------------------------------------------------------
! ====== Setting the Visco-Elastic Green Functions for Degree > 1 ======
! ----------------------------------------------------------------------
!
 do k=0, NN +1    ! Feb 2012: updated 
!
        do l=LMP1, LMAX 
!
   	den=	2.*float(l)+1. 
!
                BETAU(l,k)= (hpw(l,0)-hpw(l,k))/den
                BETAN(l,k)= (kpw(l,0)-kpw(l,k))/den
                BETAV(l,k)= (lpw(l,0)-lpw(l,k))/den				    
                BETAS(l,k)=  BETAN(l,k)-BETAU(l,k)
!
	        BETAFA(l,k)=  FLOAT(l-1)*BETAN(l,k)   						
	        BETASS(l,k)=  2.*BETAU(l,k) - FLOAT(l+1)*BETAN(l,k)   			
	        BETASS(l,k)= -BETASS(l,k)
!
        enddo
 enddo
!
!
!
!
!
! ----------------------------------- 
! ====== Reporting the results ======
! ----------------------------------- 
!
!
 open(1,file='ebs.dat',status='unknown')
 open(2,file='ebu.dat',status='unknown')
 open(3,file='ebn.dat',status='unknown')
 open(4,file='ebv.dat',status='unknown')
 open(7,file='ebfa.dat',status='unknown')
 open(8,file='ebss.dat',status='unknown')
!  
 do l=0, lmax
 	 write(1,*) l, ES(l),  (betas(l,k),  k=0, nn+1)  ! Feb 2012: updated 
 	 write(2,*) l, EU(l),  (betau(l,k),  k=0, nn+1)  ! Feb 2012: updated 
 	 write(3,*) l, EN(l),  (betan(l,k),  k=0, nn+1)	 ! Feb 2012: updated 	 
 	 write(4,*) l, EV(l),  (betav(l,k),  k=0, nn+1)  ! Feb 2012: updated 
 	 write(7,*) l, EFA(l), (betafa(l,k), k=0, nn+1)  ! Feb 2012: updated 
 	 write(8,*) l, ESS(l), (betass(l,k), k=0, nn+1)  ! Feb 2012: updated 
 enddo
!
 close(1)
 close(2)
 close(3)
 close(4) 
 close(7)
 close(8) 
!
 END PROGRAM AM 
!
!
!
!
