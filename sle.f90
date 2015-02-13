!
! This is program "SLE.F90" 
!
! GS 04-11-2008 "Intel port..."
! Modified by GS september 2008 for horizontal movements 
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
! Also reviewed in April 2010 for the implementation of degree 1 
! (remember the mode coupling issue...) 
! Revised May 26 GS for new ice SH  
! Updated GS June 2010 for Free Air & SS Gravity implementation 
! Updated GS July 2010 for the imode=5 option (see Gabriella) 
! Updated GS July 2010 for the Load implementation 
! Updated GS March 2011 for the "PHI" gravity potential implementation
! Updated DM November 2011 - Parallel execution
! Feb 2012: Implementation of the numerical derivative "in the future" 
! April 22, 2012: The eustatic vector dumped for later use (Urbino hospital)
! Aug 5, 2012; Equivalent water height (understanding the Chambers et al. 2010 results)
! Dec 21, 2012, Re-touched after the rotational counterpart has been updated
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
! SELEN is distributed in the /hope that it will be useful, but WITHOUT ANY 
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
! details. 
! 
! You should have received a copy of the GNU General Public License along 
! with SELEN.  If not, see <http://www.gnu.org/licenses/>.
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! ============================================
! Solves the SLE by the pseudo-spectral method 
! ============================================
!
! Input files:
!       - px-table.dat 
!	- shof.dat 
!	- shice.dat
!	- eb.dat
!	- ebu.dat
!       - ebv.dat 
!       - ebfa.dat  <<< new... 
! 	- sh.bin
!       - anchor.tmp  <<< new
!
! Output files:
!	- shu.bin 
!	- shn.bin 
!	- shs.bin 
!	- shg.bin     <<< new 2011  
!       - shz.bin 
!       - shv.bin 
!       - shfa.bin    <<< new... 
!       - shss.bin    <<< new... 
!       - shload*.bin <<< new 
!
!INCLUDE "harmonics.f90"
 PROGRAM SLE 
 IMPLICIT NONE 
 INCLUDE "data.inc"			    
!
 CHARACTER*12 HEADER 
 INTEGER :: NANCH 
 INTEGER I, J, K, L, P, IJ, IS, ISTEP, LJ, MJ, LI, DOM, IND
 INTEGER J_INDEX
 INTEGER, ALLOCATABLE :: LL(:), MM(:), DM(:), ANC(:), WET(:)       
!
 REAL*4 RHOE, RHOI_O_RHOE_X3, RHOW_O_RHOE_X3, RHOI_O_RHOW 
 REAL*8, ALLOCATABLE :: ALF(:,:), LONP(:), LATP(:), X(:,:)         
 REAL*8 RESH, IMSH
 REAL*8, ALLOCATABLE :: AAVV(:), BBVV(:)                           
!
 COMPLEX*16, ALLOCATABLE :: LONG_TABLE(:,:), OC(:)                 
 COMPLEX*16, ALLOCATABLE :: IIII(:,:), ZE(:,:)      
 COMPLEX*16, ALLOCATABLE :: AAAA(:,:), AAAA_MOD(:,:)               
 COMPLEX*16, ALLOCATABLE :: BBBB(:,:), BBBB_MOD(:,:)               
 COMPLEX*16, ALLOCATABLE :: HHHH(:,:), KKKK(:,:) 
 COMPLEX*16, ALLOCATABLE :: G_GRACE(:,:)                  
! 
!--------------------------------------------------
! Convolved Green functions (viscous & elastic)  (New as of Feb 2012: it was "nn") 
! 
 REAL*4, ALLOCATABLE :: BETAS(:,:),  ES(:), &  ! Sea level change  
                     BETAN(:,:),  EN(:),    &  ! Geoid undulations 
                     BETAU(:,:),  EU(:),    &  ! Vertical displacemnt
                     BETAV(:,:),  EV(:),    &  ! Horizontal displacement  
                     BETAFA(:,:), EFA(:),   &  ! Free air (FA) gravity anomaly 
                     BETASS(:,:), ESS(:)       ! Solid-surface (SS) gravity anomaly                  
!
 REAL*8, ALLOCATABLE :: KPRIME(:)                
!
 REAL*8 ERR(0:NN+1,1:SMAX)
 REAL*8 C_CONST(0:NN+1)			    
!
!--------------------------------------------------
! Output arrays & filenames (Updated Feb 2012: it was "nn")
! 
 CHARACTER*22 Out_Filename  
! - "Eustatic" Sea level change -------------------  ! New as of Apr 2012
 complex*16, ALLOCATABLE :: SE(:,:)   
 CHARACTER*22, PARAMETER :: FILEE='she.bin'   
!
! - Sea level change ------------------------------
 complex*16, allocatable :: S(:,:)             
 CHARACTER*22, PARAMETER :: FILES='shs.bin'   
!
! - Geoid undulations -----------------------------
 complex*16, allocatable :: N(:,:)   
 CHARACTER*22, PARAMETER :: FILEN='shn.bin'   
!
! - Vertical displacemnt --------------------------
 complex*16, allocatable :: U(:,:)   
 CHARACTER*22, PARAMETER :: FILEU='shu.bin'   
!
! - Normalized Gravity potential ------------------  
 complex*16, allocatable :: G(:,:)   
 CHARACTER*22, PARAMETER :: FILEG='shg.bin'   
!
! - Equivalent water heoght -----------------------  
 complex*16, allocatable :: W(:,:)   
 CHARACTER*22, PARAMETER :: FILEW='shw.bin'   
!
! - Horizontal displacement -----------------------
 complex*16, allocatable :: V(:,:)   
 CHARACTER*22, PARAMETER :: FILEV='shv.bin'  
!
! - Free air (FA) gravity anomaly -----------------
 complex*16, allocatable :: FA(:,:)   
 CHARACTER*22, PARAMETER :: FILEFA='shfa.bin'  
!
! - Solid-surface (SS) gravity anomaly ------------
 complex*16, allocatable :: SS(:,:)   
 CHARACTER*22, PARAMETER :: FILESS='shss.bin' 
!
! - Ice component of load  ------------------------
 complex*16, allocatable :: LO_I(:,:)
 CHARACTER*22, PARAMETER :: FILELI='shload_ice.bin'
!
! - Ocean component of load -----------------------  
 complex*16, allocatable :: LO_O(:,:)
 CHARACTER*22, PARAMETER :: FILELO='shload_oce.bin' 
!
! - Total load  -----------------------------------
 complex*16, allocatable :: LO_T(:,:)
 CHARACTER*22, PARAMETER :: FILELT='shload.bin'  
!
! - Sea level times the ocean function ------------  
 complex*16, allocatable :: Z(:,:,:) 
 CHARACTER*22, PARAMETER :: FILEZ='shz.bin'                       
!--------------------------------------------------
!
!
!
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! --- Getting info about the number of anchor pixels 
!
  open(1,file='anchor.tmp',status='old')
  read(1,*) nanch
  close(1)
!
! --- Allocate memory space
!
! New as of Feb 2012: 'nn' -> 'nn+1' 
 ALLOCATE( LL(JMAX), MM(JMAX), DM(JMAX), ANC(NP), WET(NP) )
 ALLOCATE( ALF(JMAX,NANCH), LONP(NP), LATP(NP), X(NP,0:NN+1) )
 ALLOCATE( AAVV(0:NN+1), BBVV(0:NN+1) ) 
 ALLOCATE( LONG_TABLE(0:LMAX,NP), OC(JMAX) )
 ALLOCATE( IIII(JMAX,0:NN+1), ZE(JMAX,0:NN+1) ) ! New as of Feb 2012
 ALLOCATE( AAAA(JMAX,0:NN+1), AAAA_MOD(JMAX,0:NN+1) )
 ALLOCATE( BBBB(JMAX,0:NN+1), BBBB_MOD(JMAX,0:NN+1) )
 ALLOCATE( HHHH(JMAX,0:NN+1), KKKK(JMAX,0:NN+1) )
 ALLOCATE( BETAS (0:LMAX,0:NN+1), ES (0:LMAX) ) ! New as of Feb 2012
 ALLOCATE( BETAN (0:LMAX,0:NN+1), EN (0:LMAX) )
 ALLOCATE( BETAU (0:LMAX,0:NN+1), EU (0:LMAX) )
 ALLOCATE( BETAV (0:LMAX,0:NN+1), EV (0:LMAX) )     
 ALLOCATE( BETAFA(0:LMAX,0:NN+1), EFA(0:LMAX) )     
 ALLOCATE( BETASS(0:LMAX,0:NN+1), ESS(0:LMAX) )     
 ALLOCATE( S(JMAX,0:NN+1), SE(JMAX,0:NN+1) )
 ALLOCATE( N(JMAX,0:NN+1), U(JMAX,0:NN+1), G(JMAX,0:NN+1) )
 ALLOCATE( V(JMAX,0:NN+1), W(JMAX,0:NN+1), FA(JMAX,0:NN+1) )
 ALLOCATE( SS(JMAX,0:NN+1), LO_I(JMAX,0:NN+1), LO_O(JMAX,0:NN+1) )
 ALLOCATE( LO_T(JMAX,0:NN+1), Z(JMAX,0:NN+1,0:SMAX) )
 ALLOCATE( KPRIME(1:LMAX) )
 ALLOCATE( G_GRACE(JMAX,0:NN+1) )
!
!
! ********************************************************
! ********************************************************
! ********************************************************
!
! --- Extracting the average Earth density from 
!     the TABOO or ALMA log files... GS July 09
!
 IF(CDE.GE. 0) IND=1 
 IF(CDE.EQ.-1) IND=2  
 IF(CDE.GE. 0.or.CDE==-1) CALL AVERAGE_EARTH_DENSITY(IND, RHOE)
! 
 IF(CDE==-2) RHOE=5511.57
! 
 RHOI_O_RHOE_X3 = 3.*RHOI/RHOE 
 RHOW_O_RHOE_X3 = 3.*RHOW/RHOE 
 RHOI_O_RHOW    =    RHOI/RHOW 
!
!
!
! --- Pre-computing 'l' and 'm' corresponding to degree 'J'
!
	do j=1, jmax 
		mm(j)=mj(j) 
		ll(j)=lj(j)
		dm(j)=2-dom(j)
	enddo	
!
! --- Reading the ALFs table from <<sh.bin>>
! 
 	Write(*,*) '    - Reading the ALFs from file sh.bin'
 	open(3,file='sh.bin',status='unknown',form='unformatted') 
		read(3)ALF
		read(3)LONG_TABLE
 	Close(3) 
!
! --- Examining the pixels table & extracting information
!
	Open(1,file='px-table.dat',status='unknown') 
	Do i=1, 4 
		Read(1,'(a12)')header
	Enddo
	Do i=1, np 
		Read (1,*) lonp(i), latp(i), anc(i), k, wet(i) 	
	Enddo
	close(1)
!
!
! --- Reading the SH OF coefficients from shof.dat 
!
 	Write(*,*) '    - Reading the SH OF coeff. from shof.dat'
 	open(3,file='shof.dat',status='unknown')
 	do j=1, jmax   
		read(3,*) k, resh, imsh 
                oc(j)=cmplx(resh, imsh)	
 	enddo
 	close(3)
!
!
! --- Reading the <<I>> array 
 		Write(*,*) "    - Reading array 'I'"
		open(1,file='shice.dat',status='unknown') 
		read(1,*) IIII
		close(1)
!
!
! --- Reading the E and BETA arrays 
 	Write(*,*) "    - Reading arrays 'E' and 'Beta'"
	open(1,file='ebs.dat', status='unknown')    
	open(2,file='ebu.dat',status='unknown')    
	open(3,file='ebn.dat',status='unknown')    
	open(4,file='ebv.dat',status='unknown')    
	open(7,file='ebfa.dat',status='unknown')    
	open(8,file='ebss.dat',status='unknown')    
!
	do li=0, lmax
		read(1,*) l, Es(l),  (betas(l,k),  k=0,nn+1) ! New as of Feb 2012: it was "nn" 
		read(2,*) l, Eu(l),  (betau(l,k),  k=0,nn+1) ! New as of Feb 2012: it was "nn" 
		read(3,*) l, En(l),  (betan(l,k),  k=0,nn+1) ! New as of Feb 2012: it was "nn" 
		read(4,*) l, Ev(l),  (betav(l,k),  k=0,nn+1) ! New as of Feb 2012: it was "nn" 
		read(7,*) l, Efa(l), (betafa(l,k), k=0,nn+1) ! New as of Feb 2012: it was "nn" 
		read(8,*) l, Ess(l), (betass(l,k), k=0,nn+1) ! New as of Feb 2012: it was "nn" 
	enddo
	close(8); close(7); close(4); close(3) ; close(2) ; close(1) 
!
!
! --- Computing the eustatic Z array...   ! New as of Feb 2012: it was "nn"
	ZE(:,:) = (0D0,0D0)  		
	DO K=0,NN+1	
		ZE(:,K) = - RHOI_O_RHOW*(IIII(1,K)/OC(1))*OC(:)
	enddo
!
! --- Computing the eustatic S array...   ! New as of Feb 2012: "SE" was dimensioned "nn"
	SE(:,:) = (0D0,0D0) 
	SE(1,:) = - RHOI_O_RHOW*(IIII(1,:)/OC(1)) 
!
!
! --- No ice load for imode==5 
!
      IF(IMODE==5) IIII(:,:)=0.0 
!
!
! --- Computing the A array... ! New as of Feb 2012: it was "nn"
	AAAA(:,:)=(0D0,0D0) 
	DO J=1, JMAX 
 	    DO K=0,NN+1
            AAAA(J,K) = ES(LL(J))*IIII(J,K)        
 	      DO P=0, K
              IF(P==0) AAAA(J,K) = AAAA(J,K)-(IIII(J,P)            )*BETAS(LL(J),K-P)
    	      IF(P/=0) AAAA(J,K) = AAAA(J,K)-(IIII(J,P)-IIII(J,P-1))*BETAS(LL(J),K-P)  
 	      ENDDO 
	      ENDDO
	ENDDO
!
! @BUG found and corrected by GS March 07 2013!!!
!
        AAAA(:,:)= RHOI_O_RHOE_X3*AAAA(:,:)
!
!
! --- Computing the ocean average of A ! New as of Feb 2012: it was "nn" 
	DO K=0, NN+1 
		AAVV(K)=0.
		DO J=1, JMAX 
		AAVV(K) = AAVV(K) + & 
		          DM(J)*REAL(OC(J)*CONJG(AAAA(J,K)))/OC(1)  
		ENDDO
	ENDDO 
!
!
! --- Computing the modified array A 
	AAAA_MOD(:,:) = AAAA(:,:)
	AAAA_MOD(1,:) = AAAA(1,:)-AAVV(:) 
!
!
! --- Computing the R-array... ! New as of Feb 2012: it was "nn"
        X(:,:) = 0D0
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J,K) &
!$OMP     SHARED(X,ALF,ANC,DM,AAAA_MOD,LONG_TABLE,MM) &
!$OMP     SCHEDULE(GUIDED)
	DO I=1, NP    
	DO K=0, NN+1
	        DO J=1, JMAX  
        	X(I,K) = X(I,K) + & 
			 ALF(J,ANC(I))*DM(J)*REAL(AAAA_MOD(J,K)*LONG_TABLE(MM(J),I)) 
		ENDDO   
	    ENDDO
	ENDDO 
!$OMP END PARALLEL DO
!      
	HHHH(:,:) = (0D0,0D0)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J,K) &
!$OMP       SHARED(WET,HHHH,X,ALF,ANC,LONG_TABLE,MM) &
!$OMP       SCHEDULE(GUIDED)
	DO J=1, JMAX 
	    DO I=1, NP 
	       IF(WET(I)==1) THEN
	          DO K=0,NN+1 ! NEW AS OF FEB 2012: IT WAS "NN"
		      HHHH(J,K) = HHHH(J,K) + X(I,K)*ALF(J,ANC(I))*CONJG(LONG_TABLE(MM(J),I))  
		  END DO
	       END IF
	    ENDDO
	ENDDO
!$OMP END PARALLEL DO
!	
! === Updating the H-array by the eustatic Z-term ===================
	HHHH(:,:)=HHHH(:,:)/FLOAT(NP) + ZE(:,:)
!
!
! --- Initializing the Z and S arrays 
	Z(:,:,0) = ZE(:,:)
	S(:,:)   = SE(:,:) 
!
!
! --- No recursion for the "Explicit approach"
	if(imode==6.or.imode==7.or.imode==3.or.SMAX==0) goto 2000     
!
!
! ----------------------- ----------------------- ----------------------- -----------------------
! ---    Recursion    --- ---	 Recursion    --- ---	 Recursion    --- ---	 Recursion    ---
! ----------------------- ----------------------- ----------------------- -----------------------
!
  Write(*,*) "    - Starting the recursion"
!
!
!
! ----------------------------
 	DO 1000 IS = 1, SMAX     
! ----------------------------
!
        write(*,*) " ======================= "
!
        write(*,'(a12,i2,a3,i2)') '     - step ', is, ' of', SMAX
!
!
! === COMPUTING THE 'B' ARRAY...  ===================================  RECURSION 
!
        BBBB(:,:)=(0D0,0D0) 
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(J,K,P) &
!$OMP     SHARED(BBBB,ES,LL,Z,IS,BETAS) &
!$OMP     SCHEDULE(GUIDED)
        DO J=1, JMAX 
	    DO K=0,NN+1
	    BBBB(J,K) = ES(LL(J))*Z(J,K,IS-1)        
	        DO P=0, K
	        IF(P==0) BBBB(J,K) = BBBB(J,K)-(Z(J,P,IS-1)-0.           )*BETAS(LL(J),K-P)
	        IF(P/=0) BBBB(J,K) = BBBB(J,K)-(Z(J,P,IS-1)-Z(J,P-1,IS-1))*BETAS(LL(J),K-P)		 
	        ENDDO 
        ENDDO
	ENDDO
!$OMP END PARALLEL DO
        BBBB(:,:) = RHOW_O_RHOE_X3*BBBB(:,:)
!
!	
! === COMPUTING THE OCEAN-AVERAGE OF ARRAY B ======================== RECURSION 
	DO K=0, NN+1 
		BBVV(K)=0D0
		DO J=1, JMAX 
		BBVV(K) = BBVV(K) + & 
		          DM(J)*REAL(OC(J)*CONJG(BBBB(J,K)))/OC(1)
		ENDDO
	ENDDO 
!
! === COMPUTING MODIFIED 'B' ARRAY ================================== RECURSION
	BBBB_MOD(:,:)=BBBB(:,:)
	BBBB_MOD(1,:)=BBBB(1,:)-BBVV(:) 
!
!
! === Computing array K...  ========================================= Recursion
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,K,J) &
!$OMP    SHARED(X,ALF,ANC,DM,BBBB_MOD,LONG_TABLE,MM) SCHEDULE(GUIDED)
	DO I=1, NP 	
 	    DO K=0, NN+1
              X(I,K)=0D0
	        DO J=1, JMAX  
        	    X(I,K) = X(I,K) + & 
			 ALF(J,ANC(I))*DM(J)*REAL(BBBB_MOD(J,K)*LONG_TABLE(MM(J),I)) 
		ENDDO   
	    ENDDO
	ENDDO 
!$OMP END PARALLEL DO
!
	KKKK(:,:) = (0D0,0D0)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J,K) &
!$OMP       SHARED(WET,KKKK,X,ALF,ANC,LONG_TABLE,MM) &
!$OMP       SCHEDULE(GUIDED)
	DO J=1, JMAX 
	    DO I=1, NP 
	        IF(WET(I)==1) THEN
                   DO K=0,NN+1          
                       KKKK(J,K) = KKKK(J,K) + X(I,K)*ALF(J,ANC(I))*CONJG(LONG_TABLE(MM(J),I))  
	           ENDDO
                ENDIF
	    ENDDO
	ENDDO
!$OMP END PARALLEL DO 
	KKKK(:,:)=KKKK(:,:)/FLOAT(NP) 
!
!
!
! --- Solving for array 'Z' st step "is" =========================== Recursion
!
	Z(:,:,is) = HHHH(:,:) + KKKK(:,:) 
!
!
! --- Solving for array 'S' ======================================== Recursion
!
	S(:,:)    = AAAA_MOD(:,:) + SE(:,:) + BBBB_MOD(:,:)

!
! ------------------------------
! ---    End of recursion    ---
! ------------------------------
!
1000 CONTINUE      
!
!
!
2000 CONTINUE       
!
!
! --- Eustatic solution: U=0, S=N, V=0  
        IF(IMODE==3.OR.SMAX==0) THEN 
				U(:,:) = (0D0,0D0)
				S(:,:) = SE(:,:)
				N(:,:) = S(:,:)
				V(:,:) = (0D0,0D0)
				GOTO 3000
				ENDIF
!
!
! ---------------------------------------------------------------------------------
! ### In PROGRESS ### In PROGRESS ### In PROGRESS ### In PROGRESS ### In PROGRESS
! ---------------------------------------------------------------------------------
!      
       Open(188,file='Error-iterations.dat',status='unknown') 
! 
       DO 11200 ISTEP=1, SMAX
		write(*,*) " Evaluating the convergence error, step: ", ISTEP 
!
              DO K=0, NN+1
       		ERR(K,ISTEP)=0D0
!
       		DO J=1,JMAX 
!     
       		ERR(K,ISTEP) = ERR(K,ISTEP) +  (ABS(Z(J,K,ISTEP+0)) - ABS(Z(J,K,ISTEP-1)))/ & 
                               ABS(Z(J,K,ISTEP-1))           
                ENDDO
             ENDDO
!
11200  CONTINUE
!
       DO K=0, NN+1    
       		WRITE (188,*) (ERR(K,ISTEP), ISTEP=1, SMAX)       
       ENDDO
!
       Close(188) 

!
! --- Array "B" for vertical displacement ===========================  Vertical displacement U  
        BBBB(:,:)=(0D0,0D0)
 	DO J=1, JMAX 
		DO K=0,NN+1
		BBBB(J,K) = EU(LL(J))*Z(J,K,SMAX)       
		DO P=0, K
		IF(P==0)  BBBB(J,K) = BBBB(J,K) - (Z(J,P,SMAX)-0.           )*BETAU(LL(J),K-P)
		IF(P/=0)  BBBB(J,K) = BBBB(J,K) - (Z(J,P,SMAX)-Z(J,P-1,SMAX))*BETAU(LL(J),K-P)
		ENDDO 
        	BBBB(J,K)= RHOW_O_RHOE_X3*BBBB(J,K)
        	ENDDO
	ENDDO
!
! --- Array "A" for vertical displacement ===========================  Vertical displacement U  
 	AAAA(:,:)=(0D0,0D0)
 	DO J=1, JMAX 
	 	DO K=0,NN+1       
 	 	AAAA(J,K) = EU(LL(J))*IIII(J,K) 
 	 	DO P=0, K
 	 	IF(P==0) AAAA(J,K) = AAAA(J,K) - (IIII(J,P)-0.         )*BETAU(LL(J),K-P) 	 
 	 	IF(P/=0) AAAA(J,K) = AAAA(J,K) - (IIII(J,P)-IIII(J,P-1))*BETAU(LL(J),K-P)
 		ENDDO 
 		AAAA(J,K)= RHOI_O_RHOE_X3*AAAA(J,K)
 		ENDDO
	ENDDO
!
! --- Vertical displacement =========================================  Vertical displacement U 
	U(:,:) = AAAA(:,:) + BBBB(:,:)
!
!
! ............
! ............
! ............
! ............
!
! --- Array "B" for horizontal displacement  ========================  Horizontal displacement V  
       BBBB(:,:)=(0D0,0D0)
 	DO J=1, JMAX 
		DO K=0,NN+1
		BBBB(J,K) = EV(LL(J))*Z(J,K,SMAX)       
		DO P=0, K
		IF(P==0)  BBBB(J,K) = BBBB(J,K) - (Z(J,P,SMAX)-0.           )*BETAV(LL(J),K-P)
		IF(P/=0)  BBBB(J,K) = BBBB(J,K) - (Z(J,P,SMAX)-Z(J,P-1,SMAX))*BETAV(LL(J),K-P)
		ENDDO 
        	BBBB(J,K)= RHOW_O_RHOE_X3*BBBB(J,K)
        	ENDDO
	ENDDO
!
! --- Array "A" for horizontal displacement  ========================  Horizontal displacement V  
       AAAA(:,:)=(0D0,0D0)
 	DO J=1, JMAX 
	 	DO K=0,NN+1       
 	 	AAAA(J,K) = EV(LL(J))*IIII(J,K) 
 	 	DO P=0, K
 	 	IF(P==0) AAAA(J,K) = AAAA(J,K) - (IIII(J,P)-0.         )*BETAV(LL(J),K-P) 	 
 	 	IF(P/=0) AAAA(J,K) = AAAA(J,K) - (IIII(J,P)-IIII(J,P-1))*BETAV(LL(J),K-P)
 		ENDDO 
 		AAAA(J,K)= RHOI_O_RHOE_X3*AAAA(J,K)
 		ENDDO
	ENDDO
!
! --- Horizontal displacement =======================================  Horizontal displacement V
	V(:,:) = AAAA(:,:) + BBBB(:,:)
!
!
! ............
! ............
! ............
! ............
!
!
! --- Array "B" for Gravity potential ============================= Potential variation G 
        BBBB(:,:)=(0D0,0D0)
 	DO J=1, JMAX 
		DO K=0,NN+1
		BBBB(J,K) = EN(LL(J))*Z(J,K,SMAX)       
		DO P=0, K
		IF(P==0)  BBBB(J,K) = BBBB(J,K) - (Z(J,P,SMAX)-0.           )*BETAN(LL(J),K-P)
		IF(P/=0)  BBBB(J,K) = BBBB(J,K) - (Z(J,P,SMAX)-Z(J,P-1,SMAX))*BETAN(LL(J),K-P)
		ENDDO 
        	BBBB(J,K)= RHOW_O_RHOE_X3*BBBB(J,K)
        	ENDDO
	ENDDO
!
!
! --- Array "A" for Gravity potential ============================= Potential variation G 
 	AAAA(:,:)=(0D0,0D0)
 	DO J=1, JMAX 
	 	DO K=0,NN+1       
 	 	AAAA(J,K) = EN(LL(J))*IIII(J,K) 
 	 	DO P=0, K
 	 	IF(P==0) AAAA(J,K) = AAAA(J,K) - (IIII(J,P)-0.         )*BETAN(LL(J),K-P) 	 
 	 	IF(P/=0) AAAA(J,K) = AAAA(J,K) - (IIII(J,P)-IIII(J,P-1))*BETAN(LL(J),K-P)
 		ENDDO 
 		AAAA(J,K)= RHOI_O_RHOE_X3*AAAA(J,K)
 		ENDDO
	ENDDO
!
!
! --- Gravity potential (loading) ===================================  Potential variation G  
	G(:,:) = AAAA(:,:) + BBBB(:,:) 
!
! --- "c constant" (total) ==========================================  The "c" constant
        C_CONST(:) =  SE(1,:) - (AAVV(:) + BBVV(:)) 
!
! --- Sea surface variation (total) =================================  Sea surface variation N
	N(:,:) = G(:,:)
	N(1,:) = N(1,:) + C_CONST(:)
!
!
! --- Normalized potential variation (GRACE) ========================  Potential variation G (Grace) 
!
        G_GRACE(:,:) = G(:,:) 
! ---
! --- <<Water height equivalent>> ==================================== Water height equivalent
! ---     
	open(112,file='k.dat',status='unknown') 
	Do i=1, 2 
		Read(112,'(a12)')header
	Enddo 	
        do L=1, LMAX 
		Read(112,*) IJ, KPRIME(L)
		KPRIME(L)=1D0+KPRIME(L)
	enddo
	close(112)
!	
        W(:,:) = (0D0,0D0)
!
        DO J=4, JMAX 
		W(J,:) = (G_GRACE(J,:)/RHOW_O_RHOE_X3)*(2D0*FLOAT(LJ(J))+1D0)/KPRIME(LJ(J))
        ENDDO
! ---
! --- <<Water height equivalent>> ==================================== Water height equivalent
! ---        
!



! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   FA Gravity DRAFT 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   FA Gravity DRAFT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   FA Gravity DRAFT
!
! --- Array "B" for FA Gravity   ! New as of Feb 2012: it was "nn" 
	bbbb(:,:)=0.
 	do j=1, jmax 
		do k=0,NN+1
		bbbb(j,k) = EFA(ll(j))*Z(j,k,SMAX)       
		do p=0, k
		if(p==0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-0.           )*BETAFA(ll(j),k-p)
		if(p/=0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-Z(j,p-1,SMAX))*BETAFA(ll(j),k-p)
		enddo 
        	bbbb(j,k)= RHOW_o_RHOE_X3*100.*1000.*(GRA_REF/RAD_REF)*bbbb(j,k)      
        	enddo
	enddo
!
! --- Array "A" for FA Gravity   ! New as of Feb 2012: it was "nn" 
 	aaaa(:,:)=0.
 	do j=1, jmax 
	 	do k=0,NN+1       
 	 	aaaa(j,k) = EFA(ll(j))*IIII(j,k) 
 	 	do p=0, k
 	 	if(p==0) aaaa(j,k) = aaaa(j,k) - (IIII(j,p)-0.         )*BETAFA(ll(j),k-p) 	 
 	 	if(p/=0) aaaa(j,k) = aaaa(j,k) - (IIII(j,p)-IIII(j,p-1))*BETAFA(ll(j),k-p)
 		enddo 
 		aaaa(j,k)= RHOI_o_RHOE_X3*100.*1000.*(GRA_REF/RAD_REF)*aaaa(j,k)
 		enddo
	enddo
!
! --- FA Gravity gravity anomaly                      (in units of "MILLIGAL")
	FA(:,:) = aaaa(:,:) + bbbb(:,:)
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   FA Gravity DRAFT 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   FA Gravity DRAFT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   FA Gravity DRAFT
!
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   SS Gravity DRAFT 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   SS Gravity DRAFT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   SS Gravity DRAFT
!
! --- Array "B" for SS Gravity   ! New as of Feb 2012: it was "nn" 
	bbbb(:,:)=0.
 	do j=1, jmax 
		do k=0,NN+1
		bbbb(j,k) = ESS(ll(j))*Z(j,k,SMAX)       
		do p=0, k
		if(p==0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-0.           )*BETASS(ll(j),k-p)
		if(p/=0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-Z(j,p-1,SMAX))*BETASS(ll(j),k-p)
		enddo 
        	bbbb(j,k)= RHOW_o_RHOE_X3*100.*1000.*(GRA_REF/RAD_REF)*bbbb(j,k)      
        	enddo
	enddo
!
! --- Array "A" for SS Gravity   ! New as of Feb 2012: it was "nn" 
 	aaaa(:,:)=0.
 	do j=1, jmax 
	 	do k=0,NN+1       
 	 	aaaa(j,k) = ESS(ll(j))*IIII(j,k) 
 	 	do p=0, k
 	 	if(p==0) aaaa(j,k) = aaaa(j,k) - (IIII(j,p)-0.         )*BETASS(ll(j),k-p) 	 
 	 	if(p/=0) aaaa(j,k) = aaaa(j,k) - (IIII(j,p)-IIII(j,p-1))*BETASS(ll(j),k-p)
 		enddo 
 		aaaa(j,k)= RHOI_o_RHOE_X3*100.*1000.*(GRA_REF/RAD_REF)*aaaa(j,k)
 		enddo
	enddo
!
! --- SS Gravity gravity anomaly          (in units of "MILLIGAL")
	SS(:,:) = aaaa(:,:) + bbbb(:,:)
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   SS Gravity DRAFT 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   SS Gravity DRAFT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   SS Gravity DRAFT
!
!
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Load components DRAFT  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Load components DRAFT  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Load components DRAFT  
!
       LO_I(:,:)  = (0d0,0d0)
       LO_O(:,:)  = (0d0,0d0)
       LO_T(:,:)  = (0d0,0d0)
!
! --- Ice  component of load    ======== RIGID EARTH 
	LO_I(:,:) = RHOI*IIII(:,:)
!
!
! --- Ocean component of load   ======== RIGID EARTH
	LO_O(:,:) = RHOW*ZE(:,:)
!
!
! --- Total load (ice + ocean)  ======== RIGID EARTH
	LO_T(:,:) = LO_I(:,:) + LO_O(:,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Load components DRAFT  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Load components DRAFT  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Load components DRAFT  
!
!
3000 CONTINUE 
!
!
!
!
!
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    - Reporting the output arrays -
!           (SH coefficients)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
! --- Sea level change ------------------------- New as of April 2012 
	Out_Filename=filee
 	open(3,file=Out_Filename,status='unknown',form='unformatted') 
        write(3) SE 
	close(3)
! --- Eustatic Sea level change ---------------- New as of April 2012 
!
!
! --- Sea level change  
	Out_Filename=files 
 	open(3,file=Out_Filename,status='unknown',form='unformatted') 
        write(3) S 
	close(3)
!
! --- Vertical displacement 
	Out_Filename=fileu 
 	open(3,file=Out_Filename,status='unknown',form='unformatted') 
        write(3) U 
	close(3) 
!
! --- Geoid undulations   
	Out_Filename=filen 
 	open(3,file=Out_Filename,status='unknown',form='unformatted') 
        write(3) N 
	close(3)
!
! --- Normalized gravity potential (new as of March 2011) ---------- 
	Out_Filename=fileg 
 	open(3,file=Out_Filename,status='unknown',form='unformatted') 
        write(3) G 
	close(3) 
!
! --- Equivalent Water height (new as of August 4, 2012) -----------  In progress.... 
	Out_Filename=filew 
 	open(3,file=Out_Filename,status='unknown',form='unformatted') 
        write(3) W 
	close(3) 
!
! --- Reduced sea level change 
	Out_Filename=filez 
 	open(3,file=Out_Filename,status='unknown',form='unformatted') 
        write(3) Z 
	close(3)
!
! --- Horizontal displacement 
	Out_Filename=filev 
 	open(3,file=Out_Filename,status='unknown',form='unformatted') 
        write(3) V 
	close(3)
!
! --- Free Air gravity anomalies  
	Out_Filename=filefa
 	open(3,file=Out_Filename,status='unknown',form='unformatted') 
        write(3) FA 
	close(3)
!
! --- Solid Surface gravity anomalies 
	Out_Filename=filess
 	open(3,file=Out_Filename,status='unknown',form='unformatted') 
        write(3) SS 
	close(3)
!
! --- Ice component of load                                              <<<< new 
    Out_Filename=fileli 
    open(3,file=Out_Filename,status='unknown',form='unformatted') 
        write(3) LO_I 
    close(3)
!
! --- Ocean component of load          				         <<<< new 
    Out_Filename=filelo 
    open(3,file=Out_Filename,status='unknown',form='unformatted') 
        write(3) LO_O 
    close(3)
!
! --- Total (ice+ocean) load             			         <<<< new 
    Out_Filename=filelt 
    open(3,file=Out_Filename,status='unknown',form='unformatted') 
        write(3) LO_T 
    close(3)
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Write(*,*) "    +++ SLE solved <3 +++" 
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
	End program SLE
!
!
!
!
!
!
