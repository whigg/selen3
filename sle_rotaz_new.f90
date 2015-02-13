!
 include "liouville.f90"
 include "gu_rot-mod.f90"
!
!
!**********************
 PROGRAM SLE 
!**********************
!
 IMPLICIT NONE 
 INCLUDE "data.inc"			    
!
 CHARACTER*12 HEADER 
 INTEGER I, J, K, L, P, IJ, IS, LJ, MJ, LI, DOM, IND, J21, NANCH, ISTEP, J_INDEX
 INTEGER, ALLOCATABLE :: LL(:), MM(:), DM(:), ANC(:), WET(:)     
!
 REAL*4 RHOE
 REAL*4 RHOI_O_RHOE_X3, & 
        RHOW_O_RHOE_X3, & 
	RHOI_O_RHOW 
 REAL*8, ALLOCATABLE :: ALF(:,:), LONP(:), LATP(:), X(:,:), KPRIME(:)  
 REAL*8, ALLOCATABLE :: AAAA_AVE(:), BBBB_AVE(:), Z_ROT_AVE(:), GMU_ROT_AVE(:)  
 REAL*8 RESH, IMSH
!
 REAL*8, PARAMETER :: CC=8.0394d37, & 
                      AA=8.0131d37, & 
		      CMA=CC-AA,    &
		      ERAD=6.371D6, & 
                      C0=(4D0*PI/3D0)*SQRT(6D0/5D0)
!
 REAL*8 ERR(0:NN+1,1:SMAX)
!
 COMPLEX*16, ALLOCATABLE :: LONG_TABLE(:,:), OC(:), IIII(:,:), ZE(:,:) 
 COMPLEX*16, ALLOCATABLE :: HHHH(:,:),    KKKK(:,:)                   
 COMPLEX*16, ALLOCATABLE :: AAAA(:,:),    AAAA_PRIME(:,:)
 COMPLEX*16, ALLOCATABLE :: BBBB(:,:),    BBBB_PRIME(:,:)
 COMPLEX*16, ALLOCATABLE :: Z_ROT(:,:),   Z_ROT_PRIME(:,:)
 COMPLEX*16, ALLOCATABLE :: PSI_LOAD(:), PM(:), DPM(:)   
! 
 COMPLEX*16, ALLOCATABLE :: GMU_ROT(:,:), GMU_ROT_PRIME(:,:) 
 COMPLEX*16, ALLOCATABLE :: U_ROT21(:),   G_ROT21(:), V_ROT21(:) 
 COMPLEX*16, ALLOCATABLE :: G_GRACE(:,:), G_DIRECT_ROT21(:) 
 COMPLEX*16, ALLOCATABLE :: U_ROT(:,:), & 
                            N_ROT(:,:), & 
			    G_ROT(:,:), &
			    V_ROT(:,:)
 COMPLEX*16, ALLOCATABLE :: U_LOAD(:,:), & 
                            G_LOAD(:,:), &
                            V_LOAD(:,:) 			     
 REAL*8 C_CONST(0:NN+1)			    

!
!--- Convolved Green functions -------------------- 
 REAL*4, ALLOCATABLE :: & 
 BETAS(:,:),  ES(:),    &  ! Sea level change  
 BETAN(:,:),  EN(:),    &  ! Geoid undulations 
 BETAU(:,:),  EU(:),    &  ! Vertical displacemnt
 BETAV(:,:),  EV(:),    &  ! Horizontal displacement  
 BETAFA(:,:), EFA(:),   &  ! Free air (FA) gravity anomaly 
 BETASS(:,:), ESS(:)       ! Solid-surface (SS) gravity anomaly                            
!
! -------------------------------------------------
! Output arrays & filenames 
! 
 CHARACTER*22 Out_Filename  
!
! - "Eustatic" Sea level change ------------------- 
 complex*16, allocatable :: SE(:,:)   
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
! - Equivalent water height -----------------------
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
! ********************************************************
! ********************************************************
! ********************************************************
!
!
! --- Getting info about the number of anchor pixels 
  open(1,file='anchor.tmp',status='old')
  read(1,*) nanch  ; close(1)
!
! --- Allocate memory space
 ALLOCATE( LL(JMAX), MM(JMAX), DM(JMAX), ANC(NP), WET(NP) )
 ALLOCATE( ALF(JMAX,NANCH), LONP(NP), LATP(NP), X(NP,0:NN+1) )
 ALLOCATE( AAAA_AVE(0:NN+1), BBBB_AVE(0:NN+1) ) 
 ALLOCATE( GMU_ROT_AVE(0:NN+1), Z_ROT_AVE(0:NN+1) )
 ALLOCATE( LONG_TABLE(0:LMAX,NP), OC(JMAX) )
 ALLOCATE( IIII(JMAX,0:NN+1), ZE(JMAX,0:NN+1), SE(JMAX,0:NN+1) )
 ALLOCATE( AAAA(JMAX,0:NN+1), AAAA_PRIME(JMAX,0:NN+1) )
 ALLOCATE( BBBB(JMAX,0:NN+1), BBBB_PRIME(JMAX,0:NN+1) )
 ALLOCATE( GMU_ROT_PRIME(JMAX,0:NN+1) )
 ALLOCATE( Z_ROT(JMAX,0:NN+1), Z_ROT_PRIME(JMAX,0:NN+1) )
 ALLOCATE( U_ROT(JMAX,0:NN+1), N_ROT(JMAX,0:NN+1), V_ROT(JMAX,0:NN+1) ) 
 ALLOCATE( G_ROT(JMAX,0:NN+1), GMU_ROT(JMAX,0:NN+1) )
 ALLOCATE( HHHH(JMAX,0:NN+1), KKKK(JMAX,0:NN+1) )
 ALLOCATE( BETAS (0:LMAX,0:NN+1), ES (0:LMAX) )
 ALLOCATE( BETAN (0:LMAX,0:NN+1), EN (0:LMAX) )
 ALLOCATE( BETAU (0:LMAX,0:NN+1), EU (0:LMAX) )
 ALLOCATE( BETAV (0:LMAX,0:NN+1), EV (0:LMAX) )     
 ALLOCATE( BETAFA(0:LMAX,0:NN+1), EFA(0:LMAX) )     
 ALLOCATE( BETASS(0:LMAX,0:NN+1), ESS(0:LMAX) )     
 ALLOCATE( S(JMAX,0:NN+1), N(JMAX,0:NN+1), U(JMAX,0:NN+1), G(JMAX,0:NN+1), V(JMAX,0:NN+1) )
 ALLOCATE( U_LOAD(JMAX,0:NN+1), V_LOAD(JMAX,0:NN+1), G_LOAD(JMAX,0:NN+1) )
 ALLOCATE( G_GRACE(JMAX,0:NN+1), G_DIRECT_ROT21(0:NN+1) )
 ALLOCATE( FA(JMAX,0:NN+1), SS(JMAX,0:NN+1) )
 ALLOCATE( LO_I(JMAX,0:NN+1), LO_O(JMAX,0:NN+1), LO_T(JMAX,0:NN+1) )
 ALLOCATE( Z(JMAX,0:NN+1,0:SMAX), W(JMAX,0:NN+1) )
 ALLOCATE( KPRIME(1:LMAX) )
 ALLOCATE( PSI_LOAD(0:NN+1), PM(0:NN+1), DPM(0:NN+1) )
 ALLOCATE( U_ROT21(0:NN+1), V_ROT21(0:NN+1), G_ROT21(0:NN+1) )
!
! ********************************************************
! ********************************************************
! ********************************************************
!
! --- Average Earth density from the TABOO or ALMA log files =====
 IF(CDE.GE. 0) IND=1 
 IF(CDE.EQ.-1) IND=2  
 IF(CDE.GE. 0.or.CDE==-1) CALL AVERAGE_EARTH_DENSITY(IND, RHOE)
! 
 IF(CDE==-2) RHOE=5511.57
!
! === Useful quantities ==========================================
 RHOI_O_RHOE_X3 = 3.*RHOI/RHOE 
 RHOW_O_RHOE_X3 = 3.*RHOW/RHOE 
 RHOI_O_RHOW    =    RHOI/RHOW 
!
! === Pre-computing 'l' and 'm' corresponding to degree 'J' ======
 do j=1, jmax 
     mm(j)=mj(j) 
     ll(j)=lj(j)
     dm(j)=2-dom(j)
 enddo	
!
! === Reading the ALFs table from <<sh.bin>> =====================
 Write(*,*) '    - Reading the ALFs from file sh.bin'
 open(3,file='sh.bin',status='unknown',form='unformatted') 
       read(3)ALF
       read(3)LONG_TABLE
 Close(3) 
!
! === Examining the pixels table & extracting information ========
 Write(*,*) '    - Reading the pixels data from the table'
 Open(1,file='px-table.dat',status='unknown') 
	Do i=1, 4 
		Read(1,'(a12)')header
	Enddo
	Do i=1, np 
		Read (1,*) lonp(i), latp(i), anc(i), k, wet(i) 	
	Enddo
 close(1)
!
! === Reading the SH OF coefficients from shof.dat ===============
 Write(*,*) '    - Reading the SH OF coeff. from shof.dat'
 	open(3,file='shof.dat',status='unknown')
 	do j=1, jmax   
		read(3,*) k, resh, imsh 
                oc(j)=cmplx(resh, imsh)	
 	enddo
 close(3)
!
! === Reading the <<I>> array ====================================
 Write(*,*) "    - Reading array 'I'"
 open(1,file='shice.dat',status='unknown') 
 read(1,*) IIII
 close(1)
!
! === Reading the E and BETA arrays ==============================
 Write(*,*) "    - Reading arrays 'E' and 'Beta'"
 open(1,file= 'ebs.dat',status='unknown')    
 open(2,file= 'ebu.dat',status='unknown')    
 open(3,file= 'ebn.dat',status='unknown')    
 open(4,file= 'ebv.dat',status='unknown')    
 open(7,file='ebfa.dat',status='unknown')    
 open(8,file='ebss.dat',status='unknown')    
!
 do li=0, lmax 
 read(1,*) l, Es(l),  ( betas(l,k), k=0,nn+1) 
 read(2,*) l, Eu(l),  ( betau(l,k), k=0,nn+1)  
 read(3,*) l, En(l),  ( betan(l,k), k=0,nn+1)  
 read(4,*) l, Ev(l),  ( betav(l,k), k=0,nn+1)  
 read(7,*) l, Efa(l), (betafa(l,k), k=0,nn+1)  
 read(8,*) l, Ess(l), (betass(l,k), k=0,nn+1)  
 enddo
  close(8)
   close(7)
    close(4)
     close(3)
      close(2)
       close(1) 
!
!
! --- Computing the eustatic Z array... ===========================
	ZE(:,:) = (0D0,0D0)  		
	DO K=0,NN+1	
		ZE(:,K) = - RHOI_O_RHOW*(IIII(1,K)/OC(1))*OC(:)
	ENDDO
!
! --- Computing the eustatic S array... ===========================
	SE(:,:) = (0D0,0D0)
	SE(1,:) = - RHOI_O_RHOW*(IIII(1,:)/OC(1)) 
!
! === Computing the A array... ====================================
!
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
	AAAA(:,:)= RHOI_O_RHOE_X3*AAAA(:,:)

!
! === Computing the ocean average of A ============================ 
	DO K=0, NN+1 
		AAAA_AVE(K)=0D0
		DO J=1, JMAX 
		AAAA_AVE(K) = AAAA_AVE(K) + & 
		          DM(J)*REAL(OC(J)*CONJG(AAAA(J,K)))/OC(1)  
		ENDDO
	ENDDO 
!
! === Modified array A ============================================
	AAAA_PRIME(:,:) = AAAA(:,:)
	AAAA_PRIME(1,:) = AAAA(1,:)-AAAA_AVE(:) 
!
! === Computing the H-array... ====================================
!
        X(:,:) = 0D0
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J,K) &
!$OMP     SHARED(X,ALF,ANC,DM,AAAA_PRIME,LONG_TABLE,MM) &
!$OMP     SCHEDULE(GUIDED)
	DO I=1, NP    
	DO K=0, NN+1
	        DO J=1, JMAX  
        	X(I,K) = X(I,K) + & 
			 ALF(J,ANC(I))*DM(J)*REAL(AAAA_PRIME(J,K)*LONG_TABLE(MM(J),I)) 
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
		    DO K=0, NN+1
			HHHH(J,K) = HHHH(J,K) + & 
			X(I,K)*ALF(J,ANC(I))*CONJG(LONG_TABLE(MM(J),I))  
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
! === Initializing the Z and S arrays =============================== 
	Z(:,:,0) = ZE(:,:) 
	S(:,:)   = SE(:,:) 
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
		BBBB_AVE(K)=0D0
		DO J=1, JMAX 
		BBBB_AVE(K) = BBBB_AVE(K) + & 
		          DM(J)*REAL(OC(J)*CONJG(BBBB(J,K)))/OC(1)
		ENDDO
	ENDDO 
!
! === COMPUTING MODIFIED 'B' ARRAY ================================== RECURSION
	BBBB_PRIME(:,:)=BBBB(:,:)
	BBBB_PRIME(1,:)=BBBB(1,:)-BBBB_AVE(:) 
!
! === Computing array K...  ========================================= Recursion
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,K,J) &
!$OMP    SHARED(X,ALF,ANC,DM,BBBB_PRIME,LONG_TABLE,MM) SCHEDULE(GUIDED)
	DO I=1, NP 
 	   DO K=0, NN+1
              X(I,K)=0D0
	      DO J=1, JMAX  
        	  X(I,K) = X(I,K) + & 
			 ALF(J,ANC(I))*DM(J)*REAL(BBBB_PRIME(J,K)*LONG_TABLE(MM(J),I)) 
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
	          DO K=0,NN+1           ! it was "1,nn+1" I do not know why... 
		      KKKK(J,K) = KKKK(J,K) + X(I,K)*ALF(J,ANC(I))*CONJG(LONG_TABLE(MM(J),I))  
		  ENDDO
	       END IF
	    ENDDO
	ENDDO
!$OMP END PARALLEL DO	 
	KKKK(:,:)=KKKK(:,:)/FLOAT(NP) 
!
! === Computing the polar motion Excitation Function (E.F.) ========= Recursion
!
       PSI_LOAD(:)=(0D0,0D0)
!
       DO K=0, NN+1 
!
       J21=J_INDEX(2,1)
!
            IF(K.EQ.0) PSI_LOAD(K) =  RHOI*IIII(J21,K) + RHOW*Z(J21,K,IS-1)
!	   
            IF(K.GE.1) PSI_LOAD(K) =  RHOI*(IIII(J21,K) -     IIII(J21,K-1))  +  &  	  
                                      RHOW*(   Z(J21,K,IS-1) -   Z(J21,K-1,IS-1)) 
!
            PSI_LOAD(K)= DCONJG(PSI_LOAD(K))*(C0*ERAD**4)/CMA      
!
       ENDDO
!
!
! === Computing the G^rot and U^rot terms of the SLE ================ Recursion
!
       G_ROT(:,:)   = (0D0,0D0)
       U_ROT(:,:)   = (0D0,0D0)
       V_ROT(:,:)   = (0D0,0D0)
!       
       GMU_ROT(:,:) = (0D0,0D0)
!
       G_ROT21(:)   =  0D0
       U_ROT21(:)   =  0D0 
       V_ROT21(:)   =  0D0 
!
       CALL GU_ROT (PSI_LOAD, 1, G_ROT21)      
       CALL GU_ROT (PSI_LOAD, 2, U_ROT21)      
       CALL GU_ROT (PSI_LOAD, 3, V_ROT21)      
!      
       J21=J_INDEX(2,1)
! 
       G_ROT(J21,:)   = G_ROT21(:)
       U_ROT(J21,:)   = U_ROT21(:)
       V_ROT(J21,:)   = V_ROT21(:)
!
       GMU_ROT(J21,:) = G_ROT(J21,:)-U_ROT(J21,:) 
!
! --- Computing the ocean average of (G-U)^{rot}  ================== Recursion  
!
	DO K=0, NN+1 
		GMU_ROT_AVE(K)=0D0
		DO J=1, JMAX 
		GMU_ROT_AVE(K) = GMU_ROT_AVE(K) + & 
		            DM(J)*REAL(OC(J)*CONJG(GMU_ROT(J,K)))/OC(1)  
		ENDDO
	ENDDO 
	write(*,*) "    - [new] ocean average of (G-U)^rot"
!
! --- Computing array (G-U)^{rot,prime} ============================ Recursion
!
	write(*,*) "    - [new] Modified (G-U)^rot array"
!
	GMU_ROT_PRIME(:,:) = GMU_ROT(:,:)
	GMU_ROT_PRIME(1,:) = GMU_ROT(1,:)-GMU_ROT_AVE(:) 
!
!
! --- Computing GMU^{rot, prime} at all pixels ===================== Recursion
!
	X(:,:)=0D0
	write(*,*) "    - [new] (G-U)^{rot,prime} at all pixels"
!	
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J,K) &
!$OMP       SHARED(X,ALF,ANC,LONG_TABLE,GMU_ROT_PRIME,MM,DM) &
!$OMP       SCHEDULE(GUIDED)
	DO I=1, NP    
		DO K=0, NN+1   	
		X(I,K)=0D0 	
		DO J=1, JMAX 	
			X(I,K) = X(I,K) + & 
                                 ALF(J,ANC(I))*DM(J)*REAL(GMU_ROT_PRIME(J,K)*LONG_TABLE(MM(J),I)) 	
		ENDDO
		ENDDO		
	ENDDO
!$OMP END PARALLEL DO
!
!
! --- Harmonics of Z^{rot} by integration over wet pixels =========== Recursion ---------------------
!
        Z_ROT(:,:)=(0D0,0D0)     
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J,K) &
!$OMP       SHARED(WET,Z_ROT,X,ALF,ANC,LONG_TABLE,MM) &
!$OMP       SCHEDULE(GUIDED)
	DO I=1,NP
	   IF(WET(I)==1) THEN
	       DO K=0, NN+1   
	           DO J=1, JMAX 
	                 Z_ROT(J,K) = Z_ROT(J,K) + X(I,K)*ALF(J,ANC(I))*CONJG(LONG_TABLE(MM(J),I))  
	           ENDDO
	       END DO
	   END IF
	ENDDO
!$OMP END PARALLEL DO	
!
	Z_ROT(:,:)=Z_ROT(:,:)/FLOAT(NP) 
!
!
! --- Solving for array 'Z' st step "is" =========================== Recursion
!
	write(*,*) "    - [new] New iterate of array Z "
!
	Z(:,:,IS) = HHHH(:,:) + KKKK(:,:) + Z_ROT(:,:)
!	
!
! --- Solving for array 'S' ======================================== Recursion
!
	S(:,:)    = AAAA_PRIME(:,:) + SE(:,:) + BBBB_PRIME(:,:) + GMU_ROT_PRIME(:,:)
!
!
! ------------------------------
! ---    End of recursion    ---
! ------------------------------
!
1000 CONTINUE  
!  
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
! ---------------------------------------------------------------------------------
! ### In PROGRESS ### In PROGRESS ### In PROGRESS ### In PROGRESS ### In PROGRESS
! ---------------------------------------------------------------------------------
! 
!
!
!
! === Computing the polar motion vector by the Liouville eqs. =======

       CALL LIOUVILLE (PSI_LOAD, PM, DPM)
!
! ............
! ............
! ............
! ............
!
!
! --- Eustatic solution: U=0, S=N, V=0 ==============================
        if(imode==3.or.smax==0) then 
				U(:,:) = 0D0
				S(:,:) = SE(:,:)
				N(:,:) = S(:,:)
				V(:,:) = 0D0 
				goto 3000
				endif
!
! ............
! ............
! ............
! ............
!
! --- Array "B" for vertical displacement ===========================  Vertical displacement U (loading) 
       BBBB(:,:)=(0D0,0D0)
       DO J=1, JMAX 
	       DO K=0,NN+1
	       BBBB(J,K) = EU(LL(J))*Z(J,K,SMAX)       
	       DO P=0, K
	       IF(P==0)  BBBB(J,K) = BBBB(J,K) - (Z(J,P,SMAX)-(0D0,0D0)   )*BETAU(LL(J),K-P)
	       IF(P/=0)  BBBB(J,K) = BBBB(J,K) - (Z(J,P,SMAX)-Z(J,P-1,SMAX))*BETAU(LL(J),K-P)
	       ENDDO 
	       ENDDO
       ENDDO
       BBBB(:,:)= RHOW_O_RHOE_X3*BBBB(:,:)
!
!
! --- Array "A" for vertical displacement ===========================  Vertical displacement U (loading) 
 	AAAA(:,:)=(0D0,0D0)
 	DO J=1, JMAX 
	 	DO K=0,NN+1       
 	 	AAAA(J,K) = EU(LL(J))*IIII(J,K) 
 	 	DO P=0, K
 	 	IF(P==0) AAAA(J,K) = AAAA(J,K) - (IIII(J,P)-(0D0,0D0)  )*BETAU(LL(J),K-P) 	 
 	 	IF(P/=0) AAAA(J,K) = AAAA(J,K) - (IIII(J,P)-IIII(J,P-1))*BETAU(LL(J),K-P)
 		ENDDO 
 		ENDDO
	ENDDO
        AAAA(:,:)= RHOI_O_RHOE_X3*AAAA(:,:)
!
!
! --- Vertical displacement from Loading ===========================  Vertical displacement U (loading) 
	U_LOAD(:,:) = AAAA(:,:) + BBBB(:,:) 

! --- Vertical displacement (total) ================================  Vertical displacement U (loading + rot) 
	U(:,:) = U_LOAD(:,:) + U_ROT(:,:) 

!
! ............
! ............
! ............
! ............
!
! --- Array "B" for horizontal displacement ===========================  Horizontal displacement V (loading) 
       BBBB(:,:)=(0D0,0D0)
       DO J=1, JMAX 
	      DO K=0,NN+1
	      BBBB(J,K) = EV(LL(J))*Z(J,K,SMAX)       
	      DO P=0, K
	      IF(P==0)  BBBB(J,K) = BBBB(J,K) - (Z(J,P,SMAX)-(0D0,0D0)   )*BETAV(LL(J),K-P)
	      IF(P/=0)  BBBB(J,K) = BBBB(J,K) - (Z(J,P,SMAX)-Z(J,P-1,SMAX))*BETAV(LL(J),K-P)
	      ENDDO 
	      ENDDO
       ENDDO
       BBBB(:,:)= RHOW_O_RHOE_X3*BBBB(:,:)
!
! --- Array "B" for horizontal displacement ===========================  Horizontal displacement V (loading)  
       AAAA(:,:)=(0D0,0D0)
       DO J=1, JMAX 
	       DO K=0,NN+1	 
	       AAAA(J,K) = EV(LL(J))*IIII(J,K) 
	       DO P=0, K
	       IF(P==0) AAAA(J,K) = AAAA(J,K) - (IIII(J,P)-(0D0,0D0)  )*BETAV(LL(J),K-P)	
	       IF(P/=0) AAAA(J,K) = AAAA(J,K) - (IIII(J,P)-IIII(J,P-1))*BETAV(LL(J),K-P)
	       ENDDO 
	       ENDDO
       ENDDO
	AAAA(:,:)= RHOI_O_RHOE_X3*AAAA(:,:)
!
! --- Vertical displacement from Loading ===========================  Horizontal displacement V (loading) 
       V_LOAD(:,:) = AAAA(:,:) + BBBB(:,:) 
!
! --- Vertical displacement (total) ================================  Horizontal displacement V (loading + rot) 
       V(:,:) = V_LOAD(:,:) + V_ROT(:,:) 

!
! ............
! ............


! ............
! ............
!
! 
! --- Array "B" for Gravity potential ===============================  Potential variation G  (loading)
       BBBB(:,:)=(0D0,0D0)
       DO J=1, JMAX 
	       DO K=0,NN+1
	       BBBB(J,K) = EN(LL(J))*Z(J,K,SMAX)       
	       DO P=0, K
	       IF(P==0)  BBBB(J,K) = BBBB(J,K) - (Z(J,P,SMAX)-(0D0,0D0)	   )*BETAN(LL(J),K-P)
	       IF(P/=0)  BBBB(J,K) = BBBB(J,K) - (Z(J,P,SMAX)-Z(J,P-1,SMAX))*BETAN(LL(J),K-P)
	       ENDDO 
	       ENDDO
       ENDDO
       BBBB(:,:)= RHOW_O_RHOE_X3*BBBB(:,:)
!
!
! --- Array "A" for Gravity potential ===============================  Potential variation G  (loading)
 	AAAA(:,:)=(0D0,0D0)
 	DO J=1, JMAX 
	 	DO K=0,NN+1       
 	 	AAAA(J,K) = EN(LL(J))*IIII(J,K) 
 	 	DO P=0, K
 	 	IF(P==0) AAAA(J,K) = AAAA(J,K) - (IIII(J,P)-(0D0,0D0)  )*BETAN(LL(J),K-P) 	 
 	 	IF(P/=0) AAAA(J,K) = AAAA(J,K) - (IIII(J,P)-IIII(J,P-1))*BETAN(LL(J),K-P)
 		ENDDO 
 		ENDDO
	ENDDO
        AAAA(:,:)= RHOI_O_RHOE_X3*AAAA(:,:)
!
!
! --- Gravity potential (loading) ===================================  Potential variation G  (loading)
	G_LOAD(:,:) = AAAA(:,:) + BBBB(:,:) 
!
! --- Gravity potential (total) =====================================  Vertical displacement U (loading + rot) 
	G(:,:) = G_LOAD(:,:)  + G_ROT(:,:) 
!
! --- "c constant" (total) ==========================================  The "c" constant (loading + rot)
        C_CONST(:) =  SE(1,:) - (AAAA_AVE(:) + BBBB_AVE(:)) - GMU_ROT_AVE(:)
!
! --- Sea surface variation (total) =================================  Sea surface variation (loading + rot) 
	N(:,:) = G(:,:)
	N(1,:) = N(1,:) + C_CONST(:)
!
!
! --- Normalized potential variation (GRACE) ========================  Potential variation G (Grace) 
!
        G_GRACE(:,:) = G(:,:) 
!
        CALL GU_ROT (PSI_LOAD, 11, G_DIRECT_ROT21)  
!
        J21=J_INDEX(2,1)
!
        G_GRACE(J21,:) = G_GRACE(J21,:) - G_DIRECT_ROT21(:) 
!
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

   

! ---
!
!
3000 CONTINUE 
!
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
!
! --- Eustatic Sea level change ---------------- New as of April 2012 
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
        write(3) G_GRACE 
	close(3) 
!
! --- Equivalent Water height (new as of August 4, 2012) -----------  
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
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Write(*,*) "    +++ SLE solved <3 +++" 
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!

	End program SLE
!
!
!
