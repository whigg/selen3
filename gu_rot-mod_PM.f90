!
!
!
!
!
!
  MODULE DATA_FOR_GU_ROT
  IMPLICIT NONE 
  INCLUDE "data.inc"
!
!  INTEGER, PARAMETER :: M=4*NV   ! Number of isostatic modes 
!  INTEGER, PARAMETER :: MP=M-1    ! Number of rotational modes 
! 
  REAL*8  KLE     ! Elastic degree 2 k Love number (loading) 
  REAL*8  KLF     ! Fluid degree 2 k Love number (loading) 
!  REAL*8  KL(M)   ! Residues of degree 2 k Love number (loading)   
  REAL*8  KTE     ! Elastic degree 2 k Love number (tidal) 
  REAL*8  KTF     ! Fluid degree 2 k Love number (tidal) 
  REAL*8  KTS   ! k tidal secular 

!  REAL*8  KT(M)   ! Residues of degree 2 k Love number (tidal)   
  REAL*8  HTE     ! Elastic degree 2 h Love number (tidal) 
  REAL*8  HTF     ! Fluid degree 2 h Love number (tidal) 
!  REAL*8  HT(M)   ! Residues of degree 2 h Love number (tidal)   
!  
! -------------------------------------------------------------------- December 2012 
  REAL*8  LTE     ! Elastic degree 2 l Love number (tidal) 
  REAL*8  LTF     ! Fluid degree 2 l Love number (tidal) 
!  REAL*8  LT(M)   ! Residues of degree 2 l Love number (tidal)   
! -------------------------------------------------------------------- December 2012 
!
!  REAL*8  S_ROOT(M)       ! Degree 2 isostatic roots  
  REAL*8  AEP             ! (Primed) elastic residue of the PMTF 
  REAL*8  ASP             ! (Primed) secular residue of the PMTF 
!  COMPLEX*16 A_ROOT(MP)   ! Rotational roots   
!  COMPLEX*16 AP(MP)       ! (Primed) rotational residues 
!
  REAL*8 BEP           ! (Primed) elastic residue of S_rot
  REAL*8 BSP           ! (Primed) secular residue of S_rot
!  COMPLEX*16 BP (MP)   ! (Primed) residue of S_rot  
!  COMPLEX*16 BPP(M)    ! (Doubly primed) secular residue of S_rot
!
  REAL*8, PARAMETER :: ERAD=6.371D6
  REAL*8, PARAMETER :: OMEGA=7.292115d-5
  REAL*8, PARAMETER :: SIGMA_21=OMEGA**2*ERAD**2/SQRT(30D0)
!
  END MODULE DATA_FOR_GU_ROT
!
!
!
!
!
!
! ********************************************
       SUBROUTINE GU_ROT (PSI, IX, GUROT)
! ********************************************
!
       USE DATA_FOR_GU_ROT
       IMPLICIT NONE 
!
!REAL*8 GTD, GTE, GTF, GTV(M) 
 
        REAL*8 GTD, GTE, GTF 
!
       INTEGER I, J, K, IP, JP, KP, IX, J21, J_INDEX
       COMPLEX*16 SFUNC
       COMPLEX*16 PSI(0:NN+1), GUROT (0:NN+1)
!COMPLEX*16 R(1:M,1:MP), G(1:M), F(1:MP)

!
! ------------------------------------------------------------------
! Given the Loading excitation function PSI, GU_ROT computes the 
! degree 2 - order 1 harmonic component of G^{rot} and U^{rot} as 
! a function of time  ---- GS November 2012 
!
!       Input:  PSI    (0:NN+1): loading excitation function (input)
!       Output: GUROT  (0:NN+1): for IX=1: Degree 2 order 1 component of G^{rot}
!                                for IX=2: Degree 2 order 1 component of U^{rot}
!                            
! ------------------------------------------------------------------ 
!
!  IF(ABS(IX).GE.3) then 
!  	Write(*,*) "Error in GU_ROT - Index is out of range"
!	stop
!	Endif
!
!  CALL READ_DATA_FOR_GU_ROT
!
! Love numbers data 
! In the previous version, here there was a call to routine "READ_DATA"
! Now, Love numbers are defined by the PREM model have been kindly 
! provided, for the TIDAL and the LOADING case, by Pascal Gegout.-  
!
! ---- Tidals of degree 2:
!
  HTE=+0.6035424003D+0
  KTE=+0.2981403969D+0 
!
! ---- Loading of degree 2: 
! HLE=-0.9915810331D+0 
  KLE=+0.2353293958D-1 
!
! The tidal secular "k" Love number is also needed. 
! We assume the observed value, according to Lambeck (1980). 
!
  KTS=0.934D0   
!
! The Elastic "A-prime" is 
!
  AEP=(1D0+KLE)*(KTS/(KTS-KTE))
!
! The Elastic "B-prime" is 
!
  BEP=(1D0+KTE-HTE)*AEP
!
!  
!
!
! Centrifugal Potential/gamma: ONLY Direct effect 
! -----------------------------------------------
  IF(IX==11) THEN
!
              GTD =  1D0
              GTE =  0D0  
              GTF =  0D0  
!
!DO  IP=1, M
!GTV(IP) = 0D0 
!ENDDO
!	    
  ENDIF	       
!
!
! Centrifugal Potential/gamma: Direct AND Indirect effect 
! ------------------------------------------------------- 
  IF(IX==1 ) THEN
!
              GTD =  1D0       
              GTE =  KTE  
              GTF =  KTF  
!
!DO  IP=1, M
!GTV(IP) = KT(IP) 
!ENDDO
!
  ENDIF	
!
! Rotational Vertical Displacement: ONLY Indirect effect 
! ------------------------------------------------------
  IF(IX==2) THEN
!
              GTD =  0D0
              GTE =  HTE  
              GTF =  HTF  
!
!DO  IP=1, M
!GTV(IP) = HT(IP) 
!ENDDO
  ENDIF	 	         
!
!
! Rotational Horizontal Displacement: ONLY Indirect effect 
! --------------------------------------------------------
  IF(IX==3) THEN
!
              GTD =  0D0
              GTE =  LTE  
              GTF =  LTF  
!
!DO  IP=1, M
!GTV(IP) = LT(IP) 
!ENDDO
  ENDIF	 	         
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
! --- Primed elastic term -------
!
  BEP = (GTD+GTE)*AEP 
!
!
! --- Primed secular term -------
!
!  BSP = (GTD+GTF)*ASP 
!
!
! --- Auxiliary array "R" -------
!
!R(:,:)=(0D0,0D0)
!DO J =1, M
!DO JP=1, MP
!R(J,JP) = GTV(J)*AP(JP)/(A_ROOT(JP)-S_ROOT(J))
!write(117,*) j, jp, r(j, jp)
!ENDDO
!ENDDO
!
! --- Vector "G(J)" -------------
!
!DO J =1, M
!G(J) = (0D0,0D0)
!DO JP=1, MP 
!G(J) = G(J) + R(J,JP)
!ENDDO
!ENDDO  
!
! --- Vector "F(JP)" ------------
!
!DO JP =1, MP
!F(JP) = (0D0,0D0)
!DO J=1, M 
!F(JP) = F(JP) + R(J,JP)
!ENDDO
!ENDDO  
!
!
!BP(:) = (0D0,0D0)
!DO JP=1, MP 
!BP(JP) = (GTD+GTE)*AP(JP) + F(JP) 
!ENDDO
!
!
!  BPP(:) = (0D0,0D0)
!  DO J=1, M 
!       BPP(J) = (AEP + ASP/S_ROOT(J))*GTV(J) - G(J) 
!  ENDDO
!
!
!  write(118,*) BEP 
!  write(118,*) BSP 
!do i=1, M
!write(118,*) i, bpp(i) 
!enddo
!do i=1, MP
!write(118,*) i, bp(i) 
!enddo




!
! 
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!    Degree 2 and order 1 term of rotational sea level 
! ///////////////////////////////////////////////////////
!
     GUROT(:)=(0D0,0D0)
!
     DO 15 I=0, NN+1 
!     
     GUROT(I)=(0D0,0D0)
!     
          DO 16 K=0, I    
          GUROT (I)  = GUROT (I) + DCONJG(PSI(K)*SFUNC(DELTA*DFLOAT(I-K)))
 16       CONTINUE
!     
     GUROT(I) = (SIGMA_21/GRA_REF)*GUROT(I) 
!     
 15  CONTINUE 
!
!
     END SUBROUTINE GU_ROT
!
!
! Work in progress! Work in progress! Work in progress! Work in progress! 
! Work in progress! Work in progress! Work in progress! Work in progress! 
! Work in progress! Work in progress! Work in progress! Work in progress! 
! Work in progress! Work in progress! Work in progress! Work in progress! 
!
!
!
!
!
!
 FUNCTION SFUNC(X)
 USE DATA_FOR_GU_ROT
 IMPLICIT NONE 
 INTEGER J, JP 
 COMPLEX*16 SFUNC
 REAL*8 X
!
! ---------------------------------------------------------------
! The function "SFUNC" - see my notes ... (those revised) 
!
! NO contribution is accounted for from A-double primed constants 
! See also the GJI benchmark paper and the work with Valentina B.
! GS August 2010 - 
! ---------------------------------------------------------------
!
  SFUNC = BEP
!
!  SFUNC = SFUNC + BSP*X   
!
!  DO J=1, MP 
!      SFUNC = SFUNC +  (BP(J)/A_ROOT(J))*(EXP(A_ROOT(J)*X)-1D0) 
!  ENDDO 
!
!  DO J=1, M 
!      SFUNC = SFUNC + (BPP(J)/S_ROOT(J))*(EXP(S_ROOT(J)*X)-1D0) 
!  ENDDO  
!
 END FUNCTION SFUNC 
!
!
!
!
!
!
