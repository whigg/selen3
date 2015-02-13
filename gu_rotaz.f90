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
  INTEGER, PARAMETER :: M=4*NV   ! Number of isostatic modes 
  INTEGER, PARAMETER :: MP=M-1    ! Number of rotational modes 
! 
  REAL*8  KLE     ! Elastic degree 2 k Love number (loading) 
  REAL*8  KLF     ! Fluid degree 2 k Love number (loading) 
  REAL*8  KL(M)   ! Residues of degree 2 k Love number (loading)   
  REAL*8  KTE     ! Elastic degree 2 k Love number (tidal) 
  REAL*8  KTF     ! Fluid degree 2 k Love number (tidal) 
  REAL*8  KT(M)   ! Residues of degree 2 k Love number (tidal)   
  REAL*8  HTE     ! Elastic degree 2 h Love number (tidal) 
  REAL*8  HTF     ! Fluid degree 2 h Love number (tidal) 
  REAL*8  HT(M)   ! Residues of degree 2 h Love number (tidal)   
!
  REAL*8  S_ROOT(M)       ! Degree 2 isostatic roots  
  REAL*8  AEP             ! (Primed) elastic residue of the PMTF 
  REAL*8  ASP             ! (Primed) secular residue of the PMTF 
  COMPLEX*16 A_ROOT(MP)   ! Rotational roots   
  COMPLEX*16 AP(MP)       ! (Primed) rotational residues 
!
  REAL*8 BEP           ! (Primed) elastic residue of S_rot
  REAL*8 BSP           ! (Primed) secular residue of S_rot
  COMPLEX*16 BP (MP)   ! (Primed) residue of S_rot  
  COMPLEX*16 BPP(M)    ! (Doubly primed) secular residue of S_rot
!
  REAL*8, PARAMETER :: ERAD=6.371D6
  REAL*8, PARAMETER :: OMEGA=7.292115d-5
  REAL*8, PARAMETER :: SIGMA_21=OMEGA**2*ERAD**2/SQRT(30D0)
!
  END MODULE DATA_FOR_GU_ROTAZ
!
!
!
!
!
! ********************************************
       SUBROUTINE GU_ROTAZ (PSI, IX, GUROT)
! ********************************************
!
       USE DATA_FOR_GU_ROTAZ
       IMPLICIT NONE 
!
       REAL*8 GTE, GTV(M) 
       INTEGER I, J, K, IP, JP, KP, IX, J21, J_INDEX
       COMPLEX*16 SFUNC, P(MP,M), E(M), F(MP)  
       COMPLEX*16 PSI(0:NN+1), GUROT (0:NN+1)
!
!
! ------------------------------------------------------------------
! Given the Loading excitation function PSI, S_ROTAZ computes the 
! degree 2 - order 1 harmonic component of G^{rot} and U^{rot} as 
! a function of time  ---- GS November 2012 
!
!       Input:  PSI    (0:NN+1): loading excitation function (input)
!       Output: GURO T (0:NN+1): for IX=1: Degree 2 order 1 component of G^{rot}
!                                for IX=2: Degree 2 order 1 component of U^{rot}
!                            
! ------------------------------------------------------------------ 
!
  IF(IX.GE.4) then 
  	Write(*,*) "Error in SUN_ROTAZ - Index is out of range"
	stop
	Endif
!
! --- Elastic  component  of the Tidal Sea Level Green function
  IF(IX==1) GTE =  1D0 + KTE    
  IF(IX==2) GTE =      + HTE
!
!
! --- Viscous components of the Tidal Sea Level Green function
  DO  IP=1, M
     IF(IX==1) GTV(IP) = KT(IP) 
     IF(IX==2) GTV(IP) = HT(IP) 
  ENDDO
!
! --- Primed elastic term 
  BEP = GTE*AEP 
!
! --- Primed secular term 
  BSP = GTE*ASP 
  DO IP=1, M 
     BSP = BSP - ASP*GTV(IP)/S_ROOT(IP) 
  ENDDO
!write(*,*) "BSP ", bsp
!
! --- Auxiliary array 
  P(:,:)=(0D0,0D0)
  DO I =1, MP 
  DO IP=1, M 
     P(I,IP) = - GTV(IP)*AP(I)/(A_ROOT(I)-S_ROOT(IP)) 
  ENDDO
  ENDDO
!
! --- Primed terms 
  DO I=1, MP 
     F(I)=(0D0,0D0)
     DO IP= 1, M 	
        F(I)=F(I) - P(I,IP)
     ENDDO
  ENDDO
!
  DO I=1, MP 
     BP(I) = GTE*AP(I) + F(I) 
  ENDDO
!
! --- Doubly-primed terms 
  DO IP=1, M 
     E(IP)=(0D0,0D0)
     DO I= 1, MP 	
        E(IP)=E(IP) + P(I,IP)
     ENDDO
  ENDDO
!
  DO IP=1, M
     BPP(IP) = AEP*GTV(IP) + E(IP) + ASP*GTV(IP)/S_ROOT(IP) 		
  ENDDO
!
!
! 
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!    Degree 2 and order 1 term of rotational sea level 
! ///////////////////////////////////////////////////////
!
     GUROT (:) = (0D0,0D0)
!
     DO 15 I=0, NN +1 
!     
     GUROT (i)= (0D0,0D0)
!     
          DO 16 K=0, I    
          GUROT (I)  = GUROT (I) + DCONJG(PSI(K)*SFUNC(DELTA*DFLOAT(I-K)))
 16       CONTINUE
!     
          GUROT (I) = (SIGMA_21/GRA_REF)*GUROT(I) 
!     
 15   CONTINUE 
!
       END SUBROUTINE GU_ROTAZ
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
 USE DATA_FOR_SUB_SUN_ROTAZ
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
  SFUNC = SFUNC + BSP*X   
!
  DO J=1, MP 
      SFUNC = SFUNC +  (BP(J)/A_ROOT(J))*(EXP(A_ROOT(J)*X)-1D0) 
  ENDDO 
!
  DO J=1, M 
      SFUNC = SFUNC + (BPP(J)/S_ROOT(J))*(EXP(S_ROOT(J)*X)-1D0) 
  ENDDO  
!
 END FUNCTION SFUNC 
!
!
!
!
!
!
  SUBROUTINE READ_DATA_FOR_GU_ROTAZ
  USE DATA_FOR_GU_ROTAZ
  IMPLICIT NONE 
!
! This subroutine has the purpose of reading data from external files
! and computing some useful combinations of Love numbers. GS AUG 2010 
!
  CHARACTER*22, PARAMETER :: FILEDEG2='Deg_2_Love_numbers.dat'  
  CHARACTER*22, PARAMETER :: FILEPMTF='PM_time_domain.dat'  
  CHARACTER*50 JUNKA
  REAL*8 REA, IMA, REAAP, IMAAP
  INTEGER :: I, J, JP, K 
!
!
  Write(*,*) '    - Reading data from external files'
! 
  Write(*,*) '    - A backup of data is reported on file <<Rotational_data.dat>>'
!
  Open(13,file='Rotational_data.dat',status='unknown') 
!  
  if(TLOVE==1) then 
!
  Write(*,*) '    - Reading *loading* degree 2 stuff from file ', trim(adjustl(FILEDEG2)) 
!                   [this file is created by the high-precision TABOO code]    
!
 	open (1,file=FILEDEG2,status='unknown') 
!
	call read_junk(12, 1)
	READ(1,*) kle
        write(13,'(A23,1x,D20.8)') '    - k loading elastic', kle
!
	call read_junk(2,  1)	
	READ(1,*) klf 
        write(13,'(A23,1x,D20.8)') '    - k loading fluid  ', klf	
!
	call read_junk(2,  1)	
        do i=1,  M   	
	read(1,*) kl(i) 
        write(13,'(A23,1x,D20.8)') '    - k loading residue', kl(i)		
	enddo
!
	call read_junk(2,    1)
	call read_junk(M+1,  1)
	call read_junk(2,    1)
        do i=1,  M   	
		read (1,*) s_root(i) 
        	write(13,'(A23,1x,D20.8)') '    - s isostatic roots', s_root(i)		                         		
	enddo
!
  Write(*,*) '    - Reading *tidal* degree 2 stuff from file ', trim(adjustl(FILEDEG2)) 
!                   [this file is created by the high-precision TABOO code]    
!
	call read_junk(2,  1)
	READ(1,*) hte
        write(13,'(A23,1x,D20.8)') '    - h tidal elastic', hte
!
	call read_junk(2,  1)	
	READ(1,*) htf 
        write(13,'(A23,1x,D20.8)') '    - h tidal fluid  ', htf	
!
	call read_junk(2,  1)	
        do i=1,  M   	
		read(1,*) ht(i) 
        	write(13,'(A23,1x,D20.8)') '    - h tidal residue', ht(i)		
	enddo
!
	call read_junk(2,  1)
	READ(1,*) kte
        write(13,'(A23,1x,D20.8)') '    - k tidal elastic', kte
!
	call read_junk(2,  1)	
	READ(1,*) ktf 
        write(13,'(A23,1x,D20.8)') '    - k tidal fluid  ', ktf	
!
	call read_junk(2,  1)	
        do i=1,  M   	
		read(1,*) kt(i) 
        	write(13,'(A23,1x,D20.8)') '    - k tidal residue', kt(i)		
	enddo
!	
        close(1) 
 else 
	Write(*, *) "The polar motion in response to loading can  only be computed"
	Write(*, *) "if file <<Deg_2_Love_numbers.dat>> exists. The program will STOP "
	Write(88,*) "The polar motion in response to loading can  only be computed"
	Write(88,*) "if file <<Deg_2_Love_numbers.dat>> exists. The program will STOP "
        CALL STOP_CONFIG
 endif
!
!
!
!
!
 Write(*,*) '    - Reading PMTF data from file ', trim(adjustl(FILEPMTF)) 
!
 if(TLOVE==1) then 
!
 	open (1,file=FILEPMTF,status='unknown') 

	call read_junk(9, 1)
!
	call read_junk(2, 1)
	READ(1,'(A4,1X,2(1X,D20.8))') JUNKA, AEP 
        write(13,'(A23,4(1x,D20.8))') '    - *Primed* Elastic R residue', AEP
!
	call read_junk(2, 1)
	READ(1,'(A4,1X,2(1X,D20.8))') JUNKA, ASP
        write(13,'(A23,4(1x,D20.8))') '    - *Primed* Secular S residue', ASP
!
	call read_junk(2, 1)
  	DO K=1, MP
		READ(1,'(I4,1X,8(1X,E20.8))')  I, REA, IMA, REAAP, IMAAP 
!		
		a_root(k) = dcmplx(rea,  ima)
		ap(k)     = dcmplx(reaap, imaap)
!
        write(13,'(A23,4(1x,D20.8))') '    - R roots & *primed* residues', a_root(k), ap(k)
!
	ENDDO	
!	
        close(1) 
 else 
	Write(*, *) "The polar motion in response to loading can  only be computed"
	Write(*, *) "if the PMTF file exists. ****** The program will STOP ****** "
	Write(88,*) "The polar motion in response to loading can  only be computed"
	Write(88,*) "if the PMTF file exists. ****** The program will STOP ****** "
        CALL STOP_CONFIG
 endif
!
  CLOSE(13)  
!
  END SUBROUTINE READ_DATA_FOR_GU_ROTAZ 
!
!
!
!
!
!














