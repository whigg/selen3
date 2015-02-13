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
! -------------------------------------------------------------------- December 2012 
  REAL*8  LTE     ! Elastic degree 2 l Love number (tidal) 
  REAL*8  LTF     ! Fluid degree 2 l Love number (tidal) 
  REAL*8  LT(M)   ! Residues of degree 2 l Love number (tidal)   
! -------------------------------------------------------------------- December 2012 
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
       REAL*8 GTD, GTE, GTF, GTV(M) 
       INTEGER I, J, K, IP, JP, KP, IX, J21, J_INDEX
       COMPLEX*16 SFUNC
       COMPLEX*16 PSI(0:NN+1), GUROT (0:NN+1)
       COMPLEX*16 R(1:M,1:MP), G(1:M), F(1:MP)

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
  CALL READ_DATA_FOR_GU_ROT
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
              DO  IP=1, M
                    GTV(IP) = 0D0 
	      ENDDO
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
              DO  IP=1, M
                    GTV(IP) = KT(IP) 
	      ENDDO
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
              DO  IP=1, M
                    GTV(IP) = HT(IP) 
	      ENDDO
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
              DO  IP=1, M
                    GTV(IP) = LT(IP) 
	      ENDDO
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
  BSP = (GTD+GTF)*ASP 
!
!
! --- Auxiliary array "R" -------
!
  R(:,:)=(0D0,0D0)
  DO J =1, M
  DO JP=1, MP
       R(J,JP) = GTV(J)*AP(JP)/(A_ROOT(JP)-S_ROOT(J))
!write(117,*) j, jp, r(j, jp)
  ENDDO
  ENDDO
!
! --- Vector "G(J)" -------------
!
  DO J =1, M
       G(J) = (0D0,0D0)
       DO JP=1, MP 
       G(J) = G(J) + R(J,JP)
       ENDDO
  ENDDO  
!
! --- Vector "F(JP)" ------------
!
  DO JP =1, MP
       F(JP) = (0D0,0D0)
       DO J=1, M 
       F(JP) = F(JP) + R(J,JP)
       ENDDO
  ENDDO  
!
!
  BP(:) = (0D0,0D0)
  DO JP=1, MP 
       BP(JP) = (GTD+GTE)*AP(JP) + F(JP) 
  ENDDO
!
!
  BPP(:) = (0D0,0D0)
  DO J=1, M 
       BPP(J) = (AEP + ASP/S_ROOT(J))*GTV(J) - G(J) 
  ENDDO
!
!
  write(118,*) BEP 
  write(118,*) BSP 
  do i=1, M
  	write(118,*) i, bpp(i) 
  enddo
  do i=1, MP
  	write(118,*) i, bp(i) 
  enddo




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

!GUROT (I)  = GUROT (I) + PSI(K)*SFUNC(DELTA*DFLOAT(I-K))


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
  SUBROUTINE READ_DATA_FOR_GU_ROT
  USE DATA_FOR_GU_ROT
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
!
!!!!!!!!!!!!!!!!!!!! New as of December 2012 (JGR with Karina) ================
!
	call read_junk(2,  1)
	READ(1,*) lte
        write(13,'(A23,1x,D20.8)') '    - l tidal elastic', lte
!
	call read_junk(2,  1)	
	READ(1,*) ltf 
        write(13,'(A23,1x,D20.8)') '    - l tidal fluid  ', ltf	
!
	call read_junk(2,  1)	
        do i=1,  M   	
		read(1,*) lt(i) 
        	write(13,'(A23,1x,D20.8)') '    - l tidal residue', lt(i)		
	enddo	
!	
!!!!!!!!!!!!!!!!!!!! New as of December 2012 (JGR with Karina) ================
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
       close(1) 
!
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
  END SUBROUTINE READ_DATA_FOR_GU_ROT 
!
!
!
!
!
!
