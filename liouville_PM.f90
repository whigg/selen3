!
!
!
!
!
!
  MODULE DATA_FOR_SUB_LIOUVILLE 
  IMPLICIT NONE 
  INCLUDE "data.inc"
!
!INTEGER, PARAMETER :: M=4*NV   ! Number of isostatic modes 
!INTEGER, PARAMETER :: MP=M-1    ! Number of rotational modes 
!
  REAL*8  KLE                      ! Elastic degree 2 k Love number 
  REAL*8  KLF                       ! Fluid degree 2 k Love number  
!
  REAL*8  KTE                      ! Elastic degree 2 k Love number 
  REAL*8  KTF                       ! Fluid degree 2 k Love number  
  REAL*8 KTS   ! k tidal secular 

!
!REAL*8  KL(M)                      ! Residues of degree 2 k Love number    
!REAL*8  S(M)                        ! Degree 2 isostatic poles 
!
  REAL*8  AEP                        ! (Primed) elastic residue of the PMTF 
  REAL*8  ASP                       ! (Primed) secular residue of the PMTF 
!
!  COMPLEX*16 A(MP)                 ! Rotational roots   
!  COMPLEX*16 AAP(MP)              ! (Primed) rotational residues 
!
  END MODULE DATA_FOR_SUB_LIOUVILLE 
!
!
!
!
!
! ****************************************
     SUBROUTINE LIOUVILLE (PSI, PM, DPM)         
! ****************************************
!
     USE DATA_FOR_SUB_LIOUVILLE 
     IMPLICIT NONE 
!
! -----------------------------------------------------------------
! Solves the Liouville equations for polar motion 
! GS November 2012 
!
!     Input:  PSI (0:NN+1): loading excitation function (input)
!     Output: PM  (0:NN+1): polar motion vector (deg)
!             PMD (0:NN+1): time-derivative of PM (deg/Ma)
! ----------------------------------------------------------------- 
!
     CHARACTER*20, PARAMETER :: FILE_PM ="m.dat"
     CHARACTER*20, PARAMETER :: FILE_PMD="m.dot"    
     REAL*8, PARAMETER :: RAD2DEG=180D0/PI
     REAL*8, PARAMETER :: ERAD=6.371D6
     CHARACTER*20 DATE, TIMC
     COMPLEX*16 PSI(0:NN+1),  PM(0:NN+1), DPM(0:NN+1) 
     COMPLEX*16 MFUNC, DMFUNC 
     INTEGER I, J, K 
     REAL*8 M1(0:NN+1),    M2(0:NN+1),    M3(0:NN+1),   & 
            M1DOT(0:NN+1), M2DOT(0:NN+1), DLOD(0:NN+1), & 
	    POLAR_DISP,  POLAR_RATE,  TIME_BP
!
! ###########################################
!
     call DATE_AND_TIME (date,timc)
!
! ---Love numbers data 
!CALL READ_DATA_FOR_SUB_LIOUVILLE 
!
!
! Love numbers data 
! In the previous version, here there was a call to routine "READ_DATA"
! Now, Love numbers are defined by the PREM model have been kindly 
! provided, for the TIDAL and the LOADING case, by Pascal Gegout.-  
!
! ---- Tidals of degree 2:
!
!  HTE=+0.6035424003D+0
  KTE=+0.2981403969D+0 
!
! ---- Loading of degree 2: 
!  HLE=-0.9915810331D+0 
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
! BEP=(1D0+KTE-HTE)*AEP
!
!  
!
     OPEN (55,file=FILE_PM, status='unknown')   
     CALL PRINT_HEADER(55)
!
     OPEN (57,file=FILE_PMD,status='unknown')   
     CALL PRINT_HEADER(57)
!
     DO 5 I=0, NN+1 
!     
!write(*,*) "Da dentro Liouville ", i, psi(i)
!     
     PM (I)=(0D0,0D0)
     DPM(I)=(0D0,0D0)     
!     
        DO 6 K=0, I 
!     
     	PM(I)  = PM(I)  + PSI(K) *  MFUNC(DELTA*DFLOAT(I-K))
!       DPM(I) = DPM(I) + PSI(K) * DMFUNC(DELTA*DFLOAT(I-K))
! 
 6      CONTINUE
!     
 5   CONTINUE 
!
!
!
     DO 7 I=0, NN+1  
!   
        M1(I) = REAL  (PM(I))
        M2(I) = AIMAG (PM(I))
!
!       M1DOT(I) = REAL (DPM(I))
!       M2DOT(I) = AIMAG(DPM(I)) 
!
! Angle between the z-axis and the rotational vector (deg)  
!
        POLAR_DISP = SQRT(M1(I)**2+M2(I)**2)*RAD2DEG
!
!
! Rate of change of the angle above (deg/Ma)   
!  
        POLAR_RATE = ((M1(I)*M1DOT(I) + M2(I)*M2DOT(I))/SQRT(M1(I)**2+M2(I)**2))*RAD2DEG*1000.
!
        POLAR_RATE = 0D0
!
        TIME_BP=FLOAT(NN+1)-FLOAT(I)*DELTA
!
        WRITE(55,'(2(4X, F8.4), 10(1X,E14.5))') FLOAT(I), & 
					     TIME_BP,  & 
	        			     M1(I)*RAD2DEG, & 
					     M2(I)*RAD2DEG, & 
					     ATAN2(M2(I),   M1(I)   )*RAD2DEG,  &
					     POLAR_DISP 
!        
	WRITE(57,'(2(4X, F8.4), 10(1X,E14.5))') FLOAT(I), & 
					     TIME_BP,  & 	
				             M1DOT(I)*RAD2DEG*1000., & 
					     M2DOT(I)*RAD2DEG*1000., & 
					     ATAN2(M2DOT(I),M1DOT(I))*RAD2DEG,  & 
					     POLAR_RATE  
!        
 7    CONTINUE 
!
 Write(55,*) "# Contact: giorgio.spada@gmail.com"
 Write(57,*) "# Contact: giorgio.spada@gmail.com"
!
!
! Work in progress! Work in progress! Work in progress! Work in progress! 
! Work in progress! Work in progress! Work in progress! Work in progress! 
! Work in progress! Work in progress! Work in progress! Work in progress! 
! Work in progress! Work in progress! Work in progress! Work in progress! 
!
!
     CLOSE(55)
     CLOSE(57)
!   
     END SUBROUTINE LIOUVILLE 
!
!
!   
!
!
!
 FUNCTION MFUNC(X)
 USE DATA_FOR_SUB_LIOUVILLE 
 IMPLICIT NONE 
 INTEGER J, JP 
 COMPLEX*16 MFUNC
 REAL*8 X
!
! ---------------------------------------------------------------
! The function "MFUNC" - see my notes ... (those revised) 
!
! NO contribution is accounted for from A-double primed constants 
! See also the GJI benchmark paper and the work with Valentina B.
! GS August 2010 - 
! ---------------------------------------------------------------
!
  MFUNC = AEP
!
!  MFUNC = MFUNC + ASP*X 
!
!  DO J=1, MP 
!      MFUNC = MFUNC + (AAP(J)/A(J))*(EXP(A(J)*X)-1D0) 
!  ENDDO  
!
 END FUNCTION MFUNC 
!
!
!
!
!
!
 FUNCTION DMFUNC(X)
 USE DATA_FOR_SUB_LIOUVILLE 
 IMPLICIT NONE 
 INTEGER J, JP 
 COMPLEX*16 DMFUNC
 REAL*8 X
!
! ---------------------------------------------------------------
! The function "DMFUNC" - see my notes ... (those revised) 
!
! NO contribution is accounted for from A-double primed constants 
! See also the GJI benchmark paper and the work with Valentina B.
! GS August 2010 - 
! ---------------------------------------------------------------
!
  DMFUNC = (0D0,0D0)
!
!  DO J=1, MP 
!      DMFUNC = DMFUNC + AAP(J)*EXP(A(J)*X) 
!  ENDDO  
!
 END FUNCTION DMFUNC 
!
!
!
!
!
!

!
!
!
!
!
     SUBROUTINE PRINT_HEADER(II)
     IMPLICIT NONE
     INTEGER II
     CHARACTER*20 DATE, TIMC
!
     call DATE_AND_TIME (date,timc)
!
     IF(II==55) THEN 
!       
! ---Header for m.dat
 Write(55,*) "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " 
 Write(55,*) '# File <<m.dat>>, created by program LIOUVILLE.F90 on ', & 
              date(1:4), '.', date(5:6), '.', date(7:8), & 
	      '      time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
 Write(55,*) "# Note: the Chandler Wobble is NOT included"
 Write(55,*) "# " 
 Write(55,*) "#     time      time BP        m1             m2           Arg (m)        Mod (m)"	
 Write(55,*) "#      ka         ka          deg            deg             deg            deg  "	  
 Write(55,*) "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " 
!     
     ENDIF
!
     IF(II==57) THEN 

! ---Header for m.dot
 Write(57,*) "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " 
 Write(57,*) '# File <<m.dot>>, created by program LIOUVILLE.F90 on ', & 
              date(1:4), '.', date(5:6), '.', date(7:8), & 
              '      time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
 Write(57,*) "# Note: the Chandler Wobble is NOT included"
 Write(57,*) "# " 
 Write(57,*) "#     time      time BP       m1_dot         m2_dot      Arg (m_dot)    Mod (m_dot)"   
 Write(57,*) "#      ka         ka          deg/Ma         deg/Ma          deg          deg/Ma   "      
 Write(57,*) "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " 
!      
     ENDIF
!   
     END SUBROUTINE PRINT_HEADER 
!
!
!
!
!
!
