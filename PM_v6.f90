!
! 
! ***************************************************************************
! Program Polar Motion - version 6 
!
! Solves the Liouville equations for GIA- COST BENCHMARK '10
! ***************************************************************************
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
! Copyright (C) 2008 Giorgio Spada, Florence Colleoni, and Paolo Stocchi 
!
! This file is part of SELEN. 
!  
! SELEN is free software: you can redistribute it and/or modify it under the 
! terms of the GNU General Public License as published by the Free Software 
! Foundation, either version 3 of the License, or at your option any later 
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
! Last update by GS on July 25, 2010 
! Some CLEANING by GS on November 22, 2012.
!
!
! WARNING: The Chandler wobble is filtered out SINCE THE ONSET of computations
!          Please act on switch "OPT" below to change this default
!
!
!
! Input files:
!    - "Deg_2_Love_numbers.dat":  Degree 2 Love numbers from TBHP.F90 
!
!
! Output files:
!    - "PM_spectrum.dat":    Table of PMTF coefficients 
!    - "PM_time_domain.dat": Table of PRIMED PMTF coefficients 
!
!
!
! **********************************************************
! **********************************************************
! **********************************************************
!
!
!
  MODULE DATA 
  IMPLICIT NONE 
  INCLUDE "data.inc"
!
  Integer, parameter  :: QP = kind(1.d0)
!
  INTEGER, PARAMETER :: LDEG=2 
  INTEGER, PARAMETER :: M=4*NV  
  INTEGER, PARAMETER :: MAX_MODES=M	
  REAL(qp), PARAMETER :: CCCC=8.0394D37, AAAA=8.0131D37, OMEGA=7.292115D-5
  REAL(qp), PARAMETER :: SECY=365.25D0*24D0*3600D0*1000D0
  REAL(qp), PARAMETER :: SR=((CCCC-AAAA)/AAAA)*OMEGA*SECY  
! --------------------------------------------
  REAL(qp) :: KL_E, KL_F, KL(M), & 
              KT_E, KT_F, KT(M), & 
	      HT_E, HT_F, HT(M), S(M), & 
	      LT_E, LT_F, LT(M)
! --------------------------------------------	      
!
  REAL(qp) :: P(0:M)
  LOGICAL FAIL  
!
  END MODULE DATA
!
!
! **********************************************************
! **********************************************************
! **********************************************************
!
!
 PROGRAM PM 
!
 USE DATA 
 IMPLICIT NONE
! 
 INTEGER I, J, K, L, I1, I2, OPT, MOPT, ILOAD, INDX(M)    
 CHARACTER*30 K_T_FILE, K_L_FILE, SPC_FILE, COE_FILE
 CHARACTER*10 CJUNK 
 CHARACTER*100 FFMMTT
 REAL(qp) COEFF(M+1)
!
 INTEGER ISO_MODES, ROT_MODES 
 COMPLEX*16 RROTS(MAX_MODES), A_RES(MAX_MODES), A_SEC, A_ELA
 CHARACTER*20 DATE, TIMC
 REAL(qp)   A_0_P2, A_E_P2 
 COMPLEX*16 A_0_P1, A_P1(M), A_P2(M-1)
 COMPLEX*16 RROOTS_P1(M), RROOTS_P2(M-1) 
 COMMON/PMT1/RROOTS_P1, A_0_P1,         A_P1
 COMMON/PMT2/RROOTS_P2, A_0_P2, A_E_P2, A_P2 
 COMPLEX*16 A_ELA_P, A_SEC_P
 COMPLEX*16 A_RES_PP(MAX_MODES), A_RES_P(MAX_MODES)  
 COMPLEX*16 CP(MAX_MODES,MAX_MODES)
 COMPLEX*16 GFACTOR(3), PMOT, DPMOT
 INTEGER, PARAMETER :: N_TIME_POINTS=6   
!
! **********************************************************
! **********************************************************
! **********************************************************
!  
!
 OPT=2 
!
!
 Write (*,*)"    - Main: - Reading the degree 2 Love numbers" 
! 
! ------------------------------------------------------------
! The Love numbers are imported from outside. They have been 
! computed by TBHP.F90 using the multi-precision FMLIB programs 
! to prevent numerical round off errors. 
! ------------------------------------------------------------
!
!
 OPEN(1,FILE='Deg_2_Love_numbers.dat',STATUS='UNKNOWN') 
!
! Reading the 10-lines header 
 	do i=1, 10 ; read(1,'(a10)')cjunk ; enddo
!
! Reading elastic loading k 
 	do i=1, 2  ; read(1,'(a10)')cjunk ; enddo
     	READ(1,*) KL_E 
!
! Reading fluid loading k 
 	do i=1, 2  ; read(1,'(a10)')cjunk ; enddo
     	READ(1,*) KL_F 
!
! Reading residues of loading k 
 	do i=1, 2  ; read(1,'(a10)')cjunk ; enddo
 	DO I=1, M 
     		READ(1,*) KL(I) 	
 	ENDDO
!
! Reading & re-arranging the coefficients of secular equation
 	do i=1, 2  ; read(1,'(a10)')cjunk ; enddo
 	DO I=1, M+1 
     		READ(1,*) COEFF(I) 		
 	ENDDO
 	DO I=0, M 
  		P(I)=COEFF(M+1-I) 
 	ENDDO
!
! Reading the s_j roots 
 	do i=1, 2  ; read(1,'(a10)')cjunk ; enddo
 	DO I=1, M 
     		READ(1,*) S(I) 	
 	ENDDO
!
! Reading elastic tidal h 
 	do i=1, 2  ; read(1,'(a10)')cjunk ; enddo
     	READ(1,*) HT_E 
!
! Reading fluid tidal h 
 	do i=1, 2  ; read(1,'(a10)')cjunk ; enddo
     	READ(1,*) HT_F 
!
! Reading residues of tidal h 
 	do i=1, 2  ; read(1,'(a10)')cjunk ; enddo
 	DO I=1, M 
     		READ(1,*) HT(I) 	
 	ENDDO
!
! ########################################################  Tidal loading (December 2012) 
!
! Reading elastic tidal l 
 	do i=1, 2  ; read(1,'(a10)')cjunk ; enddo
     	READ(1,*) LT_E 
!
! Reading fluid tidal h 
 	do i=1, 2  ; read(1,'(a10)')cjunk ; enddo
     	READ(1,*) LT_F 
!
! Reading residues of tidal h 
 	do i=1, 2  ; read(1,'(a10)')cjunk ; enddo
 	DO I=1, M 
     		READ(1,*) LT(I) 	
 	ENDDO
!
! ########################################################  Tidal loading (December 2012) 
!
!
! Reading elastic tidal k 
 	do i=1, 2  ; read(1,'(a10)')cjunk ; enddo
     	READ(1,*) KT_E 
!
! Reading fluid tidal k 
 	do i=1, 2  ; read(1,'(a10)')cjunk ; enddo
     	READ(1,*) KT_F 
!
! Reading residues of tidal k 
 	do i=1, 2  ; read(1,'(a10)')cjunk ; enddo
 	DO I=1, M 
     		READ(1,*) KT(I) 	
 	ENDDO
!
 CLOSE(1) 
!
!
!
!
! ------------------------------------------------------------
! Calling the Polar motion routine according to the option OPT. 
! After the call, the rotational roots and the residues are 
! stored in the main program for further computations. 9/22/09
! ------------------------------------------------------------
!
 IF    (OPT==1) THEN 
  	Write (*,*) "    - Main: - Calling POLAR_MOTION_1" 
 	ROT_MODES=M 
	ISO_MODES=M
 	CALL POLAR_MOTION_1 (mopt) 
	A_SEC=A_0_P1 
	A_ELA=DCMPLX(0d0,0d0)
		DO I=1, ROT_MODES 
			RROTS(I)=RROOTS_P1(I)
			A_RES(I)=A_P1(I) 
		ENDDO
 ELSEIF(OPT==2) THEN 
  	Write (*,*) "    - Main: - Calling POLAR_MOTION_2" 
 	ROT_MODES=M-1 
	ISO_MODES=M 
        CALL POLAR_MOTION_2 (mopt) 
	A_SEC=A_0_P2 
	A_ELA=A_E_P2 
		DO I=1, ROT_MODES 
			RROTS(I)=RROOTS_P2(I)
			A_RES(I)=A_P2(I) 
		ENDDO 
 ENDIF
!
!
!
! --------------------------------------------------------------
! Here we multiply the polar motion transfer function by 1 + k^L
! in the Laplace domain, decomposing the result in simple terms 
! --------------------------------------------------------------
!
 call DATE_AND_TIME (date,timc)
!
 Write(*,*) "    - pm_2: - Reporting the results on <<PM_time_domain.dat>> "
!
 open(3,file='PM_time_domain.dat',status='unknown')
!	
	Write(3,*) "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "	
        Write(3,*) "# "
  	Write(3,*) '# File <<PM_time_domain.dat>>, created by program PM_v6.f90 on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     '	    time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
	Write(3,*) "# COST GIA Test Suite "
	Write(3,*) "# Contributors: Giorgio Spada & Florence Colleoni"
	Write(3,*) "# Contact: giorgio.spada@gmail.com"
	Write(3,*) "# Number of isostatic/rotational modes: ", M, "/", M-1
	Write(3,*) "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "
	Write(3,*) "# "
!
! --- "Elastic term"
        Write(3,*) "  "
        A_ELA_P = A_ELA*(1D0+KL_E)
	WRITE(3,*) "# ELA          Re(A^_E)            Im(A^_E)"	
	Write(3,'(A4,1X,2(1X,E20.8))') 'E^ ', A_ELA_P 
!
! --- "Secular term"
        Write(3,*) "  "
        A_SEC_P = A_SEC*(1D0+KL_F)
	WRITE(3,*) "# SEC          Re(A^_S)            Im(A^_S)"	
	Write(3,'(A4,1X,2(1X,E20.8))') 'S^ ', A_SEC_P 
!
! --- CP array 
	DO I1=1, ROT_MODES
		DO I2=1, ISO_MODES
		CP (I1,I2) = A_RES(I1)*KL(I2)/(RROTS(I1)-S(I2)) 
		ENDDO
	ENDDO
!
! --- "Rotational term"
	DO J=1, ROT_MODES 
	A_RES_P(J)=DCMPLX(0D0,0D0)
		DO K=1, ISO_MODES 
		A_RES_P(J) = A_RES_P(J) + CP(J,K)
		ENDDO
		A_RES_P(J) = A_RES_P(J) + A_RES(J)*(1D0 + KL_E)
	ENDDO		
!
	WRITE(3,*) 
	WRITE(3,*) "# ROT          Re(a_i)              Im(a_i)              Re(A^_i)            Im(A^_i)"
  	DO J=1, ROT_MODES 
		WRITE(3,'(I4,1X,8(1X,E20.8))') J, DBLE(RROTS(J)), DIMAG(RROTS(J)),  &
		                                  DBLE(A_RES_P(J)),  DIMAG(A_RES_P(J))
	ENDDO
!
! --- "Isostatic term"
	DO J=1, ISO_MODES 
	A_RES_PP(J)=DCMPLX(0D0,0D0)
		DO K=1, ROT_MODES 
		A_RES_PP(J) = A_RES_PP(J) - CP(K,J)
		ENDDO
		A_RES_PP(J) = A_RES_PP(J) + A_ELA*KL(J) + A_SEC*(KL(J)/S(J))
	ENDDO	
!
	WRITE(3,*) 
	WRITE(3,*) "# ISO           s_i                Re(A^^_i)           Im(A^^_i)"
  	DO J=1, ISO_MODES 
		WRITE(3,'(I4,1X,8(1X,E20.8))') J, S(J),  & 
		                                  DBLE(A_RES_PP(J)),  DIMAG(A_RES_PP(J))
	ENDDO
!
  CLOSE(3) 
!
  END PROGRAM PM 
!
!
!
!
!
!
!#################################
 SUBROUTINE POLAR_MOTION_2  (mopt) 
!#################################
!
! Exclude the CW and set k^T_o = k^T_f 

 USE DATA
 Implicit NONE
!
 INTEGER I, J, K, L, MOPT
!
 REAL*8     ROOT_REAL(M-1), ROOT_IMAG(M-1), COEF_REAL(M), COEF_IMAG(M)
 REAL(qp) H(M-1),      B(0:M-1,M),   C(0:M-1)
 REAL(qp) NUME(0:M),   DENO(0:M-1),  NUME_P(0:M-2),     D_DENO(0:M-2)
 REAL(qp) A_0, A_E
 COMPLEX*16 NNN, DDD
 COMPLEX*16 ROOTS(M-1), A(M-1)
!
 CHARACTER*20 DATE, TIMC
 CHARACTER*20 FMT2, FMT3, FMT4
!
 COMMON/PMT2/ROOTS, A_0, A_E, A 
!
  Write (*,*)"    - pm_2: - Computing the 'h' constants"  
! Definition of the (real) amplitudes 'h' (h is REAL) 
  DO I=1, M 
      H(I) = -KT(I)/KT_F/S(I) 
  ENDDO
!
!
!
  Write (*,*)"    - pm_2: - Computing the 'b' constants"  
! Division of "p" by "s-s_i". The other possible recurrence 
! is unstable - (B is REAL)  
  DO I=1, M 
	B(M-1,I)=P(M) 
	DO K=2, M  
      	B(M-K,I)=P(M-K+1)+B(M-K+1,I)*S(I)
	ENDDO
  ENDDO
!
!
!
  Write (*,*)"    - pm_2: - Computing the 'c' constants"  
! Re-arrange the numerator after the above division. - 
! The 'c' are REAL. 
  DO K=0, M-1
 	C(K)=0D0 
 	DO I=1, M 
 	C(K)=C(K)+H(I)*B(K,I)
	ENDDO
  ENDDO
!
!
  Write (*,*)"    - pm_2: - Isolating polar motion"  
! We now turn the fraction upside down and redefine the polynomials...  
! Polynomials Nume & Deno are REAL 
  DO I=0, M
 	NUME(I)=P(I) 
  ENDDO
  DO I=0, M-1	
 	DENO(I)=C(I) 
  ENDDO
!
!
!
  Write(*,*) "    - pm_2: - Rotational roots"
! Prepares the coefficients for Subroutine CPOLY - these
! coefficients are REAL numbers.   
  DO I=1, M
          COEF_REAL(I)= DENO(M-I)
          COEF_IMAG(I)= 0d0
  ENDDO 
!
! Former call - CALL POLY2(taap,m-1,roots,hh,1000000,1)
!
! Current call to CPOLY (from http://www.netlib.org/toms/493)
! installed here with minor arrangements on september 19 2009
!  
  CALL CPOLY(COEF_REAL,COEF_IMAG,M-1,ROOT_REAL,ROOT_IMAG,FAIL)
!
  IF(     FAIL) Write(*,*) & 
          "    - pm_2: - Roots extraction (by CPOLY) UNSUCCESFULL ******* Warning"
  IF(.NOT.FAIL) Write(*,*) & 
          "    - pm_2: - Roots extraction (by CPOLY) SUCCESFUL"
!
! The Roots are COMPLEX
  DO I=1, M-1 
  	ROOTS(I)=DCMPLX(ROOT_REAL(I),ROOT_IMAG(I))
  ENDDO
!
!
!
  Write(*,*) "    - pm_2: - Decomposing the fraction"
! Decomposition of the fraction 
!
! Secular residue 
  A_0 = NUME(0)/DENO(0) 
! Write(*,*) 'Secular residue ', a_0 
!
! Elastic residue 
  A_E = NUME(M)/DENO(M-1) 
! Write(*,*) 'Elastic residue ', a_e 
!
! New numerator coefficients 
  DO K=0, M-2
  	NUME_P(K)=NUME(K+1) - A_E*DENO(K) - A_0*DENO(K+1) 
!       Write(*,*) 'New numerator coefficients ', k, nume_p(k) 	
  ENDDO 
!
!
  Write(*,*) "    - pm_2: - Derivative of the denominator"
! Denominator is s*'deno'. The derivative is thus a degree m polynomial
  DO K=0, M-2
 	D_DENO(K)=DENO(K+1)*FLOAT(K+1) 
  ENDDO
!
!
  Write(*,*) "    - pm_2: - Rotational modes"
! These are computed using the traditional formula for the inverse LT of the ratio 
! between a polynomial of degree 'm-2' and one of degree 'm-1' 
!
! Rotational modes 
  DO K=1, M-1 
  	NNN=0D0 
	DDD=0D0 
  	DO J=0, M-2 
  		NNN = NNN + NUME_P(J)*ROOTS(K)**J 
  		DDD = DDD + D_DENO(J)*ROOTS(K)**J 
 	 ENDDO	
  	A(K)=NNN/DDD
  ENDDO
!
 call DATE_AND_TIME (date,timc)
!
  Write(*,*) "    - pm_2: - Reporting the results on <<PM_spectrum.dat>> "
!
 open(3,file='PM_spectrum.dat',status='unknown')
!	
	Write(3,*) "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "	
  	Write(3,*) '# File <<PM_spectrum.dat>>, created by program PM_v6.f90 on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     '	    time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
	Write(3,*) "# COST GIA Test Suite "
	If(mopt==1) Write(3,*) "# Model: VS96"
	If(mopt==2) Write(3,*) "# Model: M03-L70-V01"
	Write(3,*) "# Contributors: Giorgio Spada & Florence Colleoni"
	Write(3,*) "# Contact: giorgio.spada@gmail.com"
	Write(3,*) "# Number of isostatic/rotational modes: ", M, "/", M-1
	Write(3,*) "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "
!
	FMT2='(I4,1X,8(1X,E20.8))'
	FMT3='(A4,1X,2(1X,E20.8))'
!
	WRITE(3,*) 
	WRITE(3,*) "# ELA           A_e"
  	WRITE(3,FMT3) 'E',  A_E 
!
	WRITE(3,*) 
	WRITE(3,*) "# SEC           A_s"
  	WRITE(3,FMT3) '0',  A_0 
!
	WRITE(3,*) 
	WRITE(3,*) "# ROT          Re(a_i)              Im(a_i)               Re(A_i)             Im(A_i)"
  	DO K=1, M-1 
		WRITE(3,FMT2) K, ROOT_REAL(K), ROOT_IMAG(K),  DBLE(A(K)),  DIMAG(A(K))
	ENDDO	
!
	WRITE(3,*) 
	WRITE(3,*) "# ISO           s(i)               Loading k             Tidal k"
  	DO K=1, M 
       		WRITE(3,FMT2) K, S(K), KL(K),  KT(K) 
  	ENDDO
!
  CLOSE(3) 
!
 END SUBROUTINE POLAR_MOTION_2
!
!
!
!
!
!
!##########################
 SUBROUTINE POLAR_MOTION_1  (mopt) 
!##########################
!
! Include the CW and set k^T_o = k^T_f 
!
 USE DATA
 IMPLICIT NONE
!
 CHARACTER*20 DATE, TIMC
 CHARACTER*20 FMT2, FMT3, FMT4
!
 INTEGER I, J, K, L, mopt 
 REAL*8     ROOT_REAL(M), ROOT_IMAG(M), COEF_REAL(M+1), COEF_IMAG(M+1) 
 REAL(qp)    R(M),         B(0:M-1,M),   C(0:M-1),       D_REAL(0:M), D_IMAG(0:M) 
 COMPLEX*16 NUME(0:M),    DENO(0:M),    D(0:M),         D_DENO(0:M)
 COMPLEX*16 TAAP(M+1),    ROOTS(M),     NNN,             DDD
 COMPLEX*16 A_0,          A_E, 		A(M)
!
 COMMON/PMT1/ROOTS, A_0, A 
!
  Write (*,*)"---> pm_1: - Computing the 'r' constants"  
! Definition of the (real) amplitudes 'r'
  DO I=1, M 
      R(I) = SR*KT(I)/KT_F/S(I) 
     	write(*,*) i, r(i)
  ENDDO
  read(*,*) 
!
!
!
  Write (*,*)"---> pm_1: - Computing the 'b' constants"  
! Division of "p" by "s-s_i". The other possible recurrence is unstable 
  DO I=1, M 
       B(M-1,I)=P(M) 
       DO K=2, M  
       B(M-K,I)=P(M-K+1)+B(M-K+1,I)*S(I)
       ENDDO
  ENDDO
!
!
!
  Write (*,*)"---> pm_1: - Computing the 'c' constants"  
! Re-arrange the numerator after the above division. - The 'c' are real. 
  DO K=0, M-1
 	C(K)=0D0 
 	DO I=1, M 
 	C(K)=C(K)+R(I)*B(K,I)
	ENDDO
  ENDDO
!
!
!
  Write (*,*)"---> pm_1: - Computing the 'd' constants"  
! Another re-organization of the numerator. The 'd' are *complex* since they come 
! from adding the (real) p polynomial to the (purely imaginary) i*c polynomial.   
!
  DO I=0, M 
 	D_REAL(I)=P(I)
  ENDDO
  DO I=0,M-1
 	D_IMAG(I)=C(I)  
  ENDDO
  D_IMAG(M)=0D0 
!
! Complex form of the numerator  
  DO I=0, M 
        D(I)=DCMPLX(D_REAL(I),D_IMAG(I))
  ENDDO 
! 
!
!
  Write (*,*)"---> pm_1: - Isolating polar motion"  
! We now turn upside down the fraction, and redefine the polynomials... the numerator 
! is real, the denominator is complex. 
  DO I=0, M
 	NUME(I)=P(I) 
 	DENO(I)=D(I) 
  ENDDO
!
! Updating the numerator... 
  DO I=0, M
 	NUME(I)= NUME(I)*DCMPLX(0.D0,-1.D0)*SR 
  ENDDO
!
!
!
  Write(*,*) "---> pm_1: - Rotational roots"
! We first rearrange the denominator then we find the rotational roots
! using rootfinder POLY2 (complex coefficients allowed in input). 
!
  DO I=1,M+1
	TAAP(I)=DENO(M+1-I) 
  ENDDO
!
! Prepares the coefficients for Subroutine CPOLY  
  DO I=1, M+1
          COEF_REAL(I)= DBLE (TAAP(I))
          COEF_IMAG(I)= DIMAG(TAAP(I))
  ENDDO 
!  
! Former call -  CALL POLY2(taap,m,roots,hh,10000,1)
!
! Current call to CPOLY (from http://www.netlib.org/toms/419)    --- 493
! installed here with minor arrangements on september 19 2009  
!
  CALL CPOLY(COEF_REAL,COEF_IMAG,M,ROOT_REAL,ROOT_IMAG,FAIL)
  IF(     FAIL) Write(*,*) & 
          "---> pm_1: - Roots extraction (by CPOLY) UNSUCCESFULL ******* Warning"
  IF(.NOT.FAIL) Write(*,*) & 
          "---> pm_1: - Roots extraction (by CPOLY) SUCCESFUL"
!
  DO I=1, M 
  	ROOTS(I)=DCMPLX(ROOT_REAL(I),ROOT_IMAG(I))
  ENDDO
!
!
!
  Write(*,*) "---> pm_1: - Derivative of the denominator"
! Denominator is s*'deno'. The derivative is thus a degree m polynomial
  DO K=0, M
 	D_DENO(K)=DENO(K)*FLOAT(K+1) 
  ENDDO
!
!
!
 Write(*,*) "---> pm_1: - Rotational modes"
! These are computed using the traditional formula for the inverse LT of the ratio 
! between a polynomial of degree 'm' and one of degree 'm+1' 
!
! Secular term 
  A_0=NUME(0)/D_DENO(0) 
! 
! Rotational modes 
  DO K=1, M 
  	NNN=0D0 
	DDD=0D0 
  	DO J=0, M 
  		NNN = NNN +  NUME(J)*ROOTS(K)**J 
  		DDD = DDD +D_DENO(J)*ROOTS(K)**J 
 	 ENDDO	
  	A(K)=NNN/DDD
  ENDDO
!
!
 	call DATE_AND_TIME (date,timc)
!
        open(3,file='PM_spectrum.dat',status='unknown')
!	
	Write(3,*) "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "	
  	Write(3,*) '# File <<PM_spectrum.dat>>, created by program PM_v6.f90 on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     '	    time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
	Write(3,*) "# COST GIA Test Suite "
	If(mopt==1) Write(3,*) "# Model: VS96"
	If(mopt==2) Write(3,*) "# Model: M03-L70-V01"
	Write(3,*) "# Contributors: Giorgio Spada & Florence Colleoni"
	Write(3,*) "# Contact: giorgio.spada@gmail.com"
	Write(3,*) "# Number of isostatic/rotational modes: ", M, "/", M
	Write(3,*) "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "
!
!	FMT2='(I4,6X,8(1X,E20.8))'
!	FMT3='(A4,6X,2(1X,E20.8))'
!

	FMT2='(I4,1X,8(1X,E20.8))'
	FMT3='(A4,1X,2(1X,E20.8))'

	A_E=0.
	WRITE(3,*) 
	WRITE(3,*) "# ELA         Re(A_e)              Im(A_E)"
  	WRITE(3,FMT3) 'E',        DBLE(A_E),DIMAG(A_E) 
!
	WRITE(3,*) 
	WRITE(3,*) "# SEC         Re(A_s)              Im(A_0)"
  	WRITE(3,FMT3) '0',        DBLE(A_0),DIMAG(A_0)
! 
	WRITE(3,*) 
	WRITE(3,*) "# ROT         Re(a_i)              Im(a_i)              Re(A_i)              Im(A_i)"
  	DO K=1, M 
		WRITE(3,FMT2) K, ROOT_REAL(K), ROOT_IMAG(K), DBLE(A(K)), DIMAG(A(K)) 
	ENDDO	
!
	WRITE(3,*) 
	WRITE(3,*) "# ISO          s(i)               Loading k             Tidal k"
  	DO K=1, M 
       		WRITE(3,FMT2) K, S(K), KL(K),  KT(K) 
  	ENDDO

!
!
END SUBROUTINE POLAR_MOTION_1 
!
!
!
     SUBROUTINE PRTC(N,P,Q)
     DOUBLE PRECISION P(50),Q(50)
     WRITE(6,10) (P(I),Q(I) ,I=1,N)
  10 FORMAT(//' COEFFICIENTS' /50(2D26.16/))
     RETURN
     END
!
     SUBROUTINE PRTZ(N,ZR,ZI)
     DOUBLE PRECISION ZR(50),ZI(50)
     WRITE(6,10) (ZR(I),ZI(I) ,I=1,N)
  10 FORMAT(//' ZEROS'/ 50(2D26.16/))
     RETURN
     END
!    
!
!
     SUBROUTINE CPOLY(OPR,OPI,DEGREE,ZEROR,ZEROI,FAIL)  	       !CPOL  10
!
! ********************************************************** 
! /// ALGORITHM 419 COLLECTED ALGORITHMS FROM ACM.       ///
! /// ALGORITHM APPEARED IN COMM. ACM, VOL. 15, NO. 02,  ///
! /// P. 097. 				                 ///
! ********************************************************** 
! FINDS THE ZEROS OF A COMPLEX POLYNOMIAL.
! OPR, OPI  -  DOUBLE PRECISION VECTORS OF REAL AND
! IMAGINARY PARTS OF THE COEFFICIENTS IN
! ORDER OF DECREASING POWERS.
! DEGREE    -  INTEGER DEGREE OF POLYNOMIAL.
! ZEROR, ZEROI  -  OUTPUT DOUBLE PRECISION VECTORS OF
! REAL AND IMAGINARY PARTS OF THE ZEROS.
! FAIL      -  OUTPUT LOGICAL PARAMETER,  TRUE  ONLY IF
! LEADING COEFFICIENT IS ZERO OR IF CPOLY
! HAS FOUND FEWER THAN DEGREE ZEROS.
! THE PROGRAM HAS BEEN WRITTEN TO REDUCE THE CHANCE OF OVERFLOW
! OCCURRING. IF IT DOES OCCUR, THERE IS STILL A POSSIBILITY THAT
! THE ZEROFINDER WILL WORK PROVIDED THE OVERFLOWED QUANTITY IS
! REPLACED BY A LARGE NUMBER.
! COMMON AREA
     COMMON/GLOBAL/PR,PI,HR,HI,QPR,QPI,QHR,QHI,SHR,SHI, SR,SI,TR,TI,PVR,PVI,ARE,MRE,ETA,INFIN,NN
     DOUBLE PRECISION SR,SI,TR,TI,PVR,PVI,ARE,MRE,ETA,INFIN,   &
 	 PR(50),PI(50),HR(50),HI(50),QPR(50),QPI(50),QHR(50),  &
 	 QHI(50),SHR(50),SHI(50)
! TO CHANGE THE SIZE OF POLYNOMIALS WHICH CAN BE SOLVED, REPLACE
! THE DIMENSION OF THE ARRAYS IN THE COMMON AREA.
     DOUBLE PRECISION XX,YY,COSR,SINR,SMALNO,BASE,XXX,ZR,ZI,BND,  &
 	 OPR(1),OPI(1),ZEROR(1),ZEROI(1), CMOD,SCALE,CAUCHY,DSQRT
     EXTERNAL SCALE
     LOGICAL FAIL,CONV
     INTEGER DEGREE,CNT1,CNT2
! INITIALIZATION OF CONSTANTS
     CALL MCON(ETA,INFIN,SMALNO,BASE)
     ARE = ETA
     MRE = 2.0D0*DSQRT(2.0D0)*ETA
     XX = .70710678
     YY = -XX
     COSR = -.060756474
     SINR = .99756405
     FAIL = .FALSE.
     NN = DEGREE+1
! ALGORITHM FAILS IF THE LEADING COEFFICIENT IS ZERO.
     IF (OPR(1) .NE. 0.0D0 .OR. OPI(1) .NE. 0.0D0) GO TO 10
 	 FAIL = .TRUE.
 	 RETURN
! REMOVE THE ZEROS AT THE ORIGIN IF ANY.
  10 IF (OPR(NN) .NE. 0.0D0 .OR. OPI(NN) .NE. 0.0D0) GO TO 20
 	 IDNN2 = DEGREE-NN+2
 	 ZEROR(IDNN2) = 0.0D0
 	 ZEROI(IDNN2) = 0.0D0
 	 NN = NN-1
 	 GO TO 10
! MAKE A COPY OF THE COEFFICIENTS.
  20 DO 30 I = 1,NN
 	 PR(I) = OPR(I)
 	 PI(I) = OPI(I)
 	 SHR(I) = CMOD(PR(I),PI(I))
  30 CONTINUE
! SCALE THE POLYNOMIAL.
     BND = SCALE (NN,SHR,ETA,INFIN,SMALNO,BASE)
     IF (BND .EQ. 1.0D0) GO TO 40
     DO 35 I = 1,NN
 	 PR(I) = BND*PR(I)
 	 PI(I) = BND*PI(I)
  35 CONTINUE
! START THE ALGORITHM FOR ONE ZERO .
  40 IF (NN.GT. 2) GO TO 50
! CALCULATE THE FINAL ZERO AND RETURN.
 	 CALL CDIVID(-PR(2),-PI(2),PR(1),PI(1),ZEROR(DEGREE), ZEROI(DEGREE))
 	 RETURN
! CALCULATE BND, A LOWER BOUND ON THE MODULUS OF THE ZEROS.
  50 DO 60 I = 1,NN
 	 SHR(I) = CMOD(PR(I),PI(I))
  60 CONTINUE
     BND = CAUCHY(NN,SHR,SHI)
! OUTER LOOP TO CONTROL 2 MAJOR PASSES WITH DIFFERENT SEQUENCES
! OF SHIFTS.
     DO 100 CNT1 = 1,2
! FIRST STAGE CALCULATION, NO SHIFT.
 	 CALL NOSHFT(5)
! INNER LOOP TO SELECT A SHIFT.
 	 DO 90 CNT2 = 1,9
! SHIFT IS CHOSEN WITH MODULUS BND AND AMPLITUDE ROTATED BY
! 94 DEGREES FROM THE PREVIOUS SHIFT
 	      XXX = COSR*XX-SINR*YY
 	      YY = SINR*XX+COSR*YY
 	      XX = XXX
 	      SR = BND*XX
 	      SI = BND*YY
! SECOND STAGE CALCULATION, FIXED SHIFT.
 	      CALL FXSHFT(10*CNT2,ZR,ZI,CONV)
 	      IF (.NOT. CONV) GO TO 80
! THE SECOND STAGE JUMPS DIRECTLY TO THE THIRD STAGE ITERATION.
! IF SUCCESSFUL THE ZERO IS STORED AND THE POLYNOMIAL DEFLATED.
 		   IDNN2 = DEGREE-NN+2
 		   ZEROR(IDNN2) = ZR
 		   ZEROI(IDNN2) = ZI
 		   NN = NN-1
 		   DO 70 I = 1,NN
 			PR(I) = QPR(I)
 			PI(I) = QPI(I)
  70		   CONTINUE
 		   GO TO 40
  80	      CONTINUE
! IF THE ITERATION IS UNSUCCESSFUL ANOTHER SHIFT IS CHOSEN.
  90	 CONTINUE
! IF 9 SHIFTS FAIL, THE OUTER LOOP IS REPEATED WITH ANOTHER
! SEQUENCE OF SHIFTS.
 100 CONTINUE
! THE ZEROFINDER HAS FAILED ON TWO MAJOR PASSES.
! RETURN EMPTY HANDED.
     FAIL = .TRUE.
     RETURN
     END
!
     SUBROUTINE  NOSHFT(L1)					       !NOSH1130
! COMPUTES  THE DERIVATIVE  POLYNOMIAL AS THE INITIAL H
! POLYNOMIAL AND COMPUTES L1 NO-SHIFT H POLYNOMIALS.
! COMMON AREA
     COMMON/GLOBAL/PR,PI,HR,HI,QPR,QPI,QHR,QHI,SHR,SHI,       &
 	   SR,SI,TR,TI,PVR,PVI,ARE,MRE,ETA,INFIN,NN
     DOUBLE PRECISION SR,SI,TR,TI,PVR,PVI,ARE,MRE,ETA,INFIN,  &
 	   PR(50),PI(50),HR(50),HI(50),QPR(50),QPI(50),QHR(50), &
 	   QHI(50),SHR(50),SHI(50)
     DOUBLE PRECISION XNI,T1,T2,CMOD
     N = NN-1
     NM1 = N-1
     DO 10 I = 1,N
 	 XNI = NN-I
 	 HR(I) = XNI*PR(I)/FLOAT(N)
 	 HI(I) = XNI*PI(I)/FLOAT(N)
  10 CONTINUE
     DO 50 JJ = 1,L1
 	 IF (CMOD(HR(N),HI(N)) .LE. ETA*10.0D0*CMOD(PR(N),PI(N))) GO TO 30
 	 CALL CDIVID(-PR(NN),-PI(NN),HR(N),HI(N),TR,TI)
 	 DO 20 I = 1,NM1
 	      J = NN-I
 	      T1 = HR(J-1)
 	      T2 = HI(J-1)
 	      HR(J) = TR*T1-TI*T2+PR(J)
 	      HI(J) = TR*T2+TI*T1+PI(J)
  20	 CONTINUE
 	 HR(1) = PR(1)
 	 HI(1) = PI(1)
 	 GO TO 50
! IF THE CONSTANT TERM IS ESSENTIALLY ZERO, SHIFT H COEFFICIENTS.
  30	 DO 40 I = 1,NM1
 	      J = NN-I
 	      HR(J) = HR(J-1)
 	      HI(J) = HI(J-1)
  40	 CONTINUE
 	 HR(1) = 0.0D0
 	 HI(1) = 0.0D0
  50 CONTINUE
     RETURN
     END
!
     SUBROUTINE FXSHFT(L2,ZR,ZI,CONV)				       !FXSH1550
! COMPUTES L2 FIXED-SHIFT H POLYNOMIALS AND TESTS FOR
! CONVERGENCE.
! INITIATES A VARIABLE-SHIFT ITERATION AND RETURNS WITH THE
! APPROXIMATE ZERO IF SUCCESSFUL.
! L2 - LIMIT OF FIXED SHIFT STEPS
! ZR,ZI - APPROXIMATE ZERO IF CONV IS .TRUE.
! CONV  - LOGICAL INDICATING CONVERGENCE OF STAGE 3 ITERATION
! COMMON AREA
     COMMON/GLOBAL/PR,PI,HR,HI,QPR,QPI,QHR,QHI,SHR,SHI, 	 &
 	    SR,SI,TR,TI,PVR,PVI,ARE,MRE,ETA,INFIN,NN	   
     DOUBLE PRECISION SR,SI,TR,TI,PVR,PVI,ARE,MRE,ETA,INFIN,	 &
 	    PR(50),PI(50),HR(50),HI(50),QPR(50),QPI(50),QHR(50), &
 	    QHI(50),SHR(50),SHI(50)
     DOUBLE PRECISION ZR,ZI,OTR,OTI,SVSR,SVSI,CMOD
 	 LOGICAL CONV,TEST,PASD,BOOL
     N = NN-1
! EVALUATE P AT S.
     CALL POLYEV(NN,SR,SI,PR,PI,QPR,QPI,PVR,PVI)
     TEST = .TRUE.
     PASD = .FALSE.
! CALCULATE FIRST T = -P(S)/H(S).
     CALL CALCT(BOOL)
! MAIN LOOP FOR ONE SECOND STAGE STEP.
     DO 50 J = 1,L2
 	 OTR = TR
 	 OTI = TI
! COMPUTE NEXT H POLYNOMIAL AND NEW T.
 	 CALL NEXTH(BOOL)
 	 CALL CALCT(BOOL)
 	 ZR = SR+TR
 	 ZI = SI+TI
! TEST FOR CONVERGENCE UNLESS STAGE 3 HAS FAILED ONCE OR THIS
! IS THE LAST H POLYNOMIAL .
 	 IF ( BOOL .OR. .NOT. TEST .OR. J .EQ. L2) GO TO 50
 	 IF (CMOD(TR-OTR,TI-OTI) .GE. .5D0*CMOD(ZR,ZI)) GO TO 40
 	      IF (.NOT. PASD) GO TO 30
! THE WEAK CONVERGENCE TEST HAS BEEN PASSED TWICE, START THE
! THIRD STAGE ITERATION, AFTER SAVING THE CURRENT H POLYNOMIAL
! AND SHIFT.
 		   DO 10 I = 1,N
 			SHR(I) = HR(I)
 			SHI(I) = HI(I)
  10		   CONTINUE
 		   SVSR = SR
 		   SVSI = SI
 		   CALL VRSHFT(10,ZR,ZI,CONV)
 		   IF (CONV) RETURN
! THE ITERATION FAILED TO CONVERGE. TURN OFF TESTING AND RESTORE
! H,S,PV AND T.
 		   TEST = .FALSE.
 		   DO 20 I = 1,N
 			HR(I) = SHR(I)
 			HI(I) = SHI(I)
  20		   CONTINUE
 		   SR = SVSR
 		   SI = SVSI
 		   CALL POLYEV(NN,SR,SI,PR,PI,QPR,QPI,PVR,PVI)
 		   CALL CALCT(BOOL)
 		   GO TO 50
  30	      PASD = .TRUE.
 	      GO TO 50
  40	 PASD = .FALSE.
  50 CONTINUE
! ATTEMPT AN ITERATION WITH FINAL H POLYNOMIAL FROM SECOND STAGE.
     CALL VRSHFT(10,ZR,ZI,CONV)
     RETURN
     END
!
     SUBROUTINE VRSHFT(L3,ZR,ZI,CONV)				       !VRSH2230
! CARRIES OUT THE THIRD STAGE ITERATION.
! L3 - LIMIT OF STEPS IN STAGE 3.
! ZR,ZI   - ON ENTRY CONTAINS THE INITIAL ITERATE, IF THE
! ITERATION CONVERGES IT CONTAINS THE FINAL ITERATE
! ON EXIT.
! CONV    -  .TRUE. IF ITERATION CONVERGES
! COMMON AREA
     COMMON/GLOBAL/PR,PI,HR,HI,QPR,QPI,QHR,QHI,SHR,SHI, 	  &
 	   SR,SI,TR,TI,PVR,PVI,ARE,MRE,ETA,INFIN,NN
     DOUBLE PRECISION SR,SI,TR,TI,PVR,PVI,ARE,MRE,ETA,INFIN,	  &
 	   PR(50),PI(50),HR(50),HI(50),QPR(50),QPI(50),QHR(50),   &
 	   QHI(50),SHR(50),SHI(50)
     DOUBLE PRECISION ZR,ZI,MP,MS,OMP,RELSTP,R1,R2,CMOD,DSQRT,ERREV,TP
     LOGICAL CONV,B,BOOL
     CONV = .FALSE.
     B = .FALSE.
     SR = ZR
     SI = ZI
! MAIN LOOP FOR STAGE THREE
     DO 60 I = 1,L3
! EVALUATE P AT S AND TEST FOR CONVERGENCE.
 	 CALL POLYEV(NN,SR,SI,PR,PI,QPR,QPI,PVR,PVI)
 	 MP = CMOD(PVR,PVI)
 	 MS = CMOD(SR,SI)
 	 IF (MP .GT. 20.0D0*ERREV(NN,QPR,QPI,MS,MP,ARE,MRE)) GO TO 10
! POLYNOMIAL VALUE IS SMALLER IN VALUE THAN A BOUND ON THE ERROR
! IN EVALUATING P, TERMINATE THE ITERATION.
 	      CONV = .TRUE.
 	      ZR = SR
 	      ZI = SI
 	      RETURN
  10	 IF (I .EQ. 1) GO TO 40
 	      IF (B .OR. MP .LT.OMP .OR. RELSTP .GE. .05D0) GO TO 30
! ITERATION HAS STALLED. PROBABLY A CLUSTER OF ZEROS. DO 5 FIXED
! SHIFT STEPS INTO THE CLUSTER TO FORCE ONE ZERO TO DOMINATE.
 		   TP = RELSTP
 		   B = .TRUE.
 		   IF (RELSTP .LT. ETA) TP = ETA
 		   R1 = DSQRT(TP)
 		   R2 = SR*(1.0D0+R1)-SI*R1
 		   SI = SR*R1+SI*(1.0D0+R1)
 		   SR = R2
 		   CALL POLYEV(NN,SR,SI,PR,PI,QPR,QPI,PVR,PVI)
 		   DO 20 J = 1,5
 			CALL CALCT(BOOL)
 			CALL NEXTH(BOOL)
  20		   CONTINUE
     OMP = INFIN
 		   GO TO 50
! EXIT IF POLYNOMIAL VALUE INCREASES SIGNIFICANTLY.
  30	      IF (MP*.1D0 .GT. OMP) RETURN
  40	 OMP = MP
! CALCULATE NEXT ITERATE.
  50	 CALL CALCT(BOOL)
 	 CALL NEXTH(BOOL)
 	 CALL CALCT(BOOL)
 	 IF (BOOL) GO TO 60
 	 RELSTP = CMOD(TR,TI)/CMOD(SR,SI)
 	 SR = SR+TR
 	 SI = SI+TI
  60 CONTINUE
     RETURN
     END
!
     SUBROUTINE CALCT(BOOL)					       !CALC2890
! COMPUTES  T = -P(S)/H(S).
! BOOL   - LOGICAL, SET TRUE IF H(S) IS ESSENTIALLY ZERO.
! COMMON AREA
     COMMON/GLOBAL/PR,PI,HR,HI,QPR,QPI,QHR,QHI,SHR,SHI, 	  & 
 	    SR,SI,TR,TI,PVR,PVI,ARE,MRE,ETA,INFIN,NN
     DOUBLE PRECISION SR,SI,TR,TI,PVR,PVI,ARE,MRE,ETA,INFIN,	  &
 	    PR(50),PI(50),HR(50),HI(50),QPR(50),QPI(50),QHR(50),  &
 	    QHI(50),SHR(50),SHI(50)
     DOUBLE PRECISION HVR,HVI,CMOD
     LOGICAL BOOL
     N = NN-1
! EVALUATE H(S).
     CALL POLYEV(N,SR,SI,HR,HI,QHR,QHI,HVR,HVI)
     BOOL = CMOD(HVR,HVI) .LE. ARE*10.0D0*CMOD(HR(N),HI(N))
     IF (BOOL) GO TO 10
 	 CALL CDIVID(-PVR,-PVI,HVR,HVI,TR,TI)
 	 RETURN
  10 TR = 0.0D0
     TI = 0.0D0
     RETURN
     END
!
     SUBROUTINE NEXTH(BOOL)					       !NEXT3110
! CALCULATES THE NEXT SHIFTED H POLYNOMIAL.
! BOOL   -  LOGICAL, IF .TRUE. H(S) IS ESSENTIALLY ZERO
! COMMON AREA
     COMMON/GLOBAL/PR,PI,HR,HI,QPR,QPI,QHR,QHI,SHR,SHI, 	&
 	   SR,SI,TR,TI,PVR,PVI,ARE,MRE,ETA,INFIN,NN
     DOUBLE PRECISION SR,SI,TR,TI,PVR,PVI,ARE,MRE,ETA,INFIN,	&
 	   PR(50),PI(50),HR(50),HI(50),QPR(50),QPI(50),QHR(50), &
 	   QHI(50),SHR(50),SHI(50)
     DOUBLE PRECISION T1,T2
     LOGICAL BOOL
     N = NN-1
     NM1 = N-1
     IF (BOOL) GO TO 20
 	 DO 10 J = 2,N
 	      T1 = QHR(J-1)
 	      T2 = QHI(J-1)
 	      HR(J) = TR*T1-TI*T2+QPR(J)
 	      HI(J) = TR*T2+TI*T1+QPI(J)
  10	 CONTINUE
 	 HR(1) = QPR(1)
 	 HI(1) = QPI(1)
 	 RETURN
! IF H(S) IS ZERO REPLACE H WITH QH.
  20 DO 30 J = 2,N
 	 HR(J) = QHR(J-1)
 	 HI(J) = QHI(J-1)
  30 CONTINUE
     HR(1) = 0.0D0
     HI(1) = 0.0D0
     RETURN
     END
!
     SUBROUTINE POLYEV(NN,SR,SI,PR,PI,QR,QI,PVR,PVI)		       !POLY3430
! EVALUATES A POLYNOMIAL  P  AT  S  BY THE HORNER RECURRENCE
! PLACING THE PARTIAL SUMS IN Q AND THE COMPUTED VALUE IN PV.
     DOUBLE PRECISION PR(NN),PI(NN),QR(NN),QI(NN), &
 		      SR,SI,PVR,PVI,T
     QR(1) = PR(1)
     QI(1) = PI(1)
     PVR = QR(1)
     PVI = QI(1)
     DO 10 I = 2,NN
 	 T = PVR*SR-PVI*SI+PR(I)
 	 PVI = PVR*SI+PVI*SR+PI(I)
 	 PVR = T
 	 QR(I) = PVR
 	 QI(I) = PVI
  10 CONTINUE
     RETURN
     END
!
     DOUBLE PRECISION FUNCTION ERREV(NN,QR,QI,MS,MP,ARE,MRE)	       !ERRE3610
! BOUNDS THE ERROR IN EVALUATING THE POLYNOMIAL BY THE HORNER
! RECURRENCE.
! QR,QI - THE PARTIAL SUMS
! MS	-MODULUS OF THE POINT
! MP	-MODULUS OF POLYNOMIAL VALUE
! ARE, MRE -ERROR BOUNDS ON COMPLEX ADDITION AND MULTIPLICATION
     DOUBLE PRECISION QR(NN),QI(NN),MS,MP,ARE,MRE,E,CMOD
     E = CMOD(QR(1),QI(1))*MRE/(ARE+MRE)
     DO 10 I = 1,NN
 	 E = E*MS+CMOD(QR(I),QI(I))
  10 CONTINUE
     ERREV = E*(ARE+MRE)-MP*MRE
     RETURN
     END
     DOUBLE PRECISION FUNCTION CAUCHY(NN,PT,Q)  		       !CAUC3760
! CAUCHY COMPUTES A LOWER BOUND ON THE MODULI OF THE ZEROS OF A
! POLYNOMIAL - PT IS THE MODULUS OF THE COEFFICIENTS.
     DOUBLE PRECISION Q(NN),PT(NN),X,XM,F,DX,DF, DABS,DEXP,DLOG
     PT(NN) = -PT(NN)
! COMPUTE UPPER ESTIMATE OF BOUND.
     N = NN-1
     X = DEXP( (DLOG(-PT(NN)) - DLOG(PT(1)))/FLOAT(N) )
     IF (PT(N).EQ.0.0D0) GO TO 20
! IF NEWTON STEP AT THE ORIGIN IS BETTER, USE IT.
 	 XM = -PT(NN)/PT(N)
 	 IF (XM.LT.X) X=XM
! CHOP THE INTERVAL (0,X) UNITL F LE 0.
  20 XM = X*.1D0
     F = PT(1)
     DO 30 I = 2,NN
 	 F = F*XM+PT(I)
  30 CONTINUE
     IF (F.LE. 0.0D0) GO TO 40
 	 X = XM
 	 GO TO 20
  40 DX = X
! DO NEWTON ITERATION UNTIL X CONVERGES TO TWO DECIMAL PLACES.
  50 IF (DABS(DX/X) .LE. .005D0) GO TO 70
 	 Q(1) = PT(1)
 	 DO 60 I = 2,NN
 	      Q(I) = Q(I-1)*X+PT(I)
  60	 CONTINUE
 	 F = Q(NN)
 	 DF = Q(1)
 	 DO 65 I = 2,N
 	      DF = DF*X+Q(I)
  65	 CONTINUE
 	 DX = F/DF
 	 X = X-DX
 	 GO TO 50
  70 CAUCHY = X
     RETURN
     END
!       
     DOUBLE PRECISION FUNCTION SCALE(NN,PT,ETA,INFIN,SMALNO,BASE)     
! RETURNS A SCALE FACTOR TO MULTIPLY THE COEFFICIENTS OF THE
! POLYNOMIAL. THE SCALING IS DONE TO AVOID OVERFLOW AND TO AVOID
! UNDETECTED UNDERFLOW INTERFERING WITH THE CONVERGENCE
! CRITERION.  THE FACTOR IS A POWER OF THE BASE.
! PT - MODULUS OF COEFFICIENTS OF P
! ETA,INFIN,SMALNO,BASE - CONSTANTS DESCRIBING THE
! FLOATING POINT ARITHMETIC.
     DOUBLE PRECISION PT(NN),ETA,INFIN,SMALNO,BASE,HI,LO,  &
 		      MAX,MIN,X,SC,DSQURT,DLOG
! FIND LARGEST AND SMALLEST MODULI OF COEFFICIENTS.
     HI = DSQRT(INFIN)
     LO = SMALNO/ETA
     MAX = 0.0D0
     MIN = INFIN
     DO 10 I = 1,NN
 	 X = PT(I)
 	 IF (X .GT. MAX) MAX = X
 	 IF (X .NE. 0.0D0 .AND. X.LT.MIN) MIN = X
  10 CONTINUE
! SCALE ONLY IF THERE ARE VERY LARGE OR VERY SMALL COMPONENTS.
     SCALE = 1.0D0
     IF (MIN .GE. LO .AND. MAX .LE. HI) RETURN
     X = LO/MIN
     IF (X .GT. 1.0D0) GO TO 20
 	 SC = 1.0D0/(DSQRT(MAX)*DSQRT(MIN))
 	 GO TO 30
  20 SC = X
     IF (INFIN/SC .GT. MAX) SC = 1.0D0
  30 L = DLOG(SC)/DLOG(BASE) + .500
     SCALE = BASE**L
     RETURN
     END     
!               
     SUBROUTINE CDIVID(AR,AI,BR,BI,CR,CI)			       !CDIV4490
! COMPLEX DIVISION C = A/B, AVOIDING OVERFLOW.
     DOUBLE PRECISION AR,AI,BR,BI,CR,CI,R,D,T,INFIN,DABS
     IF (BR .NE. 0.0D0  .OR. BI .NE. 0.0D0) GO TO 10
! DIVISION BY ZERO, C = INFINITY.
 	 CALL MCON (T,INFIN,T,T)
 	 CR = INFIN
 	 CI = INFIN
 	 RETURN
  10 IF (DABS(BR) .GE. DABS(BI)) GO TO 20
 	 R = BR/BI
 	 D = BI+R*BR
 	 CR = (AR*R+AI)/D
 	 CI = (AI*R-AR)/D
 	 RETURN
  20 R = BI/BR
     D = BR+R*BI
     CR = (AR+AI*R)/D
     CI = (AI-AR*R)/D
     RETURN
     END
     DOUBLE PRECISION FUNCTION CMOD(R,I)			       !CMOD4700
! MODULUS OF A COMPLEX NUMBER AVOIDING OVERFLOW.
     DOUBLE PRECISION R,I,AR,AI,DABS,DSQURT
     AR = DABS(R)
     AI = DABS(I)
     IF (AR .GE. AI) GO TO 10
 	 CMOD = AI*DSQRT(1.0D0+(AR/AI)**2)
 	 RETURN
  10 IF (AR .LE. AI) GO TO 20
 	 CMOD = AR*DSQRT(1.0D0+(AI/AR)**2)
 	 RETURN
  20 CMOD = AR*DSQRT(2.0D0)
     RETURN
     END
!
     SUBROUTINE MCON(ETA,INFINY,SMALNO,BASE)			       !MCON4840
! MCON PROVIDES MACHINE CONSTANTS USED IN VARIOUS PARTS OF THE
! PROGRAM. THE USER MAY EITHER SET THEM DIRECTLY OR USE THE
! STATEMENTS BELOW TO COMPUTE THEM. THE MEANING OF THE FOUR
! CONSTANTS ARE -
! ETA	    THE MAXIMUM RELATIVE REPRESENTATION ERROR
! WHICH CAN BE DESCRIBED AS THE SMALLEST POSITIVE
! FLOATING-POINT NUMBER SUCH THAT 1.0D0 + ETA IS
! GREATER THAN 1.0D0.
! INFINY    THE LARGEST FLOATING-POINT NUMBER
! SMALNO    THE SMALLEST POSITIVE FLOATING-POINT NUMBER
! BASE      THE BASE OF THE FLOATING-POINT NUMBER SYSTEM USED
! LET T BE THE NUMBER OF BASE-DIGITS IN EACH FLOATING-POINT
! NUMBER(DOUBLE PRECISION). THEN ETA IS EITHER .5*B**(1-T)
! OR B**(1-T) DEPENDING ON WHETHER ROUNDING OR TRUNCATION
! IS USED.
! LET M BE THE LARGEST EXPONENT AND N THE SMALLEST EXPONENT
! IN THE NUMBER SYSTEM. THEN INFINY IS (1-BASE**(-T))*BASE**M
! AND SMALNO IS BASE**N.
! THE VALUES FOR BASE,T,M,N BELOW CORRESPOND TO THE IBM/360.
     DOUBLE PRECISION ETA,INFINY,SMALNO,BASE
     INTEGER M,N,T
     BASE = 16.0D0
     T = 14
     M = 63
     N = -65
     ETA = BASE**(1-T)
     INFINY = BASE*(1.0D0-BASE**(-T))*BASE**(M-1)
     SMALNO = (BASE**(N+3))/BASE**3
     RETURN
     END
!
!
!
!
!
!
