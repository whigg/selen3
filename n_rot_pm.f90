!
! This is program "N_ROT_PM.F90" 
!
! Computes the contribution to sea level change "S" due to rotational variations.
! This program is "based" on "PMD_MOD.F90", with some difference in notation (!) 
!
! WARNING: For the moment assumes a reference surface gravity value from "data.inc"
!
! GS Sept 3 2010 GS for version 3.2 of SELEN 
! Re-touched November 2010 GS
! GS Prague November 2011: Missing rotational term in U and N
!                          After some tests with IDA and discussion with Zdeneck Martinek 
! March 2012: Implementation of the numerical derivative "in the future"
! Adapted to the ELASTIC multi-step case by GS on March 04, 2012 
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
! Input files:... ... ... ... 
!	- ... ... ...  
!
! Output files: ... ... ... ... 
!	- ... ... ...  
!
!
!
  MODULE DATA 
  IMPLICIT NONE 
  INCLUDE "data.inc"
!
  REAL*8 HTE   ! h tidal elastic 
  REAL*8 KTE   ! k tidal elastic 
  REAL*8 KTS   ! k tidal secular 
  REAL*8 HLE   ! h loading elastic 
  REAL*8 KLE   ! k loading elastic 
  REAL*8 AEP 
  REAL*8 BEP
!
  REAL*8, PARAMETER :: CCCC=8.0394d37, & 
  		       AAAA=8.0131d37, & 
		       OMEGA=7.292115d-5
!
  END MODULE DATA
!
!
!
!
!
!
!INCLUDE "harmonics.f90"
!
 PROGRAM SROTAZ
!
 USE DATA
 IMPLICIT NONE 
!
 CHARACTER*22, PARAMETER :: FILE_LOAD='shload_re.bin'   
! CHARACTER*22, PARAMETER :: FILE_SROT='srotaz_21.dat'
! CHARACTER*22, PARAMETER :: FILE_UROT='urotaz_21.dat'
 CHARACTER*22, PARAMETER :: FILE_NROT='nrotaz_21.dat'
!
 REAL*8, PARAMETER :: ERAD=6.371D6
 REAL*8, PARAMETER :: SIGMA_21=OMEGA**2*ERAD**2/SQRT(30D0)
!
 INTEGER I, J, K, IP, JP, KP, J21, J_INDEX
 COMPLEX*16 SFUNC
!
! ====== Modified March 2012 GS: N+1
!     
 COMPLEX*16 LOT  (JMAX,0:NN+1), DLOAD(JMAX,0:NN+1)
 COMPLEX*16 NROT (JMAX,0:NN+1), DINER_R(0:NN+1)
!
 CHARACTER*20 DATE, TIMC
! 
!
!
!
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ ! 
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !
!					 !
! Inertia variations for a rigid Earth   !
!				         !
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ ! 
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !
!
!
!
 Write(*,*) '    - Reading the harmonics of *total* surface load from file ', & 
                   trim(adjustl(FILE_LOAD)) 
!
! File with load information 
 open (1,file=FILE_LOAD,status='unknown',form='unformatted') 
 read (1) LOT
 close(1) 
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
  HLE=-0.9915810331D+0 
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
  BEP=(1D0+KTE-0D0)*AEP
!
!
! --- The following loops are unaltered 
!>>>> Defining the load variation of degree j at time k 
!
 dload(:,:)=(0.d0,0.d0)
!
 do 1 j=1, jmax 
 do 1 k=0, nn+1 	 ! - Updated 3/12
     if(k==0) then 
         dload(j,k) = lot(j,k)
	      else 
         dload(j,k) = lot(j,k) - lot(j,k-1)
     endif	 
1 continue
!
!! --- The following loops are unaltered 
!
!>>>> Time history of the "rigid" inertia, only "due to loading on a RIGID Earth" 
!     diner_r(k) are coefficients of a development in series of Heaviside functions
!
  do 4 k=0, nn+1        !  - Updated 3/12
!
! This rigid inertia is a complex variable, related with 
! polar motion. Only the degree 2 and order 1 components 
! of the load variation are involved -See my theory notes.
!
 diner_r(k) = (4d0*pi/3d0)*sqrt(6d0/5d0)*conjg(dload(j_index(2,1),k))
!
4 continue
!
! 
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!    Degree 2 and order 1 term of rotational sea level 
! ///////////////////////////////////////////////////////
!
 open (1,file=FILE_NROT,status='unknown',form='formatted') 
!
     NROT (:,:)  = (0D0,0D0)
!
     J21=J_INDEX(2,1)
!
     DO 15 I=0, NN +1 ! - Updated 3/12
!     
     NROT (J21,i)= (0D0,0D0)
!     
          DO 16 K=0, I 
!     
     	     NROT (J21,I)  = NROT (J21,I)  + DCONJG(DINER_R(K)) * SFUNC(DELTA*DFLOAT(I-K))
! 
 16       CONTINUE
!     
             NROT (J21,I) =  NROT (J21,I) * (SIGMA_21/GRA_REF)*(ERAD**4/(CCCC-AAAA))
!     
 15   CONTINUE 
!
      write(1,*)NROT
!
 close(1) 
!
!

!
! Work in progress! Work in progress! Work in progress! Work in progress! 
! Work in progress! Work in progress! Work in progress! Work in progress! 
! Work in progress! Work in progress! Work in progress! Work in progress! 
! Work in progress! Work in progress! Work in progress! Work in progress! 
!
!
 END PROGRAM SROTAZ
!
!
!
!
!
!
 FUNCTION SFUNC(X)
 USE DATA
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
! Updated for the elastic-multi step case by GS March 2012 
! Only the elastic term remains...
! ---------------------------------------------------------------
!
  SFUNC = BEP
!
 END FUNCTION SFUNC 
!
!
!
!
!
!
