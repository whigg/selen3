!
! This is program "S_ROT.F90" 
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
  INTEGER, PARAMETER :: M=4*NV   ! Number of isostatic modes 
  INTEGER, PARAMETER :: MP=M-1    ! Number of rotational modes 
! 
  REAL*8  KLE                     ! Elastic degree 2 k Love number (loading) 
  REAL*8  KLF                       ! Fluid degree 2 k Love number (loading) 
  REAL*8  KL(M)               ! Residues of degree 2 k Love number (loading)   
  REAL*8  KTE                     ! Elastic degree 2 k Love number (tidal) 
  REAL*8  KTF                       ! Fluid degree 2 k Love number (tidal) 
  REAL*8  KT(M)               ! Residues of degree 2 k Love number (tidal)   
  REAL*8  HTE                     ! Elastic degree 2 h Love number (tidal) 
  REAL*8  HTF                       ! Fluid degree 2 h Love number (tidal) 
  REAL*8  HT(M)               ! Residues of degree 2 h Love number (tidal)   
!
  REAL*8  S_ROOT(M)                  ! Degree 2 isostatic roots  
  REAL*8  AEP                        ! (Primed) elastic residue of the PMTF 
  REAL*8  ASP                        ! (Primed) secular residue of the PMTF 
  COMPLEX*16 A_ROOT(MP)              ! Rotational roots   
  COMPLEX*16 AP(MP)                  ! (Primed) rotational residues 
!
  REAL*8 BEP                     ! (Primed) elastic residue of S_rot
  REAL*8 BSP		     ! (Primed) secular residue of S_rot
  COMPLEX*16 BP (MP)		     ! (Primed) residue of S_rot  
  COMPLEX*16 BPP(M)  		     ! (Doubly primed) secular residue of S_rot
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
 CHARACTER*22, PARAMETER :: FILE_SROT='srotaz_21.dat'
 CHARACTER*22, PARAMETER :: FILE_UROT='urotaz_21.dat'
 CHARACTER*22, PARAMETER :: FILE_NROT='nrotaz_21.dat'
!
 REAL*8, PARAMETER :: ERAD=6.371D6
 REAL*8, PARAMETER :: SIGMA_21=OMEGA**2*ERAD**2/SQRT(30D0)
!
 REAL*8 GTE, GTV(M) 
 INTEGER I, J, K, IP, JP, KP, J21, J_INDEX
 COMPLEX*16 SFUNC, P(MP,M), E(M), F(MP)  
!
! ====== Modified March 2012 GS: N+1
!     
 COMPLEX*16 LOT  (JMAX,0:NN+1), DLOAD(JMAX,0:NN+1)
 COMPLEX*16 SROT (JMAX,0:NN+1), DINER_R(0:NN+1)
!
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
 CALL READ_DATA 
!
!
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
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ ! 
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !
!				       	        !
! Computing the products in the Laplace domain  !
!				                !
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ ! 
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !
!
!
! --- Elastic component of the Sea Level Green function
  GTE = 1D0 + KTE - HTE
!write(*,*) "GTE ", kte, hte, gte
!
! --- Viscous components of the Sea Level Green function
  DO  IP=1, M
     GTV(IP) = KT(IP) - HT(IP) 
!write(*,*) "GTV ", ip, gtv(ip)
  ENDDO
!
! --- Primed elastic term 
  BEP = GTE*AEP 
!write(*,*) "BEP ", bep
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
!write(*,*) ' bp ', i, bp(i) 
  ENDDO
!
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
!write(*,*) ' bpp ', ip, bpp(ip) 
  ENDDO
!
! 
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!    Degree 2 and order 1 term of rotational sea level 
! ///////////////////////////////////////////////////////
!
 open (1,file=FILE_SROT,status='unknown',form='formatted') 
!
     SROT (:,:)  = (0D0,0D0)
!
     J21=J_INDEX(2,1)
!
     DO 15 I=0, NN +1 ! - Updated 3/12
!     
     SROT (J21,i)= (0D0,0D0)
!     
          DO 16 K=0, I 
!     
     	     SROT (J21,I)  = SROT (J21,I)  + DCONJG(DINER_R(K)) * SFUNC(DELTA*DFLOAT(I-K))
! 
 16       CONTINUE
!     
        SROT (J21,I) =  SROT (J21,I) * (SIGMA_21/GRA_REF)*(ERAD**4/(CCCC-AAAA))
!     
 15   CONTINUE 
!
      write(1,*)SROT
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
  SUBROUTINE READ_DATA 
  USE DATA
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
  END SUBROUTINE READ_DATA 
!
!
!
!
!
!
