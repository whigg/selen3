!
! This is program "PMD.F90" 
!
! GS July 29 2010 for version 3.2 of SELEN 
! *** Revised GS August 2010 - g95 - Double precision implementation
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
! Input files:
!	- shload.bin 
!
!
!

  MODULE DATA 
  IMPLICIT NONE 
  INCLUDE "data.inc"
!
  INTEGER, PARAMETER :: M=4*NV   ! Number of isostatic modes 
  INTEGER, PARAMETER :: MP=M-1   ! Number of rotational modes 
  REAL*8  KLE                    ! Elastic degree 2 k Love number 
  REAL*8  KLF                    ! Fluid degree 2 k Love number  
  REAL*8  KL(M)                  ! Residues of degree 2 k Love number    
  REAL*8  S(M)                   ! Degree 2 isostatic poles 
  REAL*8  R(M)
  REAL*8  AE                     ! Elastic residue of the PMTF 
  REAL*8  AS                     ! Secular residue of the PMTF 
  COMPLEX*16 A(MP)               ! Rotational roots   
  COMPLEX*16 AA(MP)              ! Rotational residues 
  COMPLEX*16 C(1:MP,1:M)
  COMPLEX*16 P(1:M), Q(1:MP)
!
  REAL*8, PARAMETER :: CCCC=8.0394d37, AAAA=8.0131d37, OMEGA=7.292115d-5
  REAL*8, PARAMETER :: SECY=365.25d0*24Q0*3600Q0*1000Q0
  REAL*8, PARAMETER :: SR=((CCCC-AAAA)/AAAA)*OMEGA*SECY  
!
  END MODULE DATA








!INCLUDE "harmonics.f90"
 PROGRAM PMD 
 USE DATA
 IMPLICIT NONE 

 INTEGER I, J, K, JP, J_INDEX
 COMPLEX*16 LOT  (JMAX,0:NN)
 COMPLEX*16 DLOAD(JMAX,0:NN)
 COMPLEX*16 DINER_R(0:NN), DINER_D(0:NN)
! 
 CHARACTER*22, PARAMETER :: FILELOAD='shload.bin'  
 
 REAL*8, PARAMETER :: ERAD=6.371D6
 REAL*8, PARAMETER :: RAD2DEG=180D0/PI
!
 REAL*8 CNV(0:NN,0:NN)
!
 COMPLEX*16 MFUNC, DMFUNC
!
 COMPLEX*16 IIII(JMAX,0:NN)
!
 COMPLEX*16 PM(0:NN), DPM(0:NN)
!
 REAL*8 M1(0:NN), M2(0:NN), M1DOT(0:NN), M2DOT(0:NN)

 REAL*8 POLAR_DISP, POLAR_RATE
!
!
!
Write(*,*) '    - Retrieving information from existing SELEN data:' 
!
!
Write(*,*) '    - Reading the harmonics of *total* surface load from ', trim(adjustl(FILELOAD)) 
!
 open (1,file=FILELOAD,status='unknown',form='unformatted') 
 read (1) LOT
 close(1) 
!
! 
     CALL READ_DATA 



!


!
!>>>> Defining the load variation of degree j at time k 
!
 dload(:,:)=(0d0,0d0)
!
 do 1 j=1, jmax 
 do 1 k=0, nn 	  
     if(k==0) then 
         dload(j,k) = lot(j,k)
	      else 
         dload(j,k) = lot(j,k) - lot(j,k-1)
     endif	 
1 continue

 		Write(*,*) "    - Reading array 'I'"
		open(1,file='shice.dat',status='unknown') 
		read(1,*) IIII
		close(1)


 j=j_index(2,1)
 do k=0, nn 	  

 write(*,*) "ice ", j, k, iiii(j,k), lot(j,k), dload(j,k)

 enddo

 read(*,*) 



!
!
!>>>> Time history of the "rigid" inertia [unnormalized, for the moment]
!     this inertia variation is only "due to loading on a RIGID Earth" 
!     diner_r(k) are coefficients of a develpment in series of Heaviside functions
!
 do 2 k=0, nn

     diner_r(k) = (4d0*pi/3d0)*sqrt(6d0/5d0)*conjg(dload(j_index(2,1),k))

     write(*,*) 'rigid inertia ', k, diner_r(k)   !, dload(j_index(2,1),k) 

!write(*,*) 'load ', k, lot(j_index(2,1),k), dload(j_index(2,1),k)

2 continue
!
   read(*,*) 
!
!
!>>>> Time history of inertia [unnormalized, for the moment]
!     this inertia variation is only "due to loading on a Deformable Earth" 
!     diner_d(i) are inertia values at time "t_i" 
!
 do 3 i=0, nn
 do 3 k=0, i 
 	cnv(i,k) = 1d0 + klf 
 	do 3 jp=1, m 
        cnv(i,k) = cnv(i,k) + r(jp)*EXP(s(jp)*float(i-k)*delta)
3 continue
!
 do 4 i=0, nn

       diner_d(i) = (0d0,0d0)

       do 4 k=0, i 

!write(*,*) i, k, delta, cnv(i,k)

       diner_d(i) = diner_d(i) + diner_r(k)*cnv(i,k)

4 continue 
!
  do i=0, nn
     write(*,*) ' inertia ', i, diner_d(i)
  enddo
!
!>>>> Polar motion at time "t_i"
!
     do i=0, nn 
     
     PM (i)=(0d0,0d0)
     DPM(i)=(0d0,0d0)     
     
     do k=0, i 
     
     	PM(i)  = PM(i)  + diner_r(k) *  MFUNC(DELTA*DFLOAT(i-k))
     	DPM(i) = DPM(i) + diner_r(k) * DMFUNC(DELTA*DFLOAT(i-k))
 
     enddo
     
     PM (i) =  PM (i) * (ERAD**4/(CCCC-AAAA))
     DPM(i) =  DPM(i) * (ERAD**4/(CCCC-AAAA))
     
     enddo
!

 
     do i=0, nn  
   
        m1(i) = real  (pm(i))
        m2(i) = aimag (pm(i))

        m1dot(i) = real (dpm(i))
        m2dot(i) = aimag(dpm(i))
  
!
! Angle between the z-axis and the rotational vector (deg)  
     polar_disp = sqrt(m1(i)**2+m2(i)**2)*RAD2DEG
!
! Rate of change of the angle above (deg/Ma)     
     polar_rate = ((m1(i)*m1dot(i) + m2(i)*m2dot(i))/sqrt(m1(i)**2+m2(i)**2))*RAD2DEG*1000.

     write(55,'(i4,10(1x,E14.5))') i, m1(i),    m2(i),     atan2(m2(i),   m1(i)   )*RAD2DEG, & 
                                      m1dot(i), m2dot(i),  atan2(m2dot(i),m1dot(i))*RAD2DEG, & 
				      polar_disp, polar_rate 
         
     enddo
!
 END PROGRAM PMD
!
!
!
!
!
 FUNCTION MFUNC (X)
 USE DATA
 IMPLICIT NONE 
 INTEGER J, JP 
 COMPLEX*16 MA, MS, MFUNC
 REAL*8 X
!
! The function "MFUNC" - see my notes 
!
 MA = (0d0,0d0)
!
 MA = MA + AE*(1D0+KLE) + AS*(1D0+KLF)*X  
!
! DO J=1, MP 
!      MA = MA - (1D0+KLF)*(AA(J)/A(J))*(1D0-EXP(A(J)*X)) 
! ENDDO 
!
! DO JP=1, M
! DO J=1, MP
!      MA = MA + C(j,jp)*EXP(A(J)*X)
! enddo
! enddo
!
!
 MS = (0d0,0d0) 
!	 
! DO JP=1, M 
!     MS = MS - (AS*R(JP)/S(JP) + AE*R(JP)) * (1D0-EXP(S(JP)*X)) 
! ENDDO	 
! !
! DO JP=1, M
! DO J=1, MP
!       MS = MS - C(j,jp)*EXP(S(JP)*X)
! enddo
! enddo
!
!
 MFUNC = MS+MA 
!
!
 END FUNCTION MFUNC 
!
!
!
!
!
!
 FUNCTION DMFUNC (X)
 USE DATA
 IMPLICIT NONE 
 INTEGER J, JP 
 COMPLEX*16 DMA, DMS, DMFUNC
 REAL*8 X
!
!
! Time-derivative of function "MFUNC" - see my notes <<< DETAILS TO BE VERIFIED
!
!
 DMA = (0d0,0d0)
!
 DMA = DMA + AS*(1D0+KLF) 
!
! DO J=1, MP 
!     DMA = DMA + (1D0+KLF)*AA(J)*EXP(A(J)*X) 
! ENDDO 
!
! DO J=1, MP
! DO JP=1, M
!     DMA = DMA + C(J,JP)*A(J) *EXP( A(J)*X)
! enddo
! enddo
!
!
 DMS = (0d0,0d0)
!	 
! DO JP=1, M 
!      DMS = DMS + (AS*R(JP) + AE*R(JP)*S(JP)) * EXP(S(JP)*X)  
! ENDDO	 
! 
! DO JP=1, M
! DO J=1, MP
!       DMS = DMS - C(J,JP)*S(JP)*EXP(S(JP)*X)
! ENDDO
! ENDDO
!
! Write(*,*) "DMS  ", DMS
!
 DMFUNC = DMA + DMS 
!
!
 END FUNCTION DMFUNC 



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
! and computing some useful combination of Love numbers. GS AUG 2010 
!
  CHARACTER*22, PARAMETER :: FILEDEG2='degree_2_stuff.dat'  
  CHARACTER*22, PARAMETER :: FILEPMTF='PM_spectrum.dat'  
  CHARACTER*50 JUNKA
  REAL*8 REA, IMA, REAA, IMAA, KL0
  integer :: i, j, jp, k
  
  complex*16 app(m)
!
!
  Write(*,*) '    - Reading data from external files'
  Write(*,*) '    - A backup of data is reported on file <<rotational_data.bck>>'
!
  Open(13,file='rotational_data.bck',status='unknown') 
!   
  if(TLOVE==1) then 
!
  Write(*,*) '    - Reading degree 2 stuff from file ', trim(adjustl(FILEDEG2)) 
!                   [this file is created by the high-precision TABOO code]    
!
 	open (1,file=FILEDEG2,status='unknown') 

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
	read (1,*) s(i) 
        write(13,'(A23,1x,D20.8)') '    - s isostatic roots', s(i)		                         		
	enddo
!	
        close(1) 
 else 
	Write(*, *) "The polar motion in response to loading can  only be computed"
	Write(*, *) "if the <<degree 2 stuff>> file exists. The program will STOP "
	Write(88,*) "The polar motion in response to loading can  only be computed"
	Write(88,*) "if the <<degree 2 stuff>> file exists. The program will STOP "
        CALL STOP_CONFIG
 endif
!
!
!
!
!
 Write(*,*) '    - Reading PMTF data from ', trim(adjustl(FILEPMTF)) 
!
 if(TLOVE==1) then 
!
 	open (1,file=FILEPMTF,status='unknown') 

	call read_junk(7, 1)
!
	call read_junk(2, 1)
	READ(1,'(A4,1X,2(1X,D20.8))') JUNKA, AE 
        write(13,'(A23,4(1x,D20.8))') '    - Elastic R residue', AE
!
	call read_junk(2, 1)
	READ(1,'(A4,1X,2(1X,D20.8))') JUNKA, AS
        write(13,'(A23,4(1x,D20.8))') '    - Secular S residue', AS
!
	call read_junk(2, 1)
  	DO K=1, MP
		READ(1,'(I4,1X,8(1X,E20.8))')  I, REA, IMA, REAA, IMAA 
!		
		a (k) = dcmplx(rea,  ima)
		aa(k) = dcmplx(reaa, imaa)
!
        write(13,'(A23,4(1x,D20.8))') '    - R roots, residues', a(k), aa(k)
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
!
Write(*,*) '    - Computing useful constants '
!
          kl0=kle
          do jp=1, m 
               r(jp)=kl(jp)/s(jp)
	       kl0 = kl0 - r(jp)
!write(*,*) 'r ', jp, kl(jp), s(jp), r(jp)
          enddo
!write(*,*) 'k fluid ', kl0 
!read(*,*) 
!	  
          do j =  1, mp 
          do jp = 1, m 
               c(j,jp) = aa(j)*r(jp)/(a(j)-s(jp))    
!write(*,*) j, jp, c(j,jp)
          enddo
          enddo
!
!do jp=1, m 
!p(jp)=0.	
!do j=1, mp 
!p(jp) = p(jp) + c(j,jp)		  
!enddo 	
!write(*,*) 'p ', jp, p(jp)
!enddo
!
!do j=1, mp 
!q(j)=0.	
!do jp=1, m 
!q(j) = q(j) + c(j,jp)		  
!enddo 	
!write(*,*) 'q ', j, q(j)
!enddo


	 do jp=1, m 
	 
	 app(jp)= ae*r(jp) + as*r(jp)/s(jp)
	 
	 do j=1, mp 
	 
	 	app(jp) = app(jp) - c(j,jp)
	 
	 enddo

!write(*,*) 'AAAAAA ', jp, app(jp) 
	 
	 enddo

!
!
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
