!
! This is program "PMD_MOD.F90" 
!
! GS July 29 2010 for version 3.2 of SELEN 
! *** Revised GS August 2010 - g95 - Double precision implementation
! *** Revised GS August 2010 - A number of issues (see revised theory)
! Revised GS Sept 3 2010 for the implementation of "m_3"
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
  REAL*8  KLE                      ! Elastic degree 2 k Love number 
  REAL*8  KLF                       ! Fluid degree 2 k Love number  
  REAL*8  KL(M)                      ! Residues of degree 2 k Love number    
  REAL*8  S(M)                        ! Degree 2 isostatic poles 
  REAL*8  AEP                        ! (Primed) elastic residue of the PMTF 
  REAL*8  ASP                       ! (Primed) secular residue of the PMTF 
  COMPLEX*16 A(MP)                 ! Rotational roots   
  COMPLEX*16 AAP(MP)              ! (Primed) rotational residues 
!
  REAL*8, PARAMETER :: CCCC=8.0394d37, AAAA=8.0131d37, OMEGA=7.292115d-5
  REAL*8, PARAMETER :: SECY=365.25d0*24d0*3600d0*1000d0
  REAL*8, PARAMETER :: SR=((CCCC-AAAA)/AAAA)*OMEGA*SECY  
  REAL*8, PARAMETER :: LOD=86400.
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
 PROGRAM PMD 
!
 USE DATA
 IMPLICIT NONE 
!
 CHARACTER*22, PARAMETER :: FILELOAD='shload_re.bin'   
 REAL*8, PARAMETER :: ERAD=6.371D6
 REAL*8, PARAMETER :: RAD2DEG=180D0/PI
!
 INTEGER I, J, K, JP, J_INDEX
 COMPLEX*16 LOT (JMAX,0:NN), DLOAD(JMAX,0:NN), DINER_R(0:NN)
 COMPLEX*16 MFUNC, PM(0:NN), DMFUNC, DPM(0:NN)
 REAL*8 M3FUNC, DINER33_R(0:NN)
!
 REAL*8 M1(0:NN),    M2(0:NN),    M3(0:NN),   & 
        M1DOT(0:NN), M2DOT(0:NN), DLOD(0:NN), & 
	POLAR_DISP,  POLAR_RATE,  TIME_BP
!
 CHARACTER*20 DATE, TIMC
! 
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!
 Write(*,*) '    - Retrieving information from existing SELEN data:' 
!
 Write(*,*) '    - Reading the harmonics of *total* surface load from file ', & 
                   trim(adjustl(FILELOAD)) 
!
! File with load information 
 open (1,file=FILELOAD,status='unknown',form='unformatted') 
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
 do 1 k=0, nn 	  
     if(k==0) then 
         dload(j,k) = lot(j,k)
	      else 
         dload(j,k) = lot(j,k) - lot(j,k-1)
     endif	 
1 continue
!
!
!
!
!
!
!>>>> Time history of the "rigid" inertia, only "due to loading on a RIGID Earth" 
!     diner_r(k) are coefficients of a develpment in series of Heaviside functions
!
  do 4 k=0, nn
!
! This rigid inertia is a complex variable, related with 
! polar motion. Only the degree 2 and order 1 components 
! of the load variation are involved -See my theory notes.
!
 diner_r(k) =     (4d0*pi/3d0)*sqrt(6d0/5d0)*conjg(dload(j_index(2,1),k))
!
! This rigid inertia is a *real* variable, related with 
! lenght of day (LOD). Only the degree 0 and degree 2 order 
! 0 components of load are involved  -See my theory notes.

 diner33_r(k) =   (8d0*pi/3d0)   *real(dload(j_index(0,0),k)) - & 
                  (1d0/sqrt(5d0))*real(dload(j_index(2,0),k))
!
4 continue
!
!
!//////////////////////////////////////////////////////////!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
!							   !
!   Polar motion and rate of polar motion at time "t_i"    !
!							   ! 
!//////////////////////////////////////////////////////////!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
!
 call DATE_AND_TIME (date,timc)
!
 OPEN (55,file='m.dat',status='unknown')   
!
 Write(55,*) "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " 
 Write(55,*) '# File <<m.dat>>, created by program PMD_MOD.F90 on ', & 
              date(1:4), '.', date(5:6), '.', date(7:8), & 
	      '      time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
 Write(55,*) "# Note: the Chandler Wobble is NOT included"
 Write(55,*) "# " 
 Write(55,*) "#     time      time BP        m1             m2           Arg (m)        Mod (m)"	
 Write(55,*) "#      ka         ka          deg            deg             deg            deg  "	  
 Write(55,*) "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " 
!
 OPEN (57,file='m.dot',status='unknown')   
!
 Write(57,*) "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " 
 Write(57,*) '# File <<m.dot>>, created by program PMD_MOD.F90 on ', & 
              date(1:4), '.', date(5:6), '.', date(7:8), & 
              '      time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
 Write(57,*) "# Note: the Chandler Wobble is NOT included"
 Write(57,*) "# " 
 Write(57,*) "#     time      time BP       m1_dot         m2_dot      Arg (m_dot)    Mod (m_dot)"   
 Write(57,*) "#      ka         ka          deg/Ma         deg/Ma          deg          deg/Ma   "      
 Write(57,*) "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " 
!
     do 5 i=0, nn 
!     
     PM (i)=(0.d0,0.d0)
     DPM(i)=(0.d0,0.d0)     
!     
        do 6 k=0, i 
!     
     	PM(i)  = PM(i)  + diner_r(k) *  MFUNC(DELTA*DFLOAT(i-k))
     	DPM(i) = DPM(i) + diner_r(k) * DMFUNC(DELTA*DFLOAT(i-k))
! 
 6      continue
!     
        PM (i) =  PM (i) * (ERAD**4/(CCCC-AAAA))
        DPM(i) =  DPM(i) * (ERAD**4/(CCCC-AAAA))
!     
 5   continue 
!
     do 7 i=0, nn  
!   
        m1(i) = real  (pm(i))
        m2(i) = aimag (pm(i))
!
        m1dot(i) = real (dpm(i))
        m2dot(i) = aimag(dpm(i)) 
!
! Angle between the z-axis and the rotational vector (deg)  
!
        polar_disp = sqrt(m1(i)**2+m2(i)**2)*RAD2DEG
!
!
! Rate of change of the angle above (deg/Ma)   
!  
        polar_rate = ((m1(i)*m1dot(i) + m2(i)*m2dot(i))/sqrt(m1(i)**2+m2(i)**2))*RAD2DEG*1000.
!
        time_bp=float(nn)-float(i)*delta
!
        write(55,'(2(4X, F8.4), 10(1x,E14.5))') float(i), & 
					     time_bp,  & 
	        			     m1(i)*RAD2DEG, & 
					     m2(i)*RAD2DEG, & 
					     atan2(m2(i),   m1(i)   )*RAD2DEG,  &
					     POLAR_DISP 
!        
	write(57,'(2(4X, F8.4), 10(1x,E14.5))') float(i), & 
					     time_bp,  & 	
				             m1dot(i)*RAD2DEG*1000., & 
					     m2dot(i)*RAD2DEG*1000., & 
					     atan2(m2dot(i),m1dot(i))*RAD2DEG,  & 
					     POLAR_RATE  
!        
 7    continue 
!
 Write(55,*) "# Contact: giorgio.spada@gmail.com"
 Write(57,*) "# Contact: giorgio.spada@gmail.com"
!
 Close (55) 
 Close (57) 
!
!
!////////////////////////////////////////////!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
!					     !
!   "m_3" and LOD variation at time "t_i"    !
!                                            ! 
!////////////////////////////////////////////!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
!
 call DATE_AND_TIME (date,timc)
!
 OPEN (55,file='lod.dat',status='unknown')   
!
 Write(55,*) "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " 
 Write(55,*) '# File <<lod.dat>>, created by program PMD_MOD.F90 on ', & 
              date(1:4), '.', date(5:6), '.', date(7:8), & 
	      '      time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
 Write(55,*) "# "
 Write(55,*) "# " 
 Write(55,*) "#     time      time BP        m3          D(LOD) "	
 Write(55,*) "#      ka         ka          rad            ms   "	  
 Write(55,*) "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " 
!
!
 OPEN (57,file='lod.dot',status='unknown')   
!
! Write(55,*) "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " 
! Write(55,*) '# File <<m.dat>>, created by program PMD_MOD.F90 on ', & 
!              date(1:4), '.', date(5:6), '.', date(7:8), & 
!	      '      time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
! Write(55,*) "# Note: the Chandler Wobble is NOT included"
! Write(55,*) "# " 
! Write(55,*) "#     time      time BP        m1             m2           Arg (m)        Mod (m)"	
! Write(55,*) "#      ka         ka          deg            deg             deg            deg  "	  
! Write(55,*) "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " 
!
!
     do 15 i=0, nn 
!     
     M3 (i)= 0.d0
!     
        do 16 k=0, i 
!     
     	M3(i)  = M3(i)  + diner33_r(k) *  M3FUNC(DELTA*DFLOAT(i-k))
! 
 16      continue
!  
! The "-" sign here is only something more than a guess - Check the theory again... 
!   
        M3  (i)  =  - M3 (i) *ERAD**4/CCCC
	DLOD(i)  =  - M3 (i) *LOD   
!     
 15   continue 
!
!
     do 17 i=0, nn  
!   
!
!
! Rate of change of the angle above (deg/Ma)   
!  
!        polar_rate = ((m1(i)*m1dot(i) + m2(i)*m2dot(i))/sqrt(m1(i)**2+m2(i)**2))*RAD2DEG*1000.
!
        time_bp=float(nn)-float(i)*delta
!
        write(55,'(2(4X, F8.4), 10(1x,E14.5))') float(i), & 
					        time_bp,  & 
	        			        m3(i),    & 
						DLOD(i)*1000.D0 
!        
 17    continue 
!
!
!
! Work in progress! Work in progress! Work in progress! Work in progress! 
! Work in progress! Work in progress! Work in progress! Work in progress! 
! Work in progress! Work in progress! Work in progress! Work in progress! 
! Work in progress! Work in progress! Work in progress! Work in progress! 
!
!
 END PROGRAM PMD
!
!
!
!
!
!
 FUNCTION M3FUNC(X)
 USE DATA
 IMPLICIT NONE 
 INTEGER J, JP 
 REAL*8 M3FUNC, X
!
! ---------------------------------------------------------------
! The function "M3FUNC" - see my notes ... (those revised) 
! GS September 2010 - 
! ---------------------------------------------------------------
!
  M3FUNC = 1D0 + KLF 
!
  DO J=1, M 
      M3FUNC = M3FUNC + (KL(J)/S(J))*(EXP(S(J)*X)-1D0) 
  ENDDO  
!
 END FUNCTION M3FUNC 
!
!
!
!
!
!
 FUNCTION MFUNC(X)
 USE DATA
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
  MFUNC = MFUNC + ASP*X 
!
  DO J=1, MP 
      MFUNC = MFUNC + (AAP(J)/A(J))*(EXP(A(J)*X)-1D0) 
  ENDDO  
!
 END FUNCTION MFUNC 
!
!
!
!
!
!
 FUNCTION DMFUNC(X)
 USE DATA
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
  DMFUNC = ASP
!
  DO J=1, MP 
      DMFUNC = DMFUNC + AAP(J)*EXP(A(J)*X) 
  ENDDO  
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
  Write(*,*) '    - A backup of data is reported on file <<Rotational_data.bck>>'
!
  Open(13,file='Rotational_data.bck',status='unknown') 
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
		a  (k) = dcmplx(rea,  ima)
		aap(k) = dcmplx(reaap, imaap)
!
        write(13,'(A23,4(1x,D20.8))') '    - R roots & *primed* residues', a(k), aap(k)
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
