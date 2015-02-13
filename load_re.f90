!
! This is program "LOAD_RE.F90" 
!
!               Created by GS on 14 August 2010 (SELEN 3.2) 
! March 2012: Implementation of the numerical derivative "in the future"
!               ... ... ...  
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
!
! Computes the load for a RIGID EARTH - This 
!
! Input files:
!       - ... 
!	- ...  
!
! Output files:
!       - ... 
!	- ... 
!
!INCLUDE "harmonics.f90"
 PROGRAM LOAD_RE 
 IMPLICIT NONE 
 INCLUDE "data.inc"			    
!
! ====== Modified March 2012 GS: N+1
!
 COMPLEX*16 IIII(JMAX,0:NN+1), ZE(JMAX,0:NN+1)
 COMPLEX*16 OC(JMAX) 
 REAL*4 RHOI_O_RHOW
 REAL*8 RESH, IMSH
 CHARACTER*22 Out_Filename  
 INTEGER I, J, K 
! 
! - Ice component of load for a rigid Earth ------------------------ Updated 3/12
 complex*16 LO_REI(JMAX,0:NN+1)
 CHARACTER*22, PARAMETER :: FILELI='shload_rice.bin'
!
! - Ocean component of load for a rigid Earth  --------------------- Updated 3/12 
 complex*16 LO_REO(JMAX,0:NN+1)
 CHARACTER*22, PARAMETER :: FILELO='shload_roce.bin' 
!
! - Total load for a rigid Earth ----------------------------------- Updated 3/12
 complex*16 LO_RET(JMAX,0:NN+1)
 CHARACTER*22, PARAMETER :: FILELT='shload_re.bin'  
!
!
! //////////////////////////////////////////////////////////////////
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
! //////////////////////////////////////////////////////////////////
!
!
! # Ice to water density ratio
!
 RHOI_O_RHOW = RHOI/RHOW 
!
!
! # Reading the I array ...
!
 Write(*,*) "    - Reading array 'I'"
     open(1,file='shice.dat',status='unknown') 
     read(1,*) IIII
 close(1)
!
!
! # Reading the SH OF coefficients from shof.dat 
!
 Write(*,*) '    - Reading the SH OF coeff. from shof.dat'
 open(3,file='shof.dat',status='unknown')
 do j=1, jmax   
     read(3,*) k, resh, imsh 
     oc(j)=cmplx(resh, imsh)	
 enddo
 close(3)
!
!
! # Computing the eustatic Z array...-- Updated 3/12
!
 ze(:,:) = (0d0,0d0) 		
 do k=0,nn+1	
     ze(:,k) = - RHOI_O_RHOW*(iiii(1,k)/oc(1))*oc(:)
 enddo
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Load components DRAFT  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Load components DRAFT  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Load components DRAFT  
!
 LO_REI(:,:)  = (0.d0,0.d0)
 LO_REO(:,:)  = (0.d0,0.d0)
 LO_RET(:,:)  = (0.d0,0.d0)
!
! --- Ice  component of load    ======== RIGID EARTH 
 LO_REI(:,:) = RHOI*IIII(:,:)
!
!
! --- Ocean component of load   ======== RIGID EARTH
 LO_REO(:,:) = RHOW*ZE(:,:)
!
!
! --- Total load (ice + ocean)  ======== RIGID EARTH
 LO_RET(:,:) = LO_REI(:,:) + LO_REO(:,:)
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Load components DRAFT  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Load components DRAFT  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Load components DRAFT  
!
!
! --- Ice component of load                                              <<<< new 
 Out_Filename=fileli 
	open(3,file=Out_Filename,status='unknown',form='unformatted') 
        write(3) LO_REI 
 close(3)
!
! --- Ocean component of load          				         <<<< new 
 Out_Filename=filelo 
 	open(3,file=Out_Filename,status='unknown',form='unformatted') 
        write(3) LO_REO 
 close(3)
!
! --- Total (ice+ocean) load             			         <<<< new 
 Out_Filename=filelt 
 	open(3,file=Out_Filename,status='unknown',form='unformatted') 
        write(3) LO_RET 
 close(3)
!
 END PROGRAM LOAD_RE 
!
!
!
!
!
!
