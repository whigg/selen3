!
! This is program "VELSN.F90" 
!
! Created by GS September 2008 for version 2.7 
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
! *** Revised GS July 2010 - g95 - Double precision implementation
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
! SELEN is distributed in the hope that it will be useful, but WITHOUT ANY 
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
! details. 
! 
! You should have received a copy of the GNU General Public License along 
! with SELEN.  If not, see <http://www.gnu.org/licenses/>.
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!  This program first computes the 4pi-normalized, complex spherical harmonics 
!  at geodetic sites listed in file GEODETIC_DATABASE. Subsequantly, it computes 
!  the 3D displacement (UP, NORTH, EAST), and the scalars S and N at the same 
!  sites. This set of five quantites specifies completely the "geodetic state" 
!  of the site. 
!
!  Input files:
!	- A list of geodetic sites (e. g., "gps-sites.dat")
!  
!  Output files: 
!	- ...
!       - ...
!       - ...
!
!
 INCLUDE "harmonics.f90"
 PROGRAM GEO
 IMPLICIT NONE 
 INCLUDE "data.inc"
 CHARACTER*22 NAME, CJUNK 
 INTEGER I, J, NDATA, CODE(NGEOD) 
 REAL*8 AJUNK, LONS(NGEOD), LATS(NGEOD)
 COMPLEX*16 CS(JMAX,0:NN), CU(JMAX,0:NN), CN(JMAX,0:NN), CV(JMAX,0:NN)
 COMPLEX*16 Y(JMAX,NGEOD)         
!
!
!
  write(*,*) '    - velsn.f: Opening and reading the geodetic database ...'
  OPEN(1,FILE=GEODETIC_DATABASE,STATUS='unknown')
  Do i=1, ngeod
       read(1,'(a22)')cjunk 
       read(1,'(a22)')cjunk 
       read(1,*) lats(i), lons(i)  		
       read(1,'(a22)')cjunk 	
  Enddo
  Close(1)
!
!
  write(*,*) '    - velsn.f: Computing the harmonics at the geodetic sites ...'
  do i=1, NGEOD    
       If(mod(i,15)==0) write(*,*) '    - sh_rsl.f:', i, 'RSL sites of', nrsl
       call harmo(lmax, lons(i), lats(i), y(:,i)) 
  enddo
!
!
!---- Computing S-dot, U-dot, and N-dot at pixels 
!  
 Write(*,*) "    - velsn.f: Reading the harmonics from binary files"
 open (101,file='shs.bin',status='UNKNOWN',form='UNFORMATTED') 
 open (102,file='shu.bin',status='UNKNOWN',form='UNFORMATTED') 
 open (103,file='shn.bin',status='UNKNOWN',form='UNFORMATTED') 
 open (104,file='shv.bin',status='UNKNOWN',form='UNFORMATTED') 
 read (101) CS
 read (102) CU
 read (103) CN
 read (103) CV
 close(101) ; close(102) ; close(103) ; close(104) 
!
!
  write(*,*) "    - sh_rsl.f: Writing the harmonics on file 'shrsl.bin'"
  open(7,file='shgeo.bin',status='unknown',form='unformatted') 
  write(7) y; close(7) 
!
!
!
 END PROGRAM GEO
!
!
!
!
!
!
