!
! This is program "GEO_MAPS.F90" 
!
! Created by GS October 2008 for v. 2.7 [Intel]
! Reviewed GS & FC August 2009 -  "Varying coastlines" (reprise!) 
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
!     Re-touched by GS on April 2010 for te g95 implementation 
! === Revised GS & FC May 21 2011 - DELTA parameter in ice history     
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
! SELEN is distributed in the hope that it will be useful, but WITHOUT ANY 
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
! details. 
! 
! You should have received a copy of the GNU General Public License along 
! with SELEN.  If not, see <http://www.gnu.org/licenses/>.
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! 
! Produces "maps" of 3D velocity across regions specified in file "data.inc". 
! Units are mm/yr
!
! Input files: 
!	- 'shu.bin' 
!       - 'shv.bin'
!
! Output files:
!       - 'up-*' files for the up component of velocity 
!       - 'no-*' files for the north component 
!       - 'ea-*' files for the east  component 
!
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!INCLUDE "harmonics.f90"
 PROGRAM GEO_MAPS
 IMPLICIT NONE 
 INCLUDE "data.inc"

 CHARACTER*11, PARAMETER :: UN='UNKNOWN', UF='UNFORMATTED' 
 CHARACTER*4,  PARAMETER :: SUFFIX = ".pix"
 CHARACTER*13, PARAMETER :: FMT='(3(f10.6,1X))'
 REAL*4, PARAMETER :: MIN_MAG = 0.10
 CHARACTER*30 HEADER 
 INTEGER, PARAMETER :: NHEADER=5, RES_MAX=48, & 
                       NP_MAX=2*RES_MAX*(RES_MAX-1)*20+12, NMAX_3D_REGIONS=10
 INTEGER I, J, K, DOM, MJ, IREG, IMEM, IUNIT, DM(JMAX)  
 INTEGER NPX_EFF(NMAX_3D_REGIONS)
 REAL*4 LONP(NP_MAX), LATP(NP_MAX), RATE_UP(NP_MAX), RATE_NO(NP_MAX), RATE_EA(NP_MAX) 
 REAL*4 LONX, LATX, N_VEL, E_VEL, AZIMUTH, VEL_MOD  
 COMPLEX*8 CU(JMAX,0:NN), CV(JMAX,0:NN)
 COMPLEX*8 Y(JMAX,NP_MAX) 
! 
 REAL*8, PARAMETER :: DDELTA=DELTA 
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! --- Pre-computing the delta (0,m)
 do j=1, jmax 
	dm(j)=2-dom(j) 
 enddo	
!
!  
!
!Write(*,*) '    - geo_maps.f: Reading the harmonics from binary files'
 open (101,file='shu.bin',status='UNKNOWN',form='UNFORMATTED') 
 open (102,file='shv.bin',status='UNKNOWN',form='UNFORMATTED') 
 read (101) CU ; read (102) CV 
 close(101)    ; close(102)   
!
!
 iunit=100 
!
! ### Loop on the regions ### 
!
 Do 100 ireg = 1, n_3d_regions 
!
 Write(*,*) '    - Region: ', trim(adjustl(tred_regions_name(ireg)))
!
 	open(91,file='up-'//trim(adjustl(tred_regions_name(ireg)))//'.dat',status='UNKNOWN') 
 	open(92,file='no-'//trim(adjustl(tred_regions_name(ireg)))//'.dat',status='UNKNOWN') 
 	open(93,file='ea-'//trim(adjustl(tred_regions_name(ireg)))//'.dat',status='UNKNOWN') 
!
        iunit=iunit+1 
!
	open(iunit,file=trim(adjustl(tred_regions_name(ireg)))//SUFFIX,status='old') 
!
!
! --- Header lines 
	do 10 i=1, nheader 
		read(iunit,'(a30)')HEADER 
 10     continue   
!
! --- Pixels coordinates  
	npx_eff(ireg)=0
        do 20 i=1, np_max 	 
 		read(iunit,*,end=50) lonp(i), latp(i)  
		npx_eff(ireg)=npx_eff(ireg)+1
 20     continue 
!
 50 close(iunit)  
!
	do 800 imem=1, 3     
!
           do i=1, npx_eff(ireg)    
!
               if(imem==1) call             harmo (lmax, lonp(i), latp(i), y(:,i)) 
               if(imem==2) call  grad_theta_harmo (lmax, lonp(i), latp(i), y(:,i))
               if(imem==3) call grad_lambda_harmo (lmax, lonp(i), latp(i), y(:,i))
           enddo
!
! --- LOOP on the pixels for each region 
	do 30 i=1, npx_eff(ireg)
!
       	   if(imem==1)  rate_up(i)=0.  ! 'up' 
       	   if(imem==2)  rate_no(i)=0.  ! 'north' 
       	   if(imem==3)  rate_ea(i)=0.  ! 'east'  
!
        do 40 j=1, jmax 
           if(imem==1) rate_up(i) = rate_up(i) + dm(j)*real(((cu(j,nn)-cu(j,nn-1))/DDELTA)*y(j,i))  
           if(imem==2) rate_no(i) = rate_no(i) + dm(j)*real(((cv(j,nn)-cv(j,nn-1))/DDELTA)*y(j,i))  
           if(imem==3) rate_ea(i) = rate_ea(i) + dm(j)*real(((cv(j,nn)-cv(j,nn-1))/DDELTA)*y(j,i)) 
 40     continue 
!
		rate_no(i) = - rate_no(i) 
!
	if(imem==1) write(91,*) lonp(i), latp(i), rate_up(i) 
	if(imem==2) write(92,*) lonp(i), latp(i), rate_no(i) 
	if(imem==3) write(93,*) lonp(i), latp(i), rate_ea(i) 	
!
 30     continue 
!
	if(imem==1) close(91)
	if(imem==2) close(92)
	if(imem==3) close(93) 
! 
 800 continue 
!
 100 continue  
!
!
!
! ### For plotting purposes (see -SV option of psxy), we need a 4-columns file with:
!     lon, lat, azimuth (wrt north), and horizontal vector magitude. This file is created 
!     below. Notice that vectors with magnitude smaller than MIN_MAG will not be considered 
!     for plotting. This is done for cleaning a bit the resulting maps (according to my 
!     taste).    === GS  October 19, 2008  ===
!
 Do 200 ireg = 1, n_3d_regions 
!
! --- Input 
 	open(92,file='no-'//trim(adjustl(tred_regions_name(ireg)))//'.dat',status='old') 
 	open(93,file='ea-'//trim(adjustl(tred_regions_name(ireg)))//'.dat',status='old') 
!
! --- Output 
 	open(94,file='az-'//trim(adjustl(tred_regions_name(ireg)))//'.dat',status='UNKNOWN') 
!
!
	do 60 i=1, npx_eff(ireg)
!        
		read(92,*) lonx, latx, N_vel 			
		read(93,*) lonx, latx, E_vel 	

	 	vel_mod=sqrt(N_vel**2+E_vel**2)
		
		if(vel_mod >= min_mag .and. N_vel /= 0.) then
!		 
				azimuth = atan2(E_vel, N_vel)*180./pi		
!
				Write(94,*) lonx, latx, azimuth, vel_mod 
		endif   
!
!
60  continue 
!
 close(92)
 close(93)
 close(94) 
!
200 continue
!
!
!
 END PROGRAM GEO_MAPS
!
!
!
!
!
!
