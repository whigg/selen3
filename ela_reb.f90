!
! This is program "ELA_REB.F90" 
!
! Created GS July 2010 
! Elastic rebound in Greenland and in the surroundings... 
! *** Revised GS July 2010 - g95 - Double precision implementation
! *** Revised GS December 2010 - OUTPUT is NOW in MILLIMITERS (MM), not in cm. 
! === Revised GS & GG Apr 2011 - for the ANT component
! === Revised GG Apr 25 2011 - for the ANT component  
! June 30, 2012; regular grid for the one-step elastic rebound (EOS paper?)
! Modified by DM on Nov 09 2012 - dynamic memory allocation & external NANCH
! === Also modified for REAL*8 FINDPX()
! Modified DM April 05 2013 - OpenMP parallelization
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
! Produces maps of the rate of present-day sealevel change (dot S), of 
! vertical uplift (dot U), and of geoid change (dot N). Units are mm/yr
!
! Input files: 
! [For global maps, these are the harmonics on the SLE pixels]
! 	- shs.bin
!	- shu.bin
!	- shn.bin
!	- shz.bin 
!	- shf.bin 
! 	- sh.bin
! 	- px-table.dat
!
! Output files:
! 	- sdotmap.dat [to be updated...]
!	- udotmap.dat
!	- ndotmap.dat
!	- fadotmap.dat
!	- ssdotmap.dat
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!INCLUDE "harmonics.f90"
 PROGRAM ELA_REB
 IMPLICIT NONE 
 INCLUDE "data.inc"
!
 CHARACTER*13, PARAMETER :: FMT='(3(f10.6,1X))'
 CHARACTER*11, PARAMETER :: UF='UNFORMATTED', F='FORMATTED', UN='UNKNOWN' 
 REAL*8,  PARAMETER :: EPS = 0.01
 INTEGER  :: NANCH
 INTEGER, PARAMETER :: Large_Integer=1E7
 REAL*8, PARAMETER :: from_m_to_mm=1000.
!
 INTEGER I, J, JJ, K, IP, NPP, NIN, RES_MAP, DOM, MJ, NOUT
 INTEGER, ALLOCATABLE :: MM(:), ANC(:), DM(:)  
 REAL*8 LONX, LATX
 REAL*8, ALLOCATABLE :: LON(:), LAT(:), ALF(:,:)
 REAL*8, ALLOCATABLE :: RATE_S(:), RATE_U(:), RATE_N(:)
 REAL*8 SLC_AVE
 REAL*8 RESH, IMSH
 LOGICAL LOGI1, LOGI2 
 COMPLEX*16, ALLOCATABLE :: ARMOY(:), OC(:) 
 CHARACTER*30 HEADER
! REAL*4  LONXX(large_integer),  LATXX(large_integer)
! REAL*8 DLONXX(large_integer), DLATXX(large_integer)
 REAL*8, ALLOCATABLE ::  LONXX(:),  LATXX(:)
 INTEGER, ALLOCATABLE :: IDX(:)
! 
 COMPLEX*16, ALLOCATABLE :: CS(:,:),  CU(:,:), CN(:,:)
 COMPLEX*16, ALLOCATABLE :: LONG_TABLE(:,:)
!
!
! --- Data for the computations on a 1 deg x 1 deg grid - June 2012
 Integer kh 
 Integer, parameter :: delta_grid=1 
 Integer, parameter :: NGRID=180/delta_grid
 REAL*8 LONGC, LATC
 COMPLEX*16, ALLOCATABLE :: Y(:,:)
!
! --- Data for the regional analysis of Greenland 
 REAL*4, PARAMETER :: LONMIN_GR=285.0
 REAL*4, PARAMETER :: LONMAX_GR=350.0
 REAL*4, PARAMETER :: LATMIN_GR= 58.0
 REAL*4, PARAMETER :: LATMAX_GR= 86.0
! REAL*4, PARAMETER :: RES_REG_GRC = 25.   ! Now from DATA.INC
!
!
! --- Data for the regional analysis of Antarctica 
 REAL*4, PARAMETER :: LATMIN_AR= -90.
 REAL*4, PARAMETER :: LATMAX_AR=   0.
! REAL*4, PARAMETER :: RES_REG_ANC = 150.  ! Now from DATA.INC
 
!
!
	   
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!
! --- Getting info about the number of anchor pixels 
!
  open(1,file='anchor.tmp',status='old')
  read(1,*) nanch
  close(1)
!
! 
 ALLOCATE( MM(JMAX), DM(JMAX) )
 ALLOCATE( ALF(JMAX,NANCH) )
 ALLOCATE( ARMOY(JMAX), OC(JMAX) )
 ALLOCATE( CS(JMAX,0:NN),  CU(JMAX,0:NN), CN(JMAX,0:NN)  )
 ALLOCATE( LONG_TABLE(0:LMAX,NP)  )
 ALLOCATE( Y(JMAX,1) )
!
!
! >>>> Pre-computing the degree 'm' corresponding to 'J' 
 do j=1, jmax 
     mm(j)=mj(j)
     dm(j)=2-dom(j) 
 enddo	
!
!
! >>>> Reading the SH OF coefficients from shof.dat 
!
 Write(*,*) "    - Reading the SH OF coeff. from <<shof.dat>>"
 open(3,file='shof.dat',status='unknown')
 do j=1, jmax   
	read(3,*) k, resh, imsh 
        oc(j)=cmplx(resh, imsh)	
 enddo
 close(3)
!
!
! >>>> Files with SH coefficients from SLE.F90
!
 open (1,file='shs.bin',status='unknown',form=uf) 
 read (1) CS
 close(1) 
 open (1,file='shu.bin',status='unknown',form=uf) 
 read (1) CU
 close(1) 
 open (1,file='shn.bin',status='unknown',form=uf) 
 read (1) CN
 close(1) 
!
! //////////////////////////////////////////////////////////////////  Global
! //////////////////////////////////////////////////////////////////  Global 
!
! >> For global analyses, I need to input the spherical harmonics at
!    pixels (note: the same set of pixels used for solving the SLE). 
!
! //////////////////////////////////////////////////////////////////  Global
! //////////////////////////////////////////////////////////////////  Global  
!
 IF(ELAREB_GG==1.or.ELAREB_AG==1.or.ELAREB_SM==1) then 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         IF(ELAREB_GRID==1) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!---- Reading the pixels table file... -----------------------------  Global
!
    nout = np
    allocate( lon(np), lat(np), anc(np) )
    allocate( rate_s(np), rate_u(np), rate_n(np) )
!    
 	open(2,file='px-table.dat',status='unknown')
 	do i=1, 4 
    		read(2,'(a30)') HEADER 	
 	enddo 
 	do i=1, np
    		read(2,*) lon(i), lat(i), anc(i)
 	enddo 
 	close(2)  
!
!---- Reading the table of spherical harmonics from <<sh.bin>> -----  Global
 	Write(*,*) "    - Reading ALFs & TRIGs from sh.bin"
 	open(3,file='sh.bin',status='unknown',form='unformatted')  
 		read(3)ALF
 		read(3)LONG_TABLE
 	Close(3) 
 	DO J=1, JMAX 
		ALF(J,:)=ALF(J,:)*DM(J)
 	ENDDO
!	
!!!!!!!!!!!!!!!!!!!!!!!!
	ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!	
!
!---- Average SLC 
!
	 slc_ave=0.
!
 	 do j=1, jmax 
 	       slc_ave = slc_ave + dm(j)*real(oc(j)*conjg(cs(j,nn)))/oc(1)  
 	 enddo
         write(*,*) "    - AVERAGE Sea Level Change (units are mm)= ", slc_ave*1000. 
         write(*,*) "    - [this is computed by a surface integral over the oceans]"
!
!
!
!
         IF(ELAREB_GRID==1) THEN      
!
		RATE_S   = 0.
		RATE_U   = 0.
		RATE_N   = 0.
!
        Write(*,*) "    - Computing S, U & N at ", np, " pixels "
!
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP PRIVATE(I,J) SHARED(RATE_S,RATE_U,RATE_N,ALF,ANC,LONG_TABLE,MM,CS,CU,CN) &
!$OMP SCHEDULE(GUIDED)
 	    DO I=1, NP 
	    DO J=1, JMAX 
  	  	  rate_s(i) = rate_s(i)   + ALF(j,anc(i))*real((cs(j,nn))*LONG_TABLE(MM(J),I))   
		  rate_u(i) = rate_u(i)   + ALF(j,anc(i))*real((cu(j,nn))*LONG_TABLE(MM(J),I)) 
		  rate_n(i) = rate_n(i)   + ALF(j,anc(i))*real((cn(j,nn))*LONG_TABLE(MM(J),I)) 
        END DO
        END DO
!$OMP END PARALLEL DO
!
      DEALLOCATE( ANC )
!
      ENDIF
!
!
! IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - 
! IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - 
!
        IF(ELAREB_GRID==2) THEN
!
! Define grid pixels
!
        nout = 2*ngrid*ngrid
        allocate( lon(2*ngrid*ngrid), lat(2*ngrid*ngrid) )
        allocate( rate_s(2*ngrid*ngrid), rate_u(2*ngrid*ngrid), rate_n(2*ngrid*ngrid) )
!
        k = 0
        do i=1, 2*NGRID 
	       longc = float(delta_grid)/2.*(2*i-1)
	       do j=1,NGRID
  		      latc = - 90. + float(delta_grid)/2.*(2*j-1)
              k = k + 1
              lon(k) = longc
              lat(k) = latc
           end do
        end do
!
! From the south pole (-90^) to the north pole (+90^)... 
! [modified July 15, 2012]
!
!
! --- Call HARMO once outside the parallel loop to initialize 
!     the SAVEd arrays
!
  call HARMO(LMAX,lon(1),lat(1),ARMOY)  
!	
!$OMP PARALLEL DO DEFAULT(PRIVATE) &
!$OMP SHARED(RATE_S,RATE_U,RATE_N,LON,LAT,CS,CU,CN) &
!$OMP SCHEDULE(GUIDED)
         do k=1,2*NGRID*NGRID
             if( .not. allocated(armoy) ) allocate(armoy(jmax))
        	 call harmo (lmax, lon(k), lat(k), armoy) 	
         end do
!$OMP END PARALLEL DO
!
         ENDIF
!
!
! IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - 
! IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - 
!
!
!
!
!
!
!---- Opening the files for output. For SLC, a normalized global map is also produced...  
!
	IF(ELAREB_GG==1) then 
		open (101,file= 'smap_gree_glob.dat',       status='unknown') 
		open (102,file= 'umap_gree_glob.dat',       status='unknown') 
		open (103,file= 'nmap_gree_glob.dat',       status='unknown') 
		open (104,file= 'smap_gree_glob_norm.dat',  status='unknown') 
	Endif
!
	IF(ELAREB_AG==1) then 
		open (101,file= 'smap_anta_glob.dat',       status='unknown') 
		open (102,file= 'umap_anta_glob.dat',       status='unknown') 
		open (103,file= 'nmap_anta_glob.dat',       status='unknown') 
		open (104,file= 'smap_anta_glob_norm.dat',  status='unknown') 
	Endif
!
	IF(ELAREB_SM==1) then 
		open (101,file= 'smap_glac_glob.dat',       status='unknown') 
		open (102,file= 'umap_glac_glob.dat',       status='unknown') 
		open (103,file= 'nmap_glac_glob.dat',       status='unknown') 
		open (104,file= 'smap_glac_glob_norm.dat',  status='unknown') 
	Endif
!
!
! 
    DO I=1,NOUT
!
!
! --- South pole  ------------------------  Global
            If (lat(i)==-90.) then  
       		latx=lat(i)+eps
       		do k=1, 4 
                	lonx=(k-1)*90.   		
     	        	write(101,FMT) lonx, latx, rate_s(i)*from_m_to_mm 
	        	write(102,FMT) lonx, latx, rate_u(i)*from_m_to_mm  
	        	write(103,FMT) lonx, latx, rate_n(i)*from_m_to_mm  
     	        	write(104,FMT) lonx, latx, rate_s(i)/slc_ave  						
       		enddo
!
! --- North pole ------------------------  Global
            elseif (lat(i)==+90.) then   
       		latx=lat(i)-eps
       		do k=1, 4 
			 lonx=(k-1)*90.   			
            	         write(101,FMT) lonx, latx, rate_s(i)*from_m_to_mm  
            	         write(102,FMT) lonx, latx, rate_u(i)*from_m_to_mm 
            	         write(103,FMT) lonx, latx, rate_n(i)*from_m_to_mm  
     	        	 write(104,FMT) lonx, latx, rate_s(i)/slc_ave 
       		enddo
!
! --- Elsewhere -------------------------  Global
            elseif (lat(i)/=-90.and.lat(i)/=90.) then 
       		latx=lat(i) 
		lonx=lon(i) 
			 write(101,FMT) lonx, latx, rate_s(i)*from_m_to_mm  
			 write(102,FMT) lonx, latx, rate_u(i)*from_m_to_mm  
			 write(103,FMT) lonx, latx, rate_n(i)*from_m_to_mm  
     	        	 write(104,FMT) lonx, latx, rate_s(i)/slc_ave 
      	    Endif	       
!
      END DO
! 
 deallocate( lon, lat, rate_s, rate_u, rate_n )
!
 Close(101) ; Close(102) ; Close(103) ; Close(104) 
!
 ENDIF ! --------------------------------- Global  
!
!
! //////////////////////////////////////////////////////////////////  / Regional Greenland
! //////////////////////////////////////////////////////////////////  / Regional Greenland  
!
!  For regional analyses, I need to input the spherical harmonics at
!  pixels (note: a NEW set of pixels is used here for solving the SLE,
!  generated according to the spatial resolution set in config.dat..). 
!
! //////////////////////////////////////////////////////////////////  / Regional Greenland 
! //////////////////////////////////////////////////////////////////  / Regional Greenland   
!
  IF(ELAREB_GR==1) then 
!
!
! Tegmark resolution according to the grid spacing 
  RES_MAP=AINT(0.4*(6.371E3/RES_REG_GRC))
! write(*,*) res_map				
!
! Number of pixels corresponding to this resolution 
  NPP=2*RES_MAP*(RES_MAP-1)*20+12
! write(*,*) NPP
!
! Allocate memory space
!
  ALLOCATE( LONXX(NPP), LATXX(NPP), IDX(NPP) )
!
! # Pixelization of the sphere 			
    call FINDPX(RES_MAP, NPP, LONXX, LATXX)	
!
! Pixels within the Greenland rectangle 
  nin=0
  DO I=1, NPP 
  	LOGI1=(lonxx(i).ge.lonmin_gr.and.lonxx(i).le.lonmax_gr)
  	LOGI2=(latxx(i).ge.latmin_gr.and.latxx(i).le.latmax_gr)   
  	if( LOGI1.AND.LOGI2) THEN 
  	   nin=nin+1
  	   IDX(NIN)=i
  	endif  
  END DO 
!
  Write(*,*) '    - There are ', nin, ' pixels with size ', & 
              RES_REG_GRC, ' km in Greenland' 
!
!
! # Extract the regional pixels
!
  allocate( lon(nin), lat(nin) )
  lon = lonxx( idx(1:nin) )
  lat = latxx( idx(1:nin) )
  deallocate( lonxx, latxx, idx )
  allocate( rate_s(nin), rate_u(nin), rate_n(nin) )
!
  RATE_S   = 0.
  RATE_U   = 0.
  RATE_N   = 0.
!
!
! --- Call HARMO once outside the parallel loop to initialize 
!     the SAVEd arrays
!
  call HARMO(LMAX,lon(1),lat(1),ARMOY)  
!
!
!$OMP PARALLEL DO DEFAULT(PRIVATE) &
!$OMP    SHARED(NIN,LON,LAT,RATE_S,RATE_U,RATE_N,DM,CS,CU,CN) &
!$OMP    SCHEDULE(GUIDED)
  DO I=1, NIN
       IF( .NOT. ALLOCATED(ARMOY) ) ALLOCATE(ARMOY(JMAX))
       call HARMO(LMAX, lon(i), lat(i), ARMOY)
       DO J=1, JMAX  
          rate_s(i) = rate_s(i)  + dm(j)*real(cs(j,nn)*ARMOY(J))   
          rate_u(i) = rate_u(i)  + dm(j)*real(cu(j,nn)*ARMOY(J)) 
          rate_n(i) = rate_n(i)  + dm(j)*real(cn(j,nn)*ARMOY(J)) 
       END DO
  END DO
!$OMP END PARALLEL DO  

! 
! --- Note that north & south poles are "never" included 
!     in the regional map of Greenland... Rev GS 25/04/11.
!
!
   open(101,file= 'smap_gree_reg.dat',  status='unknown') 
   open(102,file= 'umap_gree_reg.dat',  status='unknown') 
   open(103,file= 'nmap_gree_reg.dat',  status='unknown') 
   do i=1,nin
          write(101,FMT) lon(i), lat(i), rate_s(i)*from_m_to_mm  
          write(102,FMT) lon(i), lat(i), rate_u(i)*from_m_to_mm 
          write(103,FMT) lon(i), lat(i), rate_n(i)*from_m_to_mm 
   end do

   close(101); close(102); close(103) 
!
!
   DEALLOCATE( LON, LAT, RATE_S, RATE_U, RATE_N )
!
!
   ENDIF   ! / Regional Greenland   
!
!
!

!
! //////////////////////////////////////////////////////////////////  / Regional Antarctica 
! //////////////////////////////////////////////////////////////////  / Regional Antarctica  
!
!  For regional analyses, I need to input the spherical harmonics at
!  pixels (note: a NEW set of pixels is used here for solving the SLE,
!  generated according to the spatial resolution set in config.dat..). 
!
! //////////////////////////////////////////////////////////////////  / Regional Antarctica  
! //////////////////////////////////////////////////////////////////  / Regional Antarctica    
!
  IF(ELAREB_AR==1) then 
!
! Tegmark resolution according to the grid spacing 
  RES_MAP=AINT(0.4*(6.371E3/RES_REG_ANC))
! write(*,*) res_map				
!
! Number of pixels corresponding to this resolution 
  NPP=2*RES_MAP*(RES_MAP-1)*20+12
! write(*,*) NPP
  allocate( lonxx(npp), latxx(npp), idx(npp) )
!
! # Pixelization of the sphere 			
    call FINDPX(RES_MAP, NPP, LONXX, LATXX)	
!
! Pixels within the Antarctic "cap" 
  nin=0
  DO I=1, NPP 
  	LOGI2=(latxx(i).ge.latmin_ar.and.latxx(i).le.latmax_ar)   
  	if(LOGI2) then
  	   nin=nin+1  
  	   idx(nin)=i
  	endif
  end do
!
  Write(*,*) '    - There are ', nin, ' pixels with size ', RES_REG_ANC, ' km in Antarctica' 
!
!
!
! # Extract the regional pixels
!
  allocate( lon(nin), lat(nin) )
  lon = lonxx( idx(1:nin) )
  lat = latxx( idx(1:nin) )
  deallocate( lonxx, latxx, idx )
  allocate( rate_s(nin), rate_u(nin), rate_n(nin) )
!
  RATE_S   = 0.
  RATE_U   = 0.
  RATE_N   = 0.
!
!
!
!
! --- Call HARMO once outside the parallel loop to initialize 
!     the SAVEd arrays
!
  call HARMO(LMAX,lon(1),lat(1),ARMOY)  
!
!
!$OMP PARALLEL DO DEFAULT(PRIVATE) &
!$OMP    SHARED(NIN,LON,LAT,RATE_S,RATE_U,RATE_N,DM,CS,CU,CN) &
!$OMP    SCHEDULE(GUIDED)
  DO I=1, NIN
       IF( .NOT. ALLOCATED(ARMOY) ) ALLOCATE(ARMOY(JMAX))
       call HARMO(LMAX, lon(i), lat(i), ARMOY)
       DO J=1, JMAX  
          rate_s(i) = rate_s(i)  + dm(j)*real(cs(j,nn)*ARMOY(J))   
          rate_u(i) = rate_u(i)  + dm(j)*real(cu(j,nn)*ARMOY(J)) 
          rate_n(i) = rate_n(i)  + dm(j)*real(cn(j,nn)*ARMOY(J)) 
       END DO
  END DO
!$OMP END PARALLEL DO
!
!
  open(101,file= 'smap_anta_reg.dat',  status='unknown') 
  open(102,file= 'umap_anta_reg.dat',  status='unknown') 
  open(103,file= 'nmap_anta_reg.dat',  status='unknown') 
!
  do i=1,nin
      latx=lat(i) 
      lonx=lon(i) 
!
! --- South pole  ------------------------  Global
            If (lat(i)==-90.) then  
       		latx=lat(i)+eps
       		do k=1, 4 
                	lonx=(k-1)*90.   		
     	        	write(101,FMT) lonx, latx, rate_s(i)*from_m_to_mm 
	        	    write(102,FMT) lonx, latx, rate_u(i)*from_m_to_mm  
	        	    write(103,FMT) lonx, latx, rate_n(i)*from_m_to_mm  
       		enddo
!
! --- Elsewhere -------------------------  Global
            elseif (lat(i)/=-90.and.lat(i)/=90.) then 
      		latx=lat(i) 
      		lonx=lon(i) 
			 write(101,FMT) lonx, latx, rate_s(i)*from_m_to_mm  
			 write(102,FMT) lonx, latx, rate_u(i)*from_m_to_mm  
			 write(103,FMT) lonx, latx, rate_n(i)*from_m_to_mm  
      	    Endif	
!
  end do
!
   close(101); close(102); close(103) 
!
!
   deallocate( lon,lat,rate_s,rate_u,rate_n )
!
   write(*,*) '    - Number of pixels in Antarctica: ', nin
!
   close(101); close(102); close(103) 
!
   ENDIF   ! / Regional Antarctica   
!
! --- Release memory
!
 DEALLOCATE( MM, DM )
 DEALLOCATE( ALF )
 DEALLOCATE( ARMOY, OC )
 DEALLOCATE( CS,  CU, CN  )
 DEALLOCATE( LONG_TABLE  )
 DEALLOCATE( Y )
! 
 END PROGRAM ELA_REB
!
!
!
!
!
!
