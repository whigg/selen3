!
! This is program "GMAPS.F90" 
!
! Last modified GS 04-11-2008 [Intel port]
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! Reviewed GS June 2010 - FREE AIR & SS ANOMALIES
! Reviewed GS July 2010 - Implementation of the "LOAD" 
! *** Revised GS July 2010 - g95 - Double precision implementation
! *** Revised GS March 011 - for the \Phi/\gamma implementation...
! === Revised GS & FC May 21 2011 - DELTA parameter in ice history     
! === Revised DM Nov 30 2011 - Dynamic memory & loop parallelization execution
! Feb 2012: Implementation of the numerical derivative "in the future"    
! July 15, 2012; regular grid for the global maps at presen time...
! Aug 5, 2012; Equivalent water height (understanding the Chambers et al. 2010 results)
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
! 	- shs.bin
!	- shu.bin
!	- shn.bin
!	- shg.bin
!	- shz.bin 
!	- shf.bin 
!	- shice.dat 
! 	- sh.bin
! 	- px-table.dat
!       - anchor.tmp      [NEW]
!
! Output files:
! 	- sdotmap.dat 
!	- udotmap.dat
!	- ndotmap.dat
!	- gdotmap.dat
!       - wdotmap.dat     ! ---- New as of Aug 5, 2012: Equivalent Water Height 
!	- fadotmap.dat
!	- ssdotmap.dat
!       - loadimap.dat    ! load (ice)
!       - loadomap.dat    ! load (ocean)  
!       - loadtmap.dat    ! load (ice + ocean) 
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!INCLUDE "harmonics.f90"
 PROGRAM GMAPS
 IMPLICIT NONE 
 INCLUDE "data.inc"
!
 CHARACTER*30 HEADER
 CHARACTER*13, PARAMETER :: FMT='(3(f10.6,1X))'
 CHARACTER*11, PARAMETER :: UF='UNFORMATTED', F='FORMATTED', UN='UNKNOWN' 
 REAL*8, PARAMETER :: EPS = 0.01
!
 INTEGER :: NANCH
 INTEGER I, J, K, DOM, MJ, IT1, IT2, NBUF
 REAL*8 :: DT
 INTEGER, ALLOCATABLE :: MM(:), ANC(:), DM(:)  
!
 REAL*8 LONX, LATX
 REAL*8, ALLOCATABLE :: LON(:), LAT(:)
 REAL*8, ALLOCATABLE :: ALF(:,:)
 REAL*8, ALLOCATABLE :: RATE_S(:), RATE_U(:), RATE_N(:), RATE_FA(:), RATE_SS(:)
 REAL*8, ALLOCATABLE :: RATE_LOI(:), RATE_LOO(:), RATE_LOT(:), RATE_G(:), RATE_W(:)
! 
 REAL*8, PARAMETER :: DDELTA=DELTA 
!
 COMPLEX*16, ALLOCATABLE ::  CS(:,:),  CU(:,:),  & 
                             CN(:,:),  CFA(:,:), &
                             CSS(:,:), CLI(:,:), &
                             CLO(:,:), CLT(:,:), & 
	                         LONG_TABLE(:,:), CG(:,:), CW(:,:) 
!
!
! --- Data for the computations on a 1 deg x 1 deg grid - July 2012
 Integer kh, jj 
 Integer, parameter :: delta_grid=1 
 Integer, parameter :: NGRID=180/delta_grid
 REAL*8, ALLOCATABLE :: LONGC(:), LATC(:)
 COMPLEX*16, ALLOCATABLE :: Y(:,:)	   
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
! --- Allocating memory
!     Updated February 2012 with "nn+1"
!
  if (gmaps_grid==1)   nbuf=np               
  if (gmaps_grid==2)   nbuf=2*ngrid*ngrid
!
  allocate( mm(jmax) )
  allocate( dm(jmax) )
  if (gmaps_grid==1)  allocate( anc(np) )
  if (gmaps_grid==1)  allocate( lon(np), lat(np) )
  if (gmaps_grid==1)  allocate( alf(jmax,nanch) )
  if (gmaps_grid==1)  allocate( long_table(0:lmax,np) )
  if (glob_i_s == 1)  allocate( cs(jmax,0:nn+1),  rate_s(nbuf) )
  if (glob_i_u == 1)  allocate( cu(jmax,0:nn+1),  rate_u(nbuf) )
  if (glob_i_n == 1)  allocate( cn(jmax,0:nn+1),  rate_n(nbuf) )
  if (glob_i_g == 1)  allocate( cg(jmax,0:nn+1),  rate_g(nbuf) )
  if (glob_i_w == 1)  allocate( cw(jmax,0:nn+1),  rate_w(nbuf) )
  if (glob_i_fa== 1)  allocate( cfa(jmax,0:nn+1), rate_fa(nbuf) )
  if (glob_i_ss== 1)  allocate( css(jmax,0:nn+1), rate_ss(nbuf) )
  if (glob_i_loi== 1) allocate( cli(jmax,0:nn+1), rate_loi(nbuf) )
  if (glob_i_loo== 1) allocate( clo(jmax,0:nn+1), rate_loo(nbuf) )
  if (glob_i_lot== 1) allocate( clt(jmax,0:nn+1), rate_lot(nbuf) )
  if (gmaps_grid==2)  allocate( Y(JMAX,2*NGRID*NGRID) )
  if (gmaps_grid==2)  allocate( longc(2*NGRID*NGRID), latc(2*NGRID*NGRID ) )
!
!
! --- Pre-computing the degree 'm' corresponding to 'J'
!
 do j=1, jmax 
        mm(j)=mj(j)
	dm(j)=2-dom(j) 
 enddo	
!
!
 IF(GMAPS_GRID==1) THEN ! +++++++++++++++++++++++++++++++++++
!
!
!---- Reading the pixels table file... 
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
!
!---- Reading the table of spherical harmonics from <<sh.bin>>
!
 Write(*,*) "    - Reading ALFs & TRIGs from sh.bin"
 open(3,file='sh.bin',status='unknown',form='unformatted')  
 	read(3)ALF
 	read(3)LONG_TABLE
 Close(3) 
 DO J=1, JMAX 
	ALF(J,:)=ALF(J,:)*DM(J)
 ENDDO
!
 ENDIF ! +++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!
!---- Computing S-dot, U-dot, N-dot and other at pixels...  
!
!
 if(glob_i_s==1) then 
 	    open (101,file= 'sdotmap.dat',  status='unknown') 
 	    open (1,file='shs.bin',       status='unknown',form=uf) 
 	    read (1) CS ; close(1) 
 endif
 if(glob_i_u==1) then 
 	    open(102,file= 'udotmap.dat',  status='unknown') 
 	    open (1,file='shu.bin',       status='unknown',form=uf) 
 	    read (1) CU ; close(1) 
 endif
 if(glob_i_n==1) then 
 	    open(103,file= 'ndotmap.dat',  status='unknown') 
 	    open (1,file='shn.bin',       status='unknown',form=uf) 
 	    read (1) CN ; close(1) 
 endif
 if(glob_i_g==1) then 
 	    open(109,file= 'gdotmap.dat',  status='unknown') 
 	    open (1,file='shg.bin',       status='unknown',form=uf) 
 	    read (1) CG ; close(1) 
 endif
 if(glob_i_w==1) then 
 	    open(111,file= 'wdotmap.dat',  status='unknown') 
 	    open (1,file='shw.bin',       status='unknown',form=uf) 
 	    read (1) CW ; close(1) 
 endif
 if(glob_i_fa==1) then 
 	    open(104,file= 'fadotmap.dat', status='unknown') 
 	    open (1,file='shfa.bin',      status='unknown',form=uf) 
 	    read (1) CFA ; close(1) 
 endif
 if(glob_i_ss==1) then 
 	    open(105,file= 'ssdotmap.dat', status='unknown') 
 	    open (1,file='shss.bin',      status='unknown',form=uf) 
 	    read (1) CSS ; close(1) 
 endif
 if(glob_i_loi==1) then 
 	    open(106,file= 'lidotimap.dat',status='unknown')
 	    open (1,file='shload_ice.bin',status='unknown',form=uf) 
 	    read (1) CLI ; close(1) 
 endif
 if(glob_i_loo==1) then 
 	    open(107,file= 'lodotimap.dat',status='unknown')
 	    open (1,file='shload_oce.bin',status='unknown',form=uf) 
 	    read (1) CLO ; close(1) 
 endif
 if(glob_i_lot==1) then 
 	    open(108,file= 'ltdotimap.dat',status='unknown')
 	    open (1,file='shload.bin',    status='unknown',form=uf) 
 	    read (1) CLT ; close(1) 
 endif
!
!
! 
 IF(GMAPS_GRID==1) THEN ! +++++++++++++++++++++++++++++++++++
!
!
!
! IF( IDER==1 ) THEN
!    IT1 = NN-1
!    IT2 = NN
!    DT  = DDELTA
! ENDIF
! IF( IDER==2 ) THEN
!    IT1 = NN-1
!    IT2 = NN+1
!    DT  = 2.D0*DDELTA
! ENDIF
 
 IF( IDER==1 ) THEN
    IT1 = NN-1
    IT2 = NN
    DT  = DDELTA
 ENDIF
 IF( IDER==2 ) THEN
    IT1 = NN-1
    IT2 = NN+1
    DT  = 2.D0*DDELTA
 ENDIF 
    
!
!
 Write(*,*) "    - Computing rates at ", np, " pixels "
!
!
 if(glob_i_s  ==1) then
    Write(*,*) "      --> Sealevel change"
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J) &
!$OMP    SHARED(RATE_S,ALF,ANC,CS,LONG_TABLE,MM,IT1,IT2,DT) &
!$OMP    SCHEDULE(GUIDED)
    do i=1,np
       rate_s(i) = 0.
       do j=1,jmax
           rate_s(i) = rate_s(i) + ALF(j,anc(i))*REAL(((cs(j,it2)  -cs(j,it1))/dt)*LONG_TABLE(MM(J),I))   
       end do
    end do
!$OMP END PARALLEL DO
 end if
!
 if(glob_i_u  ==1) then
     Write(*,*) "      --> Vertical displacement"
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J) &
!$OMP    SHARED(RATE_U,ALF,ANC,CU,LONG_TABLE,MM,IT1,IT2,DT) &
!$OMP    SCHEDULE(GUIDED)
    do i=1,np
       rate_u(i) = 0.
       do j=1,jmax
           rate_u(i) = rate_u(i) + ALF(j,anc(i))*REAL(((cu(j,it2)  -cu(j,it1))/dt)*LONG_TABLE(MM(J),I))   
       end do
    end do
!$OMP END PARALLEL DO
 end if
!
 if(glob_i_n  ==1) then
     Write(*,*) "      --> Geoid height variation"
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J) &
!$OMP    SHARED(RATE_N,ALF,ANC,CN,LONG_TABLE,MM,IT1,IT2,DT) &
!$OMP    SCHEDULE(GUIDED)
    do i=1,np
       rate_n(i) = 0.
       do j=1,jmax
           rate_n(i) = rate_n(i) + ALF(j,anc(i))*REAL(((cn(j,it2)  -cn(j,it1))/dt)*LONG_TABLE(MM(J),I))   
       end do
    end do
!$OMP END PARALLEL DO
 end if
!
 if(glob_i_g  ==1) then
     Write(*,*) "      --> Normalized gravity potential"
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J) &
!$OMP    SHARED(RATE_G,ALF,ANC,CG,LONG_TABLE,MM,IT1,IT2,DT) &
!$OMP    SCHEDULE(GUIDED)
    do i=1,np
       rate_g(i) = 0.
       do j=1,jmax
           rate_g(i) = rate_g(i) + ALF(j,anc(i))*REAL(((cg(j,it2)  -cg(j,it1))/dt)*LONG_TABLE(MM(J),I))   
       end do
    end do
!$OMP END PARALLEL DO
 end if
!
 if(glob_i_w  ==1) then
     Write(*,*) "      --> Equivalent water heigth"
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J) &
!$OMP    SHARED(RATE_W,ALF,ANC,CW,LONG_TABLE,MM,IT1,IT2,DT) &
!$OMP    SCHEDULE(GUIDED)
    do i=1,np
       rate_w(i) = 0.
       do j=1,jmax
           rate_w(i) = rate_w(i) + ALF(j,anc(i))*REAL(((cw(j,it2)  -cw(j,it1))/dt)*LONG_TABLE(MM(J),I))   
       end do
    end do
!$OMP END PARALLEL DO
 end if 
!
 if(glob_i_fa  ==1) then
     Write(*,*) "      --> Free air gravity anomaly"
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J) &
!$OMP    SHARED(RATE_fa,ALF,ANC,Cfa,LONG_TABLE,MM,IT1,IT2,DT) &
!$OMP    SCHEDULE(GUIDED)
    do i=1,np
       rate_fa(i) = 0.
       do j=1,jmax
           rate_fa(i) = rate_fa(i) + ALF(j,anc(i))*REAL(((cfa(j,it2)  -cfa(j,it1))/dt)*LONG_TABLE(MM(J),I))   
       end do
    end do
!$OMP END PARALLEL DO
 end if
!
 if(glob_i_ss  ==1) then
     Write(*,*) "      --> Solid surface gravity anomaly"
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J) &
!$OMP    SHARED(RATE_SS,ALF,ANC,CSS,LONG_TABLE,MM,IT1,IT2,DT) &
!$OMP    SCHEDULE(GUIDED)
    do i=1,np
       rate_ss(i) = 0.
       do j=1,jmax
           rate_ss(i) = rate_ss(i) + ALF(j,anc(i))*REAL(((css(j,it2)  -css(j,it1))/dt)*LONG_TABLE(MM(J),I))   
       end do
    end do
!$OMP END PARALLEL DO
 end if
!
 if(glob_i_loi  ==1) then
     Write(*,*) "      --> Ice load"
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J) &
!$OMP    SHARED(RATE_LOI,ALF,ANC,CLI,LONG_TABLE,MM,IT1,IT2,DT) &
!$OMP    SCHEDULE(GUIDED)
    do i=1,np
       rate_loi(i) = 0.
       do j=1,jmax
           rate_loi(i) = rate_loi(i) + ALF(j,anc(i))*REAL(((cli(j,it2)  -cli(j,it1))/dt)*LONG_TABLE(MM(J),I))   
       end do
    end do
!$OMP END PARALLEL DO
 end if
!
 if(glob_i_loo  ==1) then
     Write(*,*) "      --> Ocean load"
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J) &
!$OMP    SHARED(RATE_LOO,ALF,ANC,CLO,LONG_TABLE,MM,IT1,IT2,DT) &
!$OMP    SCHEDULE(GUIDED)
    do i=1,np
       rate_loo(i) = 0.
       do j=1,jmax
           rate_loo(i) = rate_loo(i) + ALF(j,anc(i))*REAL(((clo(j,it2)  -clo(j,it1))/dt)*LONG_TABLE(MM(J),I))   
       end do
    end do
!$OMP END PARALLEL DO
 end if
!
 if(glob_i_lot  ==1) then
     Write(*,*) "      --> Total load"
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J) &
!$OMP    SHARED(RATE_LOT,ALF,ANC,CLT,LONG_TABLE,MM,IT1,IT2,DT) &
!$OMP    SCHEDULE(GUIDED)
    do i=1,np
       rate_lot(i) = 0.
       do j=1,jmax
           rate_lot(i) = rate_lot(i) + ALF(j,anc(i))*REAL(((clt(j,nn)  -clt(j,it1))/dt)*LONG_TABLE(MM(J),I))   
       end do
    end do
!$OMP END PARALLEL DO
 end if
!
!
!
! --- Write outputs
!
 do i=1,np
!
!
! --- South pole
       If     (lat(i)==-90.) then  
       		latx=lat(i)+eps
       		do k=1, 4 
			lonx=(k-1)*90.   		
            		 if(glob_i_s   ==1) write(101,FMT) lonx, latx, rate_s(i) 
            		 if(glob_i_u   ==1) write(102,FMT) lonx, latx, rate_u(i) 
            		 if(glob_i_n   ==1) write(103,FMT) lonx, latx, rate_n(i) 
            		 if(glob_i_g   ==1) write(109,FMT) lonx, latx, rate_g(i)      ! r 
            		 if(glob_i_w   ==1) write(111,FMT) lonx, latx, rate_w(i)      ! r 
            		 if(glob_i_fa  ==1) write(104,FMT) lonx, latx, rate_fa(i) 
            		 if(glob_i_ss  ==1) write(105,FMT) lonx, latx, rate_ss(i) 
            		 if(glob_i_loi ==1) write(106,FMT) lonx, latx, rate_loi(i)
            		 if(glob_i_loo ==1) write(107,FMT) lonx, latx, rate_loo(i) 
            		 if(glob_i_lot ==1) write(108,FMT) lonx, latx, rate_lot(i) 
       		enddo
!
! --- North pole
       elseif (lat(i)==+90.) then   
       		latx=lat(i)-eps
       		do k=1, 4 
			lonx=(k-1)*90.   			
            	         if(glob_i_s   ==1) write(101,FMT) lonx, latx, rate_s(i) 
            	         if(glob_i_u   ==1) write(102,FMT) lonx, latx, rate_u(i) 
            	         if(glob_i_n   ==1) write(103,FMT) lonx, latx, rate_n(i) 
             		     if(glob_i_g   ==1) write(109,FMT) lonx, latx, rate_g(i)	   ! r 
                 		 if(glob_i_w   ==1) write(111,FMT) lonx, latx, rate_w(i)       ! r 
            	         if(glob_i_fa  ==1) write(104,FMT) lonx, latx, rate_fa(i) 
            	         if(glob_i_ss  ==1) write(105,FMT) lonx, latx, rate_ss(i) 
            	         if(glob_i_loi ==1) write(106,FMT) lonx, latx, rate_loi(i)
            	         if(glob_i_loo ==1) write(107,FMT) lonx, latx, rate_loo(i) 
            	         if(glob_i_lot ==1) write(108,FMT) lonx, latx, rate_lot(i) 
       		enddo
!
! --- Elsewhere 
       elseif (lat(i)/=-90.and.lat(i)/=90.) then 
       		latx=lat(i) 
		lonx=lon(i) 
            	         if(glob_i_s   ==1) write(101,FMT) lonx, latx, rate_s(i) 
            	         if(glob_i_u   ==1) write(102,FMT) lonx, latx, rate_u(i) 
            	         if(glob_i_n   ==1) write(103,FMT) lonx, latx, rate_n(i) 
            		     if(glob_i_g   ==1) write(109,FMT) lonx, latx, rate_g(i)	   ! r
                		 if(glob_i_w   ==1) write(111,FMT) lonx, latx, rate_w(i)       ! r 
            	         if(glob_i_fa  ==1) write(104,FMT) lonx, latx, rate_fa(i) 
            	         if(glob_i_ss  ==1) write(105,FMT) lonx, latx, rate_ss(i) 
            	         if(glob_i_loi ==1) write(106,FMT) lonx, latx, rate_loi(i)
            	         if(glob_i_loo ==1) write(107,FMT) lonx, latx, rate_loo(i)
            	         if(glob_i_lot ==1) write(108,FMT) lonx, latx, rate_lot(i) 
       Endif	       
!
end do
!
 do k=101, 109
     close(k)
 enddo
!
 ENDIF ! On gmaps_grid =1 ------------------------------------------------------
!
!
!
!
! IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - 
! IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - 
!
! 
 IF(GMAPS_GRID==2) THEN ! +++++++++++++++++++++++++++++++++++
!
! --- Define grid nodes
!
 k = 0
!
! From the south pole (-90^) to the north pole (+90^)... 
! [modified July 15, 2012]
!
 do i=1,2*NGRID
 do j=1,NGRID
   k=k+1
   longc(k) = float(delta_grid)/2.*(2*i-1)
   latc(k)  = - 90.+ float(delta_grid)/2.*(2*j-1)
 end do
 end do
!
!
! --- Choose derivation type
!
 IF( IDER==1 ) THEN
    IT1 = NN-1
    IT2 = NN
    DT  = DDELTA
 ENDIF
 IF( IDER==2 ) THEN
    IT1 = NN-1
    IT2 = NN+1
    DT  = 2.D0*DDELTA
 ENDIF     
!
!
! --- Compute SHs
!
 do k=1,2*ngrid*ngrid
    call harmo( lmax, longc(k), latc(k), y(:,1) )
 end do
!
 Write(*,*) "    - Computing rates on a 1 deg x 1 deg grid "
!
!
!
 if(glob_i_s  ==1) then
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J) &
!$OMP    SHARED(RATE_S,Y,CS,IT1,IT2,DT) &
!$OMP    SCHEDULE(GUIDED)
    do i=1,2*ngrid*ngrid
       rate_s(i) = 0.
       do j=1,jmax
           rate_s(i) = rate_s(i) + (2.0-dom(j))*REAL((( cs(j,it2) -  cs(j,it1))/dt)*Y(J,I)) 
       end do
    end do
!$OMP END PARALLEL DO
 end if
!
 if(glob_i_u  ==1) then
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J) &
!$OMP    SHARED(RATE_U,Y,CU,IT1,IT2,DT) &
!$OMP    SCHEDULE(GUIDED)
    do i=1,2*ngrid*ngrid
       rate_u(i) = 0.
       do j=1,jmax
           rate_u(i) = rate_u(i) + (2.0-dom(j))*REAL((( cu(j,it2) -  cu(j,it1))/dt)*Y(J,I)) 
       end do
    end do
!$OMP END PARALLEL DO
 end if
!
 if(glob_i_n  ==1) then
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J) &
!$OMP    SHARED(RATE_N,Y,CN,IT1,IT2,DT) &
!$OMP    SCHEDULE(GUIDED)
    do i=1,2*ngrid*ngrid
       rate_n(i) = 0.
       do j=1,jmax
           rate_n(i) = rate_n(i) + (2.0-dom(j))*REAL((( cn(j,it2) -  cn(j,it1))/dt)*Y(J,I)) 
       end do
    end do
!$OMP END PARALLEL DO
 end if
!
 if(glob_i_g  ==1) then
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J) &
!$OMP    SHARED(RATE_G,Y,CG,IT1,IT2,DT) &
!$OMP    SCHEDULE(GUIDED)
    do i=1,2*ngrid*ngrid
       rate_g(i) = 0.
       do j=1,jmax
           rate_g(i) = rate_g(i) + (2.0-dom(j))*REAL((( cg(j,it2) -  cg(j,it1))/dt)*Y(J,I)) 
       end do
    end do
!$OMP END PARALLEL DO
 end if
!
 if(glob_i_w  ==1) then
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J) &
!$OMP    SHARED(RATE_W,Y,CW,IT1,IT2,DT) &
!$OMP    SCHEDULE(GUIDED)
    do i=1,2*ngrid*ngrid
       rate_w(i) = 0.
       do j=1,jmax
           rate_w(i) = rate_w(i) + (2.0-dom(j))*REAL((( cw(j,it2) -  cw(j,it1))/dt)*Y(J,I)) 
       end do
    end do
!$OMP END PARALLEL DO
 end if
!
 if(glob_i_fa ==1) then
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J) &
!$OMP    SHARED(RATE_FA,Y,CFA,IT1,IT2,DT) &
!$OMP    SCHEDULE(GUIDED)
    do i=1,2*ngrid*ngrid
       rate_fa(i) = 0.
       do j=1,jmax
           rate_fa(i) = rate_fa(i) + (2.0-dom(j))*REAL((( cfa(j,it2) -  cfa(j,it1))/dt)*Y(J,I)) 
       end do
    end do
!$OMP END PARALLEL DO
 end if
!
 if(glob_i_ss ==1) then
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J) &
!$OMP    SHARED(RATE_SS,Y,CSS,IT1,IT2,DT) &
!$OMP    SCHEDULE(GUIDED)
    do i=1,2*ngrid*ngrid
       rate_ss(i) = 0.
       do j=1,jmax
           rate_ss(i) = rate_ss(i) + (2.0-dom(j))*REAL((( css(j,it2) -  css(j,it1))/dt)*Y(J,I)) 
       end do
    end do
!$OMP END PARALLEL DO
 end if
!
 if(glob_i_loi==1) then
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J) &
!$OMP    SHARED(RATE_LOI,Y,CLI,IT1,IT2,DT) &
!$OMP    SCHEDULE(GUIDED)
    do i=1,2*ngrid*ngrid
       rate_loi(i) = 0.
       do j=1,jmax
           rate_loi(i) = rate_loi(i) + (2.0-dom(j))*REAL((( cli(j,it2) -  cli(j,it1))/dt)*Y(J,I)) 
       end do
    end do
!$OMP END PARALLEL DO
 end if
!
 if(glob_i_loo==1) then
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J) &
!$OMP    SHARED(RATE_LOO,Y,CLO,IT1,IT2,DT) &
!$OMP    SCHEDULE(GUIDED)
    do i=1,2*ngrid*ngrid
       rate_loo(i) = 0.
       do j=1,jmax
           rate_loo(i) = rate_loo(i) + (2.0-dom(j))*REAL((( clo(j,it2) -  clo(j,it1))/dt)*Y(J,I)) 
       end do
    end do
!$OMP END PARALLEL DO
 end if
!
 if(glob_i_lot==1) then
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J) &
!$OMP    SHARED(RATE_LOT,Y,CLT,IT1,IT2,DT) &
!$OMP    SCHEDULE(GUIDED)
    do i=1,2*ngrid*ngrid
       rate_lot(i) = 0.
       do j=1,jmax
           rate_lot(i) = rate_lot(i) + (2.0-dom(j))*REAL((( clt(j,it2) -  clt(j,it1))/dt)*Y(J,I)) 
       end do
    end do
!$OMP END PARALLEL DO
 end if
!
! --- Write outputs
!
!
 do i=1,2*ngrid*ngrid
      IF (glob_i_s==1)    write(101,FMT) longc(i), latc(i), rate_s(i) 
      IF (glob_i_u==1)    write(102,FMT) longc(i), latc(i), rate_u(i) 
      IF (glob_i_n==1)    write(103,FMT) longc(i), latc(i), rate_n(i) 
      IF (glob_i_g==1)    write(109,FMT) longc(i), latc(i), rate_g(i)   
      IF (glob_i_w==1)    write(111,FMT) longc(i), latc(i), rate_w(i)
      IF (glob_i_fa==1)   write(104,FMT) longc(i), latc(i), rate_fa(i) 
      IF (glob_i_ss==1)   write(105,FMT) longc(i), latc(i), rate_ss(i) 
      IF (glob_i_loi==1)  write(106,FMT) longc(i), latc(i), rate_loi(i)
      IF (glob_i_loo==1)  write(107,FMT) longc(i), latc(i), rate_loo(i) 
      IF (glob_i_lot==1)  write(108,FMT) longc(i), latc(i), rate_lot(i) 
 end do
!
!
 do k=101, 109  
     close(k)
 enddo
!
 ENDIF ! On gmaps_grid =2 ------------------------------------------------------
!
!
! IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - 
! IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - IN PROGRESS - 
!
!
!
! ---- Release memory
!
  deallocate( mm )
  deallocate( dm )
  if( gmaps_grid == 2 ) deallocate( y, longc, latc )
  if( gmaps_grid == 1 ) deallocate( anc, lon, lat, alf, long_table )
  if (glob_i_s == 1)  deallocate( cs,  rate_s   )
  if (glob_i_u == 1)  deallocate( cu,  rate_u   )
  if (glob_i_n == 1)  deallocate( cn,  rate_n   )
  if (glob_i_g == 1)  deallocate( cg,  rate_g   )
  if (glob_i_fa== 1)  deallocate( cfa, rate_fa  )
  if (glob_i_ss== 1)  deallocate( css, rate_ss  )
  if (glob_i_loi== 1) deallocate( cli, rate_loi )
  if (glob_i_loo== 1) deallocate( clo, rate_loo )
  if (glob_i_lot== 1) deallocate( clt, rate_lot )
!
!
!
 END PROGRAM GMAPS
!
!
!
!
!
!
