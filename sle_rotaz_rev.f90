!
! This is program "SLE_ROTAZ.F90" 
!
! GS 04-11-2008 "Intel port..."
! Modified by GS september 2008 for horizontal movements 
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
! Also reviewed in April 2010 for the implementation of degree 1 
! (remember the mode coupling issue...) 
! Revised May 26 GS for new ice SH  
! Updated GS June 2010 for Free Air & SS Gravity implementation 
! Updated GS July 2010 for the imode=5 option (see Gabriella) 
! Updated GS July 2010 for the Load implementation 
! November 2010: Copied from "SLE.F90" for the implementation 
! ... of rotational feedback  --- GS 2010 
! November 2010: Implementation of the rotational feedback  
! Updated GS March 2011 for the "PHI" gravity potential implementation
! GS Prague November 2011: Missing rotational term in U and N
!   Thir program is for "N" based on that for "S"
!   After some tests with IDA and discussion with Zdeneck Martinek 
! Revised DM November 2011 - parallel execution & dynamic memory
! Feb 2012: Implementation of the numerical derivative "in the future" 
! April 22, 2012: The eustatic vector dumped for later use (Urbino hospita)
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
! ============================================
! Solves the SLE by the pseudo-spectral method 
! ============================================
!
! Input files:
!       - px-table.dat 
!	- shof.dat 
!	- shice.dat
!	- ebs.dat
!	- ebu.dat
!	- ebn.dat
!       - ebv.dat 
!       - ebfa.dat  <<< new... 
! 	- sh.bin
!       - anchor.tmp  <<< new
!
! Output files:
!	- shu.bin 
!	- shn.bin 
!	- shs.bin 
!	- shg.bin     <<< new 2011  
!       - shz.bin 
!       - shv.bin 
!       - shfa.bin    <<< new... 
!       - shss.bin    <<< new... 
!       - shload*.bin <<< new 
!
!INCLUDE "harmonics.f90"
 PROGRAM SLE 
 IMPLICIT NONE 
 INCLUDE "data.inc"			    
!
 CHARACTER*12 HEADER 
 INTEGER :: NANCH
 INTEGER I, J, K, L, P, IJ, IS, LJ, MJ, LI, DOM, IND
 INTEGER J_INDEX
 INTEGER, ALLOCATABLE :: LL(:), MM(:), DM(:), ANC(:), WET(:)     
!
 REAL*8 RESH, IMSH
 REAL*4 RHOE, RHOI_O_RHOE_X3, RHOW_O_RHOE_X3, RHOI_O_RHOW 
 REAL*8, ALLOCATABLE :: ALF(:,:), LONP(:), LATP(:), X(:,:)         
 REAL*8, ALLOCATABLE :: AAVV(:), BBVV(:), ZROTVV(:), SROTVV(:)                           
!
 COMPLEX*16, ALLOCATABLE :: LONG_TABLE(:,:), OC(:)                 
 COMPLEX*16, ALLOCATABLE :: IIII(:,:), ZE(:,:)        
 COMPLEX*16, ALLOCATABLE :: AAAA(:,:), AAAA_MOD(:,:)               
 COMPLEX*16, ALLOCATABLE :: ZROT(:,:), ZROT_MOD(:,:) 
 COMPLEX*16, ALLOCATABLE :: SROT(:,:), SROT_MOD(:,:) 
 COMPLEX*16, ALLOCATABLE :: UROT(:,:), NROT(:,:)
 COMPLEX*16, ALLOCATABLE :: BBBB(:,:), BBBB_MOD(:,:)               
 COMPLEX*16, ALLOCATABLE :: HHHH(:,:), KKKK(:,:)                   
! 
!--------------------------------------------------
! Convolved Green functions (viscous & elastic) 
! 
 REAL*4, ALLOCATABLE :: BETAS(:,:),  ES(:),    &  ! Sea level change  
                        BETAN(:,:),  EN(:),    &  ! Geoid undulations 
                        BETAU(:,:),  EU(:),    &  ! Vertical displacemnt
                        BETAV(:,:),  EV(:),    &  ! Horizontal displacement  
                        BETAFA(:,:), EFA(:),   &  ! Free air (FA) gravity anomaly 
                        BETASS(:,:), ESS(:)       ! Solid-surface (SS) gravity anomaly                            
!
!--------------------------------------------------
! Degree 2 order 1 rotationally-induced sea level change (input)
!
 CHARACTER*22, PARAMETER :: FILE_SROT='srotaz_21.dat'
 CHARACTER*22, PARAMETER :: FILE_UROT='urotaz_21.dat'
 CHARACTER*22, PARAMETER :: FILE_NROT='nrotaz_21.dat'
!
 REAL*8, ALLOCATABLE :: KPRIME(:)                
!
!
!
! Output arrays & filenames 
! 
 CHARACTER*22 Out_Filename  
!
! - "Eustatic" Sea level change -------------------  ! New as of Apr 2012
 complex*16, allocatable :: SE(:,:)   
 CHARACTER*22, PARAMETER :: FILEE='she.bin'   
!
! - Sea level change ------------------------------
 complex*16, allocatable :: S(:,:), SNR(:,:)            
 CHARACTER*22, PARAMETER :: FILES='shs.bin'   
!
! - Geoid undulations -----------------------------
 complex*16, allocatable :: N(:,:)   
 CHARACTER*22, PARAMETER :: FILEN='shn.bin'   
!
! - Vertical displacemnt --------------------------
 complex*16, allocatable :: U(:,:)   
 CHARACTER*22, PARAMETER :: FILEU='shu.bin'   
!
! - Normalized Gravity potential ------------------ \frac{\Phi}{\gamma})
 complex*16, allocatable :: G(:,:)   
 CHARACTER*22, PARAMETER :: FILEG='shg.bin'   
!
! - Equivalent water height --------   ! New as of Aug 2012: NEW NEW NEW NEW NEW NEW ---------
 complex*16, allocatable :: W(:,:)   
 CHARACTER*22, PARAMETER :: FILEW='shw.bin'   
!
! - Horizontal displacement -----------------------
 complex*16, allocatable :: V(:,:)   
 CHARACTER*22, PARAMETER :: FILEV='shv.bin'  
!
! - Free air (FA) gravity anomaly -----------------
 complex*16, allocatable :: FA(:,:)   
 CHARACTER*22, PARAMETER :: FILEFA='shfa.bin'  
!
! - Solid-surface (SS) gravity anomaly ------------
 complex*16, allocatable :: SS(:,:)   
 CHARACTER*22, PARAMETER :: FILESS='shss.bin' 
!
! - Ice component of load  ------------------------
 complex*16, allocatable :: LO_I(:,:)
 CHARACTER*22, PARAMETER :: FILELI='shload_ice.bin'
!
! - Ocean component of load -----------------------  
 complex*16, allocatable :: LO_O(:,:)
 CHARACTER*22, PARAMETER :: FILELO='shload_oce.bin' 
!
! - Total load  -----------------------------------
 complex*16, allocatable :: LO_T(:,:)
 CHARACTER*22, PARAMETER :: FILELT='shload.bin'  
!
! - Sea level times the ocean function ------------  
 complex*16, allocatable :: Z(:,:,:)                      
 CHARACTER*22, PARAMETER :: FILEZ='shz.bin'  
!--------------------------------------------------
!
!
!
!
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! --- Getting info about the number of anchor pixels 
!
  open(1,file='anchor.tmp',status='old')
  read(1,*) nanch
  close(1)
!
! --- Allocate memory space (New as of Feb 2012: it was "nn")
!
 ALLOCATE( LL(JMAX), MM(JMAX), DM(JMAX), ANC(NP), WET(NP) )
 ALLOCATE( ALF(JMAX,NANCH), LONP(NP), LATP(NP), X(NP,0:NN+1) )
 ALLOCATE( AAVV(0:NN+1), BBVV(0:NN+1) ) 
 ALLOCATE( SROTVV(0:NN+1), ZROTVV(0:NN+1) )
 ALLOCATE( LONG_TABLE(0:LMAX,NP), OC(JMAX) )
 ALLOCATE( IIII(JMAX,0:NN+1), ZE(JMAX,0:NN+1), SE(JMAX,0:NN+1) )
 ALLOCATE( AAAA(JMAX,0:NN+1), AAAA_MOD(JMAX,0:NN+1) )
 ALLOCATE( BBBB(JMAX,0:NN+1), BBBB_MOD(JMAX,0:NN+1) )
 ALLOCATE( SROT(JMAX,0:NN+1), SROT_MOD(JMAX,0:NN+1) )
 ALLOCATE( ZROT(JMAX,0:NN+1), ZROT_MOD(JMAX,0:NN+1) )
 ALLOCATE( UROT(JMAX,0:NN+1), NROT(JMAX,0:NN+1) )
 ALLOCATE( HHHH(JMAX,0:NN+1), KKKK(JMAX,0:NN+1) )
 ALLOCATE( BETAS (0:LMAX,0:NN+1), ES (0:LMAX) )
 ALLOCATE( BETAN (0:LMAX,0:NN+1), EN (0:LMAX) )
 ALLOCATE( BETAU (0:LMAX,0:NN+1), EU (0:LMAX) )
 ALLOCATE( BETAV (0:LMAX,0:NN+1), EV (0:LMAX) )     
 ALLOCATE( BETAFA(0:LMAX,0:NN+1), EFA(0:LMAX) )     
 ALLOCATE( BETASS(0:LMAX,0:NN+1), ESS(0:LMAX) )     
 ALLOCATE( S(JMAX,0:NN+1),   N(JMAX,0:NN+1), U(JMAX,0:NN+1), G(JMAX,0:NN+1) )
 ALLOCATE( SNR(JMAX,0:NN+1), V(JMAX,0:NN+1), FA(JMAX,0:NN+1), SS(JMAX,0:NN+1) )
 ALLOCATE( LO_I(JMAX,0:NN+1), LO_O(JMAX,0:NN+1), LO_T(JMAX,0:NN+1) )
 ALLOCATE( Z(JMAX,0:NN+1,0:SMAX), W(JMAX,0:NN+1) )
 ALLOCATE( KPRIME(1:LMAX) )
!
!
!
!
! ********************************************************
! ********************************************************
! ********************************************************
!
! --- Extracting the average Earth density from 
!     the TABOO or ALMA log files... GS July 09
!
 IF(CDE.GE. 0) IND=1 
 IF(CDE.EQ.-1) IND=2  
 IF(CDE.GE. 0.or.CDE==-1) CALL AVERAGE_EARTH_DENSITY(IND, RHOE)
! 
 IF(CDE==-2) RHOE=5511.57
! 
 RHOI_O_RHOE_X3 = 3.*RHOI/RHOE 
 RHOW_O_RHOE_X3 = 3.*RHOW/RHOE 
 RHOI_O_RHOW    =    RHOI/RHOW 
!
!
!
! --- Pre-computing 'l' and 'm' corresponding to degree 'J'
!
	do j=1, jmax 
		mm(j)=mj(j) 
		ll(j)=lj(j)
		dm(j)=2-dom(j)
	enddo	
!
!##################################################################  new
!##################################################################  new
!
! --- Reading the degree 2 order 1 rotationally-induced sealevel 
! 
 	Write(*,*) '    - Reading the degree 2 order 1 sea level'
 	open(3,file=FILE_SROT,status='unknown') 
		read(3,*)SROT
 	Close(3)


! 	Write(*,*) '    - Reading the degree 2 order 1 vertical displacement'
! 	open(3,file=FILE_UROT,status='unknown') 
!		read(3,*)UROT
! 	Close(3)
! 	Write(*,*) '    - Reading the degree 2 order 1 sea surface variation'
! 	open(3,file=FILE_NROT,status='unknown') 
!		read(3,*)NROT
! 	Close(3) 
!
!##################################################################  new
!##################################################################  new
!
!
! --- Reading the ALFs table from <<sh.bin>>
! 
 	Write(*,*) '    - Reading the ALFs from file sh.bin'
 	open(3,file='sh.bin',status='unknown',form='unformatted') 
		read(3)ALF
		read(3)LONG_TABLE
 	Close(3) 
!
! --- Examining the pixels table & extracting information
!
	Open(1,file='px-table.dat',status='unknown') 
	Do i=1, 4 
		Read(1,'(a12)')header
	Enddo
	Do i=1, np 
		Read (1,*) lonp(i), latp(i), anc(i), k, wet(i) 	
	Enddo
	close(1)
!
!
! --- Reading the SH OF coefficients from shof.dat 
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
! --- Reading the <<I>> array 
!
 		Write(*,*) "    - Reading array 'I'"
		open(1,file='shice.dat',status='unknown') 
		read(1,*) IIII
		close(1)
!
!
! --- Reading the E and BETA arrays 
 	Write(*,*) "    - Reading arrays 'E' and 'Beta'"
	open(1,file='ebs.dat', status='unknown')    
	open(2,file='ebu.dat',status='unknown')    
	open(3,file='ebn.dat',status='unknown')    
	open(4,file='ebv.dat',status='unknown')    
	open(7,file='ebfa.dat',status='unknown')    
	open(8,file='ebss.dat',status='unknown')    
!
	do li=0, lmax
		read(1,*) l, Es(l),  (betas(l,k),  k=0,nn+1) ! New as of Feb 2012: it was "nn" 
		read(2,*) l, Eu(l),  (betau(l,k),  k=0,nn+1) ! New as of Feb 2012: it was "nn" 
		read(3,*) l, En(l),  (betan(l,k),  k=0,nn+1) ! New as of Feb 2012: it was "nn" 
		read(4,*) l, Ev(l),  (betav(l,k),  k=0,nn+1) ! New as of Feb 2012: it was "nn" 
		read(7,*) l, Efa(l), (betafa(l,k), k=0,nn+1) ! New as of Feb 2012: it was "nn" 
		read(8,*) l, Ess(l), (betass(l,k), k=0,nn+1) ! New as of Feb 2012: it was "nn" 
	enddo
	close(8); close(7); close(4); close(3) ; close(2) ; close(1) 
!
!
! --- Computing the eustatic Z array... ! New as of Feb 2012: it was "nn"
	ze(:,:) = 0.  		
	do k=0,nn+1	
		ze(:,k) = - rhoi_o_rhow*(iiii(1,k)/oc(1))*oc(:)
	enddo
!
! --- Computing the eustatic S array... ! New as of Feb 2012: it was "nn"
	se(:,:) = (0D0,0D0)
	se(1,:) = - rhoi_o_rhow*(iiii(1,:)/oc(1)) 
!
!
! --- No ice load for imode==5 
!
      if(imode==5) IIII(:,:)=0.0 
!
!
! --- Computing the A array... ! New as of Feb 2012: it was "nn"
	aaaa(:,:)=0.
	do j=1, jmax 
 	    do k=0,NN+1
            aaaa(j,k) = ES(ll(j))*IIII(j,k)        
 	      do p=0, k
              if(p==0) aaaa(j,k) = aaaa(j,k)-(IIII(j,p)            )*BETAS(ll(j),k-p)
    	      if(p/=0) aaaa(j,k) = aaaa(j,k)-(IIII(j,p)-IIII(j,p-1))*BETAS(ll(j),k-p)  
 	      enddo 
	      aaaa(j,k)= RHOI_o_RHOE_x3*aaaa(j,k)
	      enddo
	enddo
!
!
! --- Computing the ocean average of A ! New as of Feb 2012: it was "nn" 
	do k=0, NN+1 
		aavv(k)=0.
		do j=1, jmax 
		aavv(k) = aavv(k) + & 
		          dm(j)*real(oc(j)*conjg(aaaa(j,k)))/oc(1)  
		enddo
	enddo 
!
!
! --- Computing the modified array A 
	aaaa_mod(:,:) = aaaa(:,:)
	aaaa_mod(1,:) = aaaa(1,:)-aavv(:) 
!
!
! --- Computing the H-array... ! New as of Feb 2012: it was "nn"
        x = 0.0
	hhhh = 0.0
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J,K) &
!$OMP     SHARED(X,ALF,ANC,DM,AAAA_MOD,LONG_TABLE,MM) &
!$OMP     SCHEDULE(GUIDED)
	do i=1, np    
	do k=0, NN+1
	        do j=1, jmax  
        	x(i,k) = x(i,k) + & 
			 ALF(j,anc(i))*dm(j)*real(aaaa_mod(j,k)*long_table(mm(j),i)) 
		enddo   
	    enddo
	enddo 
!$OMP END PARALLEL DO
!
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J,K) &
!$OMP       SHARED(WET,HHHH,X,ALF,ANC,LONG_TABLE,MM) &
!$OMP       SCHEDULE(GUIDED)
	do j=1, jmax 
	    do i=1, np 
	    	if(wet(i)==1) then
		    do k=0, NN+1
			hhhh(j,k) = hhhh(j,k) + x(i,k)*alf(j,anc(i))*conjg(long_table(mm(j),i))  
		    end do
		end if
	    enddo
	enddo
!$OMP END PARALLEL DO
!
! --- Updating the H-array by the eustatic Z-term 

	hhhh(:,:)=hhhh(:,:)/float(np) + ze(:,:)
!
!
! ///////////////////////////////////////////////////////////  new
! --- NEW for the implementation of the rotational feedback 
! ///////////////////////////////////////////////////////////  new 
!







! --- Computing S^{rot} at all pixels 
!     [Note: Using variable "x" in order to save memory] 
!
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J,K) &
!$OMP       SHARED(X,ALF,ANC,LONG_TABLE,SROT,MM,DM) &
!$OMP       SCHEDULE(GUIDED)
	do i=1, np    
		do k=0, NN+1   ! New as of Feb 2012: it was "nn"	
		x(i,k)=0d0 	
		do j=1, jmax 	
			x(i,k) = x(i,k) + & 
                                 ALF(j,anc(i))*dm(j)*real(SROT(j,k)*long_table(mm(j),i)) 	
		enddo
		enddo		
	enddo
!$OMP END PARALLEL DO
	write(*,*) "    - [new] S^rot at all pixels"
!
! --- Computing the harmonics of Z^{rot} by integration over wet pixels 
!
        zrot(:,:)=0d0     
!
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J,K) &
!$OMP       SHARED(WET,ZROT,X,ALF,ANC,LONG_TABLE,MM) &
!$OMP       SCHEDULE(GUIDED)
	do i=1,np
	   if(wet(i)==1) then
	       do k=0, NN+1   ! New as of Feb 2012: it was "nn"
	           do j=1, jmax 
	                 zrot(j,k) = zrot(j,k) + x(i,k)*alf(j,anc(i))*conjg(long_table(mm(j),i))  
	           enddo
	       end do
	   end if
	enddo
!$OMP END PARALLEL DO	
!
	zrot(:,:)=zrot(:,:)/float(np) 
	write(*,*) "    - [new] Harmonics of Z^rot"
!
! --- Computing the ocean average of S^{rot}  ! New as of Feb 2012: it was "nn"  
!
	do k=0, NN+1 
		srotvv(k)=0.
		do j=1, jmax 
		srotvv(k) = srotvv(k) + & 
		            dm(j)*real(oc(j)*conjg(srot(j,k)))/oc(1)  
		enddo
	enddo 
	write(*,*) "    - [new] ocean average of S^rot"
!
! --- Computing the modified array S^{rot} 
!
	srot_mod(:,:) = srot(:,:)
	srot_mod(1,:) = srot(1,:)-srotvv(:) 
	write(*,*) "    - [new] Modified S^rot array"
!	
!
! --- Computing the ocean average of Z^{rot}  ! New as of Feb 2012: it was "nn"  
!
	do k=0, NN+1 
		zrotvv(k)=0.
		do j=1, jmax 
		zrotvv(k) = zrotvv(k) + & 
		            dm(j)*real(oc(j)*conjg(zrot(j,k)))/oc(1)  
		enddo
	enddo 
	write(*,*) "    - [new] Ocean average of Z^rot"
!
! --- Computing the modified array Z^{rot} 
!
	zrot_mod(:,:) = zrot(:,:)
	zrot_mod(1,:) = zrot(1,:)-zrotvv(:) 
	write(*,*) "    - [new] Modified Z^rot array"
!
!
! --- Updating the H-array by the (modified) rotational Z-term 

	hhhh(:,:)=hhhh(:,:) + zrot_mod(:,:)
	write(*,*) "    - [new] Updating H by Z^rot modified"
!
!
! /////////////////////////////////////////////////////////// new
! --- NEW for the implementation of the rotational feedback 
! /////////////////////////////////////////////////////////// new 
!
!
!
!
!
!
!
! --- Initializing the Z and S arrays 
	Z(:,:,0) = ZE(:,:) 
	S(:,:)   = SE(:,:) 
!
!
! --- No recursion for the "Explicit approach"
	if(imode==6.or.imode==7.or.imode==3.or.SMAX==0) goto 2000     
!
!
!
! -----------------------
! ---    Recursion    ---
! -----------------------
!
	Write(*,*) "    - Starting the recursion"
!
 	do 1000 is = 1, SMAX     
!
!
        write(*,'(a12,i2,a3,i2)') '     - step ', is, ' of', SMAX
!
!
! --- Computing the 'B' array...  ! New as of Feb 2012: it was "nn"
        bbbb(:,:)=0.
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(J,K,P) &
!$OMP     SHARED(BBBB,ES,LL,Z,IS,BETAS) &
!$OMP     SCHEDULE(GUIDED)
        do j=1, jmax 
	    do k=0,NN+1
	    bbbb(j,k) = ES(ll(j))*Z(j,k,is-1)        
	        do p=0, k
	        if(p==0) bbbb(j,k) = bbbb(j,k)-(Z(j,p,is-1)-0.           )*BETAS(ll(j),k-p)
	        if(p/=0) bbbb(j,k) = bbbb(j,k)-(Z(j,p,is-1)-Z(j,p-1,is-1))*BETAS(ll(j),k-p)		 
	        enddo 
        enddo
	enddo
!$OMP END PARALLEL DO
        bbbb = RHOW_o_RHOE_X3 * bbbb
!	
!
! --- Computing the ocean-average of array B   ! New as of Feb 2012: it was "nn"
	bbvv(:)=0.
	do k=0, NN+1 
		bbvv(k)=0.
		do j=1, jmax 
		bbvv(k) = bbvv(k) + & 
		          dm(j)*real(oc(j)*conjg(bbbb(j,k)))/oc(1)
		enddo
	enddo 
!
!
! --- Computing modified 'B' array
	bbbb_mod(:,:)=bbbb(:,:)
	bbbb_mod(1,:)=bbbb(1,:)-bbvv(:) 
!
!
! --- Computing array K...  ! New as of Feb 2012: it was "nn"
	x = 0.0
	kkkk = 0.0
!
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,K,J) &
!$OMP    SHARED(X,ALF,ANC,DM,BBBB_MOD,LONG_TABLE,MM) SCHEDULE(GUIDED)
	do i=1, np 
 	   do k=0, NN+1
	      do j=1, jmax  
        	  x(i,k) = x(i,k) + & 
			 ALF(j,anc(i))*dm(j)*real(bbbb_mod(j,k)*long_table(mm(j),i)) 
              enddo   
	   enddo
	enddo 
!$OMP END PARALLEL DO
!
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J,K) &
!$OMP       SHARED(WET,KKKK,X,ALF,ANC,LONG_TABLE,MM) &
!$OMP       SCHEDULE(GUIDED)
	do j=1, jmax 
	    do i=1, np 
	       if(wet(i)==1) then
	          do k=1,NN+1
		      kkkk(j,k) = kkkk(j,k) + x(i,k)*alf(j,anc(i))*conjg(long_table(mm(j),i))  
		  enddo
	       end if
	    enddo
	enddo
!$OMP END PARALLEL DO
!	 
	kkkk(:,:)=kkkk(:,:)/float(np) 
!
!
! --- Solving for arrays 'Z' and 'S' 
	Z(:,:,is) = HHHH(:,:) + KKKK(:,:) 
!	
! /////////////////////////////////////////////////////////// new
! --- NEW for the implementation of the rotational feedback 
! /////////////////////////////////////////////////////////// new 

	S(:,:)    = AAAA_MOD(:,:) + SE(:,:) + BBBB_MOD(:,:) + SROT_MOD(:,:)
!
	SNR(:,:)  = S(:,:)     !   -  SROT_MOD(:,:)
!
! /////////////////////////////////////////////////////////// new
! --- NEW for the implementation of the rotational feedback 
! /////////////////////////////////////////////////////////// new 
!
! ------------------------------
! ---    End of recursion    ---
! ------------------------------
!
1000 CONTINUE      
!
!
!
2000 CONTINUE       
!
!
! --- Eustatic solution: U=0, S=N, V=0  
        if(imode==3.or.smax==0) then 
				U(:,:) = 0.
				S(:,:) = SE(:,:)
				N(:,:) = S(:,:)
				V(:,:) = 0. 
				goto 3000
				endif
!
! --- Array "B" for vertical displacement ! New as of Feb 2012: it was "nn"
!	bbbb(:,:)=0.
! 	do j=1, jmax 
!		do k=0,NN+1
!		bbbb(j,k) = EU(ll(j))*Z(j,k,SMAX)       
!		do p=0, k
!		if(p==0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-0.           )*BETAU(ll(j),k-p)
!		if(p/=0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-Z(j,p-1,SMAX))*BETAU(ll(j),k-p)
!		enddo 
!        	bbbb(j,k)= RHOW_o_RHOE_X3*bbbb(j,k)
!        	enddo
!	enddo
!
! --- Array "B" for vertical displacement ! New as of Feb 2012: it was "nn"
	bbbb(:,:)=0.
 	do j=1, jmax 
		do k=0,NN+1
		bbbb(j,k) = EU(ll(j))*SNR(j,k)       
		do p=0, k
		if(p==0)  bbbb(j,k) = bbbb(j,k) - (SNR(j,p)-0.      )*BETAU(ll(j),k-p)
		if(p/=0)  bbbb(j,k) = bbbb(j,k) - (SNR(j,p)-SNR(j,p-1))*BETAU(ll(j),k-p)
		enddo 
        	bbbb(j,k)= RHOW_o_RHOE_X3*bbbb(j,k)
        	enddo
	enddo
!
! --- Array "A" for vertical displacement ! New as of Feb 2012: it was "nn" 
 	aaaa(:,:)=0.
 	do j=1, jmax 
	 	do k=0,NN+1       
 	 	aaaa(j,k) = EU(ll(j))*IIII(j,k) 
 	 	do p=0, k
 	 	if(p==0) aaaa(j,k) = aaaa(j,k) - (IIII(j,p)-0.         )*BETAU(ll(j),k-p) 	 
 	 	if(p/=0) aaaa(j,k) = aaaa(j,k) - (IIII(j,p)-IIII(j,p-1))*BETAU(ll(j),k-p)
 		enddo 
 		aaaa(j,k)= RHOI_o_RHOE_X3*aaaa(j,k)
 		enddo
	enddo
!
!
! --- Vertical displacement
	U(:,:) = aaaa(:,:) + bbbb(:,:)
!
! --- Adding the Rotational effect -------- NEW NEW NEW NEW NEW 	
    U(:,:) = U(:,:) + UROT(:,:)
!
! --- Array "B" for Geoid heigth   ! New as of Feb 2012: it was "nn"
!	bbbb(:,:)=0.
! 	do j=1, jmax 
!		do k=0,NN+1
!		bbbb(j,k) = EN(ll(j))*Z(j,k,SMAX)       
!		do p=0, k
!		if(p==0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-0.           )*BETAN(ll(j),k-p)
!		if(p/=0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-Z(j,p-1,SMAX))*BETAN(ll(j),k-p)
!		enddo 
!        	bbbb(j,k)= RHOW_o_RHOE_X3*bbbb(j,k)
!        	enddo
!	enddo
!
	bbbb(:,:)=0.
 	do j=1, jmax 
		do k=0,NN+1
		bbbb(j,k) = EN(ll(j))*SNR(j,k)       
		do p=0, k
		if(p==0)  bbbb(j,k) = bbbb(j,k) - (SNR(j,p)-0.        )*BETAN(ll(j),k-p)
		if(p/=0)  bbbb(j,k) = bbbb(j,k) - (SNR(j,p)-SNR(j,p-1))*BETAN(ll(j),k-p)
		enddo 
        	bbbb(j,k)= RHOW_o_RHOE_X3*bbbb(j,k)
        	enddo
	enddo
!
! --- Array "A" for Geoid heigth   ! New as of Feb 2012: it was "nn" 
 	aaaa(:,:)=0.
 	do j=1, jmax 
	 	do k=0,NN+1     
 	 	aaaa(j,k) = EN(ll(j))*IIII(j,k) 
 	 	do p=0, k
 	 	if(p==0) aaaa(j,k) = aaaa(j,k) - (IIII(j,p)-0.         )*BETAN(ll(j),k-p) 	 
 	 	if(p/=0) aaaa(j,k) = aaaa(j,k) - (IIII(j,p)-IIII(j,p-1))*BETAN(ll(j),k-p)
 		enddo 
 		aaaa(j,k)= RHOI_o_RHOE_X3*aaaa(j,k)
 		enddo
	enddo
!
! --- Geoid undulations 
	N(:,:) = aaaa(:,:) + bbbb(:,:)
!
! --- Adding the Rotational effect -------- NEW NEW NEW NEW NEW 	
    N(:,:) = N(:,:) + NROT(:,:)
!
! --- Normalized gravity potential is "N" without the additional "c" constant... 
!     New as March 2011
!
	G(:,:) = (0d0,0d0)
        G(:,:) = N(:,:)
!
! ---
! --- <<Water height equivalent>>, new as of August 2012 (Chambers et al. 2010 paper) ========= in progress
! ---     
	open(112,file='k.dat',status='unknown') 
	Do i=1, 2 
		Read(112,'(a12)')header
	Enddo 	
        do L=1, LMAX 
		Read(112,*) IJ, KPRIME(L)
		KPRIME(L)=1D0+KPRIME(L)
	enddo
	close(112)
!	
        W(:,:) = (0d0,0d0)
!
        do j=4, jmax 
		W(j,:) = (G(j,:)/RHOW_o_RHOE_X3)*(2D0*FLOAT(LJ(J))+1D0)/KPRIME(LJ(J))
        enddo
! ---
! --- <<Water height equivalent>>, new as of August 2012 (Chambers et al. 2010 paper) ========= in progress
! ---        
!
! --- Adding a constant to geoid undulations 
	N(1,:) = N(1,:) +  SE(1,:) - AAVV(:) - BBVV(:)  
!
! --- Geoid undulations (previous formulation) 
!       N(:,:) = S(:,:) + U(:,:)
!
!
! --- Array "B" for horizontal displacement ! New as of Feb 2012: it was "nn" 
	bbbb(:,:)=0.
 	do j=1, jmax 
		do k=0,NN+1
		bbbb(j,k) = EV(ll(j))*Z(j,k,SMAX)       
		do p=0, k
		if(p==0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-0.           )*BETAV(ll(j),k-p)
		if(p/=0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-Z(j,p-1,SMAX))*BETAV(ll(j),k-p)
		enddo 
        	bbbb(j,k)= RHOW_o_RHOE_X3*bbbb(j,k)
        	enddo
	enddo
!
! --- Array "A" for horizontal displacement ! New as of Feb 2012: it was "nn" 
 	aaaa(:,:)=0.
 	do j=1, jmax 
	 	do k=0,NN+1       
 	 	aaaa(j,k) = EV(ll(j))*IIII(j,k) 
 	 	do p=0, k
 	 	if(p==0) aaaa(j,k) = aaaa(j,k) - (IIII(j,p)-0.         )*BETAV(ll(j),k-p) 	 
 	 	if(p/=0) aaaa(j,k) = aaaa(j,k) - (IIII(j,p)-IIII(j,p-1))*BETAV(ll(j),k-p)
 		enddo 
 		aaaa(j,k)= RHOI_o_RHOE_X3*aaaa(j,k)
 		enddo
	enddo
!
! --- Horizontal displacement
	V(:,:) = aaaa(:,:) + bbbb(:,:)
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   FA Gravity DRAFT 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   FA Gravity DRAFT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   FA Gravity DRAFT
!
! --- Array "B" for FA Gravity ! New as of Feb 2012: it was "nn" 
	bbbb(:,:)=0.
 	do j=1, jmax 
		do k=0,NN+1
		bbbb(j,k) = EFA(ll(j))*Z(j,k,SMAX)       
		do p=0, k
		if(p==0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-0.           )*BETAFA(ll(j),k-p)
		if(p/=0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-Z(j,p-1,SMAX))*BETAFA(ll(j),k-p)
		enddo 
        	bbbb(j,k)= RHOW_o_RHOE_X3*100.*1000.*(GRA_REF/RAD_REF)*bbbb(j,k)      
        	enddo
	enddo
!
! --- Array "A" for FA Gravity ! New as of Feb 2012: it was "nn" 
 	aaaa(:,:)=0.
 	do j=1, jmax 
	 	do k=0,NN+1       
 	 	aaaa(j,k) = EFA(ll(j))*IIII(j,k) 
 	 	do p=0, k
 	 	if(p==0) aaaa(j,k) = aaaa(j,k) - (IIII(j,p)-0.         )*BETAFA(ll(j),k-p) 	 
 	 	if(p/=0) aaaa(j,k) = aaaa(j,k) - (IIII(j,p)-IIII(j,p-1))*BETAFA(ll(j),k-p)
 		enddo 
 		aaaa(j,k)= RHOI_o_RHOE_X3*100.*1000.*(GRA_REF/RAD_REF)*aaaa(j,k)
 		enddo
	enddo
!
! --- FA Gravity gravity anomaly                      (in units of "MILLIGAL")
	FA(:,:) = aaaa(:,:) + bbbb(:,:)
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   FA Gravity DRAFT 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   FA Gravity DRAFT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   FA Gravity DRAFT
!
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   SS Gravity DRAFT 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   SS Gravity DRAFT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   SS Gravity DRAFT
!
! --- Array "B" for SS Gravity ! New as of Feb 2012: it was "nn" 
	bbbb(:,:)=0.
 	do j=1, jmax 
		do k=0,NN+1
		bbbb(j,k) = ESS(ll(j))*Z(j,k,SMAX)       
		do p=0, k
		if(p==0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-0.           )*BETASS(ll(j),k-p)
		if(p/=0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-Z(j,p-1,SMAX))*BETASS(ll(j),k-p)
		enddo 
        	bbbb(j,k)= RHOW_o_RHOE_X3*100.*1000.*(GRA_REF/RAD_REF)*bbbb(j,k)      
        	enddo
	enddo
!
! --- Array "A" for SS Gravity ! New as of Feb 2012: it was "nn" 
 	aaaa(:,:)=0.
 	do j=1, jmax 
	 	do k=0,NN+1       
 	 	aaaa(j,k) = ESS(ll(j))*IIII(j,k) 
 	 	do p=0, k
 	 	if(p==0) aaaa(j,k) = aaaa(j,k) - (IIII(j,p)-0.         )*BETASS(ll(j),k-p) 	 
 	 	if(p/=0) aaaa(j,k) = aaaa(j,k) - (IIII(j,p)-IIII(j,p-1))*BETASS(ll(j),k-p)
 		enddo 
 		aaaa(j,k)= RHOI_o_RHOE_X3*100.*1000.*(GRA_REF/RAD_REF)*aaaa(j,k)
 		enddo
	enddo
!
! --- SS Gravity gravity anomaly          (in units of "MILLIGAL")
	SS(:,:) = aaaa(:,:) + bbbb(:,:)
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   SS Gravity DRAFT 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   SS Gravity DRAFT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   SS Gravity DRAFT
!
!
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Load components DRAFT  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Load components DRAFT  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Load components DRAFT  
!
       LO_I(:,:)  = (0d0,0d0)
       LO_O(:,:)  = (0d0,0d0)
       LO_T(:,:)  = (0d0,0d0)
!
! --- Ice  component of load    ======== RIGID EARTH 
	LO_I(:,:) = RHOI*IIII(:,:)
!
!
! --- Ocean component of load   ======== RIGID EARTH
	LO_O(:,:) = RHOW*ZE(:,:)
!
!
! --- Total load (ice + ocean)  ======== RIGID EARTH
	LO_T(:,:) = LO_I(:,:) + LO_O(:,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Load components DRAFT  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Load components DRAFT  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Load components DRAFT  
!
!
3000 CONTINUE 
!
!
!
!
!
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    - Reporting the output arrays -
!           (SH coefficients)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
! --- Eustatic Sea level change ---------------- New as of April 2012 
	Out_Filename=filee
 	open(3,file=Out_Filename,status='unknown',form='unformatted') 
        write(3) SE 
	close(3)
! --- Eustatic Sea level change ---------------- New as of April 2012 
!
!
! --- Sea level change  
	Out_Filename=files 
 	open(3,file=Out_Filename,status='unknown',form='unformatted') 
        write(3) S 
	close(3)
!
! --- Vertical displacement 
	Out_Filename=fileu 
 	open(3,file=Out_Filename,status='unknown',form='unformatted') 
        write(3) U 
	close(3) 
!
! --- Geoid undulations   
	Out_Filename=filen 
 	open(3,file=Out_Filename,status='unknown',form='unformatted') 
        write(3) N 
	close(3) 
!
! --- Normalized gravity potential (new as of March 2011) ---------- 
	Out_Filename=fileg 
 	open(3,file=Out_Filename,status='unknown',form='unformatted') 
        write(3) G 
	close(3) 
!
! --- Equivalent Water height (new as of August 4, 2012) -----------  In progress.... 
	Out_Filename=filew 
 	open(3,file=Out_Filename,status='unknown',form='unformatted') 
        write(3) W 
	close(3) 
!
! --- Reduced sea level change 
	Out_Filename=filez 
 	open(3,file=Out_Filename,status='unknown',form='unformatted') 
        write(3) Z 
	close(3)
!
! --- Horizontal displacement 
	Out_Filename=filev 
 	open(3,file=Out_Filename,status='unknown',form='unformatted') 
        write(3) V 
	close(3)
!
! --- Free Air gravity anomalies  
	Out_Filename=filefa
 	open(3,file=Out_Filename,status='unknown',form='unformatted') 
        write(3) FA 
	close(3)
!
! --- Solid Surface gravity anomalies 
	Out_Filename=filess
 	open(3,file=Out_Filename,status='unknown',form='unformatted') 
        write(3) SS 
	close(3)
!
! --- Ice component of load                                              <<<< new 
!Out_Filename=fileli 
!open(3,file=Out_Filename,status='unknown',form='unformatted') 
!write(3) LO_I 
!close(3)
!
! --- Ocean component of load          				         <<<< new 
!Out_Filename=filelo 
!open(3,file=Out_Filename,status='unknown',form='unformatted') 
!write(3) LO_O 
!close(3)
!
! --- Total (ice+ocean) load             			         <<<< new 
!Out_Filename=filelt 
!open(3,file=Out_Filename,status='unknown',form='unformatted') 
!write(3) LO_T 
!close(3)
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Write(*,*) "    +++ SLE solved <3 +++" 
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
! --- Free up memory space
!
 DEALLOCATE( LL, MM, DM, ANC, WET )
 DEALLOCATE( ALF, LONP, LATP, X )
 DEALLOCATE( AAVV, BBVV ) 
 DEALLOCATE( SROTVV, ZROTVV )
 DEALLOCATE( LONG_TABLE, OC )
 DEALLOCATE( IIII, ZE, SE )
 DEALLOCATE( AAAA, AAAA_MOD )
 DEALLOCATE( BBBB, BBBB_MOD )
 DEALLOCATE( SROT, SROT_MOD )
 DEALLOCATE( ZROT, ZROT_MOD )
 DEALLOCATE( UROT, NROT ) 
 DEALLOCATE( HHHH, KKKK )
 DEALLOCATE( BETAS , ES  )
 DEALLOCATE( BETAN , EN  )
 DEALLOCATE( BETAU , EU  )
 DEALLOCATE( BETAV , EV  )     
 DEALLOCATE( BETAFA, EFA )     
 DEALLOCATE( BETASS, ESS )     
 DEALLOCATE( S, N, U, G )
 DEALLOCATE( V, FA, SS )
 DEALLOCATE( LO_I, LO_O, LO_T )
 DEALLOCATE( Z )
 DEALLOCATE( SNR, W, KPRIME )
!
!
!
	End program SLE
!
!
!
!
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
 SUBROUTONE S_ROTAZ(SROT)
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

















