!!
! This is program "SLE_VAROC.F90" 
!    ... 
!    GS & FC June 22 2008 "Variable coastlines"
!    Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
!    *** Reviewed GS & FC November 2009 - Porting under gfortran 
!    Also revised by GS April 21, 2010 - degree 1 implementation by ALMA
!    Also revised by GS April 25, 2010 - Cleaning and new geoid implememtation
!    Revised May 26 GS for new ice SH  
!    Updated GS July 2010 for the imode=5 option (see Gabriella) 
!    Updated GS June 2010 for Free Air & SS Gravity implementation 
! ***  Revised GS August 2010 - g95 - Double precision implementation
!    Revised DM November 2011 - Loop parallelization
!    Feb 2012: Implementation of the numerical derivative "in the future" 
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
!       - px-table-XX.dat 
!	- shof-XX.dat 
!	- shice.dat
!	- ebs.dat
!	- ebu.dat
!	- ebn.dat
!	- ebv.dat
! 	- sh.bin
!       - anchor.tmp   [NEW]
!
! Output files:
!	- shs.bin 
!	- shu.bin 
!	- shn.bin 
!	- shv.bin 
!       - shz.bin 
!       - shfa.bin <<< new... 
!       - shss.bin <<< new... 
!
!
!INCLUDE "harmonics.f90"
 PROGRAM SLE_VAROC
 IMPLICIT NONE 
 INCLUDE "data.inc" 
!
 CHARACTER*12 HEADER 
 CHARACTER*2 LABCHAR(0:NN+1)
 INTEGER :: NANCH
 INTEGER I, J, K, L, P, IS, LJ, MJ, LI, KT, DOM, IND 
 INTEGER, ALLOCATABLE :: LL(:), MM(:), DM(:), ANC(:), WET(:,:)   
!
 REAL*4 RHOE, RHOI_O_RHOE_X3, RHOW_O_RHOE_X3, RHOI_O_RHOW 
!
 REAL*8, ALLOCATABLE :: ALF(:,:), X(:,:)    
 REAL*8 JUNK, RESH, IMSH
 REAL*8, ALLOCATABLE :: AAVV(:), BBVV(:)  
!
 COMPLEX*16, ALLOCATABLE :: LONG_TABLE(:,:), OC(:,:)	   
 COMPLEX*16, ALLOCATABLE :: IIII(:,:), ZE(:,:), SE(:,:)    
 COMPLEX*16, ALLOCATABLE :: AAAA(:,:), AAAA_MOD(:,:)	   
 COMPLEX*16, ALLOCATABLE :: BBBB(:,:), BBBB_MOD(:,:)	   
 COMPLEX*16, ALLOCATABLE :: HHHH(:,:), KKKK(:,:)	   
! 
!--------------------------------------------------
! Convolved Green functions (viscous & elastic) 
! 
 REAL*4, ALLOCATABLE :: BETAS(:,:),  ES(:), &  ! Sea level change  
                     BETAN(:,:),  EN(:),    &  ! Geoid undulations 
                     BETAU(:,:),  EU(:),    &  ! Vertical displacemnt
                     BETAV(:,:),  EV(:),    &  ! Horizontal displacement  
                     BETAFA(:,:), EFA(:),   &  ! Free air (FA) gravity anomaly 
                     BETASS(:,:), ESS(:)       ! Solid-surface (SS) gravity anomaly      !
!--------------------------------------------------
! Output arrays & filenames 
! 
 CHARACTER*22 Out_Filename  
!
! - Sea level change ------------------------------
 complex*16, allocatable :: S(:,:)   
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
 complex*16, allocatable ::  Z(:,:,:)
 CHARACTER*22, PARAMETER :: FILEZ='shz.bin'  
!--------------------------------------------------
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
! --- Allocate memory space
!
 ALLOCATE( LL(JMAX), MM(JMAX), DM(JMAX), ANC(NP), WET(NP,0:NN+1) )
 ALLOCATE( ALF(JMAX,NANCH), X(NP,0:NN+1) )
 ALLOCATE( AAVV(0:NN+1), BBVV(0:NN+1) )  
 ALLOCATE( LONG_TABLE(0:LMAX,NP), OC(JMAX,0:NN+1) )
 ALLOCATE( IIII(JMAX,0:NN+1), ZE(JMAX,0:NN+1), SE(JMAX,0:NN+1) )
 ALLOCATE( AAAA(JMAX,0:NN+1), AAAA_MOD(JMAX,0:NN+1) )
 ALLOCATE( BBBB(JMAX,0:NN+1), BBBB_MOD(JMAX,0:NN+1) )
 ALLOCATE( HHHH(JMAX,0:NN+1), KKKK(JMAX,0:NN+1) )
 ALLOCATE( BETAS (0:LMAX,0:NN+1), ES (0:LMAX) )
 ALLOCATE( BETAN (0:LMAX,0:NN+1), EN (0:LMAX) )
 ALLOCATE( BETAU (0:LMAX,0:NN+1), EU (0:LMAX) )
 ALLOCATE( BETAV (0:LMAX,0:NN+1), EV (0:LMAX) )     
 ALLOCATE( BETAFA(0:LMAX,0:NN+1), EFA(0:LMAX) )     
 ALLOCATE( BETASS(0:LMAX,0:NN+1), ESS(0:LMAX) )     
 ALLOCATE( S(JMAX,0:NN+1), N(JMAX,0:NN+1), U(JMAX,0:NN+1), G(JMAX,0:NN+1) )
 ALLOCATE( V(JMAX,0:NN+1), FA(JMAX,0:NN+1), SS(JMAX,0:NN+1) )
 ALLOCATE( LO_I(JMAX,0:NN+1), LO_O(JMAX,0:NN+1), LO_T(JMAX,0:NN+1) )
 ALLOCATE( Z(JMAX,0:NN+1,0:SMAX) )
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
! --- A time-dependent character label ! Updated february 2012 
!
	Do kt=0, nn+1 	
		write( labchar(kt), '(i2.2)' ) kt
		!open(3,file='junk.dat',status='unknown') 
		!if(kt<=9) write(3,'(a1,i1)') '0',kt  
		!if(kt> 9) write(3,'(i2)')        kt
		!close(3)
        	!open(3,file='junk.dat',status='unknown')
		!read(3,'(a2)') labchar(kt)
		!close(3)
	Enddo
!
!
! --- Pre-computing 'l' and 'm' corresponding to degree 'J'
	do j=1, jmax 
		mm(j)=mj(j) 
		ll(j)=lj(j)
		dm(j)=2-dom(j)
	enddo	
!
! --- Reading the ALFs table from <<sh.bin>>
 	open(3,file='sh.bin',status='unknown',form='unformatted') 
		read(3)ALF
		read(3)LONG_TABLE
 	Close(3) 
!
! --- Reading the tables with ocean function from previous step ! Updated february 2012 	
        Do kt =0, nn+1         
 	if(kt==0.or.kt==nn+1)Write(*,*) "    - Time ", kt, " of", nn+1
		Open(1,file='px-table-'//labchar(kt)//'.dat',status='unknown') 
		Do i=1, 4 ; Read(1,'(a12)')header ; Enddo
			Do i=1, np 
			Read (1,*) JUNK, JUNK, anc(i), k, wet(i,kt) 	
			Enddo
		close(1)
	Enddo
!
! --- Reading the SH OF coefficients from shof-XX.dat ! Updated february 2012 
 	Write(*,*) "    - sle.f: Reading the SH OF "
        Do kt =0, nn +1
 	if(kt==0.or.kt==nn+1)Write(*,*) "    - Time ", kt, " of", nn
 		open(3,file='shof-'//labchar(kt)//'.dat',status='unknown')
 		do j=1, jmax   
			read(3,*) k, resh, imsh 
                	oc(j,kt)=cmplx(resh, imsh)	
 		enddo
 		close(3)
	Enddo
!
! --- Reading the <<I>> array 
!IIII(:,:)=0. 
!If(imode/=5) then 
		open(1,file='shice.dat',status='unknown') 
		read(1,*) IIII
		close(1)
!Endif
!
!
! --- Reading the E and BETA arrays 
	open(1,file= 'ebs.dat', status='unknown')    
	open(2,file= 'ebu.dat',status='unknown')    
	open(3,file= 'ebn.dat',status='unknown')    
	open(4,file= 'ebv.dat',status='unknown')  
	open(7,file='ebfa.dat',status='unknown')    
	open(8,file='ebss.dat',status='unknown')      
!
	do li=0, lmax
		read(1,*) l, Es(l), (betas(l,k), k=0,nn+1)  ! Updated february 2012 
		read(2,*) l, Eu(l), (betau(l,k), k=0,nn+1)  ! Updated february 2012 
		read(3,*) l, En(l), (betan(l,k), k=0,nn+1)  ! Updated february 2012 
		read(4,*) l, Ev(l), (betav(l,k), k=0,nn+1)  ! Updated february 2012 
		read(7,*) l, Efa(l), (betafa(l,k), k=0,nn+1)  ! Updated february 2012 
		read(8,*) l, Ess(l), (betass(l,k), k=0,nn+1)  ! Updated february 2012 
	enddo
	close(8); close(7); close(4); close(3) ; close(2) ; close(1) 
!
!
! --- Computing the eustatic Z array...  ! Updated february 2012
	ze(:,:) = 0.  		
	do k=0,nn+1	
		ze(:,k) = - rhoi_o_rhow*(iiii(1,k)/oc(1,k))*oc(:,k)
	enddo
!
!
! --- Computing the eustatic S array...  ! Updated february 2012
 	se(:,:) = 0.
	se(1,:) = - rhoi_o_rhow*(iiii(1,:)/oc(1,:)) 
!	
!
! --- No ice load for imode==5 
!
      if(imode==5) IIII(:,:)=0.0 
!
!
! --- Computing the A array...  ! Updated february 2012
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
! --- Computing the ocean average of A  ! Updated february 2012 
	do k=0, NN+1 
		aavv(k)=0.
		do j=1, jmax 
		aavv(k) = aavv(k) + & 
		          dm(j)*real(oc(j,k)*conjg(aaaa(j,k)))/oc(1,k)  
		enddo
	enddo 
!
!
! --- Computing the modified ocean average of A 
	aaaa_mod(:,:) = aaaa(:,:)
	aaaa_mod(1,:) = aaaa(1,:)-aavv(:) 
!
!
! --- Computing the R-array...  ! Updated february 2012
        x = 0.
        hhhh = 0.	
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
	    	do k=0, NN+1  
	            if(wet(i,k)==1) hhhh(j,k) = hhhh(j,k) + x(i,k)*alf(j,anc(i))*conjg(long_table(mm(j),i))  
	        enddo
	    enddo
	enddo
!$OMP END PARALLEL DO
!	 
	hhhh(:,:)=hhhh(:,:)/float(np) + ze(:,:)
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
! --- Computing the 'B' array...  ! Updated february 2012
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
!
        bbbb = RHOW_o_RHOE_X3 * bbbb

!	
!
! --- Computing the ocean-average of array B   ! Updated february 2012
	bbvv(:)=0.
	do k=0, NN+1 
		bbvv(k)=0.
		do j=1, jmax 
		bbvv(k) = bbvv(k) + & 
		          dm(j)*real(oc(j,k)*conjg(bbbb(j,k)))/oc(1,k)
		enddo
	enddo 
!
!
! --- Computing modified 'B' array
	bbbb_mod(:,:)=bbbb(:,:)
	bbbb_mod(1,:)=bbbb(1,:)-bbvv(:) 
!
!
! --- Computing array K...  ! Updated february 2012
!
        x = 0.
	kkkk = 0.
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
		do k=0, NN+1
	    if(wet(i,k)==1) kkkk(j,k) = & 
			  kkkk(j,k) + x(i,k)*alf(j,anc(i))*conjg(long_table(mm(j),i))  
	    enddo
	enddo
	enddo 
!$OMP END PARALLEL DO
!
	kkkk(:,:)=kkkk(:,:)/float(np) 
!
!
! --- Solving for arrays 'Z' and 'S' 
	Z(:,:,is) = HHHH(:,:) + KKKK(:,:) 
	S(:,:)    = AAAA_MOD(:,:) + SE(:,:) + BBBB_MOD(:,:)
!
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
! --- Eustatic solution: U=0, S=N.  
        if(imode==3.or.smax==0) then 
				U(:,:) = 0.
				S(:,:) = SE(:,:)
				N(:,:) = S(:,:)
 				V(:,:) = 0. 
				goto 3000
				endif
!
! --- Array "B" for vertical displacement  ! Updated february 2012 
 	bbbb(:,:)=0.
 	do j=1, jmax 
		do k=0,NN+1
		bbbb(j,k) = EU(ll(j))*Z(j,k,SMAX)       
		do p=0, k
		if(p==0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-0.           )*BETAU(ll(j),k-p)
		if(p/=0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-Z(j,p-1,SMAX))*BETAU(ll(j),k-p)
		enddo 
        	bbbb(j,k)= RHOW_o_RHOE_X3*bbbb(j,k)
        	enddo
	enddo
!
! --- Array "A" for vertical displacement  ! Updated february 2012 
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
!
! --- Array "B" for Geoid heigth  ! Updated february 2012   
	bbbb(:,:)=0.
 	do j=1, jmax 
		do k=0,NN+1
		bbbb(j,k) = EN(ll(j))*Z(j,k,SMAX)       
		do p=0, k
		if(p==0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-0.           )*BETAN(ll(j),k-p)
		if(p/=0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-Z(j,p-1,SMAX))*BETAN(ll(j),k-p)
		enddo 
        	bbbb(j,k)= RHOW_o_RHOE_X3*bbbb(j,k)
        	enddo
	enddo
!
! --- Array "A" for Geoid heigth   ! Updated february 2012 
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
! --- Normalized gravity potential is "N" without the additional "c" constant... 
!     New as March 2011
!
	G(:,:) = (0d0,0d0)
        G(:,:) = N(:,:)
!
! --- Adding a constant to geoid undulations 
	N(1,:) = N(1,:) +  SE(1,:) - AAVV(:) - BBVV(:)  
!
! --- Geoid undulations (previous formulation) 
!       N(:,:) = S(:,:) + U(:,:)
!
!
! --- Array "B" for horizontal displacement  ! Updated february 2012 
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
! --- Array "A" for horizontal displacement  ! Updated february 2012 
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
!
! --- Horizontal displacement
	V(:,:) = aaaa(:,:) + bbbb(:,:)
!
!
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   FA Gravity DRAFT 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   FA Gravity DRAFT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   FA Gravity DRAFT
!
! --- Array "B" for FA Gravity   ! Updated february 2012
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
! --- Array "A" for FA Gravity  ! Updated february 2012 
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
! --- Array "B" for SS Gravity  ! Updated february 2012 
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
! --- Array "A" for SS Gravity  ! Updated february 2012 
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
! --- SS Gravity gravity anomaly                      (in units of "MILLIGAL")
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
! --- Ice  component of load 
	LO_I(:,:)=RHOI*IIII(:,:)
!
!
! --- Ocean component of load 
	LO_O(:,:)=RHOW*Z(:,:,SMAX)
!
!
! --- Total load (ice + ocean)  
	LO_T(:,:)= LO_I(:,:) + LO_O(:,:)
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
!
!
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    - Reporting the output arrays -
!           (SH coefficients)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
	Out_Filename=fileli 
 	open(3,file=Out_Filename,status='unknown',form='unformatted') 
        write(3) LO_I 
	close(3)
!
! --- Ocean component of load          				         <<<< new 
	Out_Filename=filelo 
 	open(3,file=Out_Filename,status='unknown',form='unformatted') 
        write(3) LO_O 
	close(3)
!
! --- Total (ice+ocean) load             				 <<<< new 
	Out_Filename=filelt 
 	open(3,file=Out_Filename,status='unknown',form='unformatted') 
        write(3) LO_T 
	close(3)
!
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
 DEALLOCATE( ALF, X )
 DEALLOCATE( AAVV, BBVV )  
 DEALLOCATE( LONG_TABLE, OC )
 DEALLOCATE( IIII, ZE, SE )
 DEALLOCATE( AAAA, AAAA_MOD )
 DEALLOCATE( BBBB, BBBB_MOD )
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
!
	End program SLE_VAROC
!
!
!
!
!
!
