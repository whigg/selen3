!
! This is program "SLE.F90" 
!
!               GS 04-11-2008 "Intel port..."
! 		Modified by GS september 2008 for horizontal movements 
! 		Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! 		*** Reviewed GS & FC November 2009 - Porting under gfortran 
! 		Also reviewed in April 2010 for the implementation of degree 1 
! 		(remember the mode coupling issue...) 
! 		Revised May 26 GS for new ice SH  
!               Updated GS June 2010 for Free Air & SS Gravity implementation 
!               Updated GS July 2010 for the imode=5 option (see Gabriella) 
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
!	- eb.dat
!	- ebu.dat
!       - ebv.dat 
!       - ebfa.dat  <<< new... 
! 	- sh.bin
!
! Output files:
!	- shu.bin 
!	- shn.bin 
!	- shs.bin 
!       - shz.bin 
!       - shv.bin 
!       - shfa.bin    <<< new... 
!       - shss.bin    <<< new... 
!       - shload*.bin <<< new 
!
 INCLUDE "harmonics.f90"
 PROGRAM SLE 
 IMPLICIT NONE 
 INCLUDE "data.inc"			    
!
 CHARACTER*12 HEADER 
 INTEGER, PARAMETER :: NANCH=FLOAT(NP)/7. 
 INTEGER I, J, K, L, P, IS, LJ, MJ, LI, DOM, IND
 INTEGER LL(JMAX), MM(JMAX), DM(JMAX), ANC(NP), WET(NP) 
!
 REAL*4 RHOE, RHOI_O_RHOE_X3, RHOW_O_RHOE_X3, RHOI_O_RHOW 
 REAL*4 ALF(JMAX,NANCH), LONP(NP), LATP(NP), X(NP,0:NN)  
 REAL*4 RESH, IMSH, AAVV(0:NN), BBVV(0:NN)
!
 COMPLEX*8 LONG_TABLE(0:LMAX,NP), OC(JMAX) 
 COMPLEX*8 IIII(JMAX,0:NN), ZE(JMAX,0:NN), SE(JMAX,0:NN)
 COMPLEX*8 AAAA(JMAX,0:NN), AAAA_MOD(JMAX,0:NN)
 COMPLEX*8 BBBB(JMAX,0:NN), BBBB_MOD(JMAX,0:NN) 
 COMPLEX*8 HHHH(JMAX,0:NN), KKKK(JMAX,0:NN) 
! 
!--------------------------------------------------
! Convolved Green functions (viscous & elastic) 
! 
 REAL*4 BETAS(0:LMAX,0:NN),  ES(0:LMAX),   &   ! Sea level change  
        BETAN(0:LMAX,0:NN),  EN(0:LMAX),   &   ! Geoid undulations 
        BETAU(0:LMAX,0:NN),  EU(0:LMAX),   &   ! Vertical displacemnt
        BETAV(0:LMAX,0:NN),  EV(0:LMAX),   &   ! Horizontal displacement  
        BETAFA(0:LMAX,0:NN), EFA(0:LMAX),  &   ! Free air (FA) gravity anomaly 
        BETASS(0:LMAX,0:NN), ESS(0:LMAX)       ! Solid-surface (SS) gravity anomaly                  
!
!--------------------------------------------------
! Output arrays & filenames 
! 
 CHARACTER*22 Out_Filename  
!
! - Sea level change ------------------------------
 COMPLEX*8 S(JMAX,0:NN)   
 CHARACTER*22, PARAMETER :: FILES='shs.bin'   
!
! - Geoid undulations -----------------------------
 COMPLEX*8 N(JMAX,0:NN)   
 CHARACTER*22, PARAMETER :: FILEN='shn.bin'   
!
! - Vertical displacemnt --------------------------
 COMPLEX*8 U(JMAX,0:NN)   
 CHARACTER*22, PARAMETER :: FILEU='shu.bin'   
!
! - Horizontal displacement -----------------------
 COMPLEX*8 V(JMAX,0:NN)   
 CHARACTER*22, PARAMETER :: FILEV='shv.bin'  
!
! - Free air (FA) gravity anomaly -----------------
 COMPLEX*8 FA(JMAX,0:NN)   
 CHARACTER*22, PARAMETER :: FILEFA='shfa.bin'  
!
! - Solid-surface (SS) gravity anomaly ------------
 COMPLEX*8 SS(JMAX,0:NN)   
 CHARACTER*22, PARAMETER :: FILESS='shss.bin' 
!
! - Ice component of load  ------------------------
 COMPLEX*8 LO_I(JMAX,0:NN)
 CHARACTER*22, PARAMETER :: FILELI='shload_ice.bin'
!
! - Ocean component of load -----------------------  
 COMPLEX*8 LO_O(JMAX,0:NN)
 CHARACTER*22, PARAMETER :: FILELO='shload_oce.bin' 
!
! - Total load  -----------------------------------
 COMPLEX*8 LO_T(JMAX,0:NN)
 CHARACTER*22, PARAMETER :: FILELT='shload.bin'  
!
! - Sea level times the ocean function ------------  
 COMPLEX*8  Z(JMAX,0:NN,0:SMAX)
 CHARACTER*22, PARAMETER :: FILEZ='shz.bin'  
!--------------------------------------------------

!
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
 CALL AVERAGE_EARTH_DENSITY(IND, RHOE)
 RHOI_O_RHOE_X3 = 3.*RHOI/RHOE 
 RHOW_O_RHOE_X3 = 3.*RHOW/RHOE 
 RHOI_O_RHOW    =    RHOI/RHOW   
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
!IIII(:,:)=0. 
!If(imode/=5) then 
 		Write(*,*) "    - Reading array 'I'"
		open(1,file='shice.dat',status='unknown') 
		read(1,*) IIII
		close(1)
!Endif
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
	open(9,file='ebt.dat' .....)    
!
	do li=0, lmax
		read(1,*) l, Es(l),  (betas(l,k),  k=0,nn) 
		read(2,*) l, Eu(l),  (betau(l,k),  k=0,nn) 
		read(3,*) l, En(l),  (betan(l,k),  k=0,nn) 
		read(4,*) l, Ev(l),  (betav(l,k),  k=0,nn) 
		read(7,*) l, Efa(l), (betafa(l,k), k=0,nn) 
		read(8,*) l, Ess(l), (betass(l,k), k=0,nn) 
		..... 
	enddo
	close(8); close(7); close(4); close(3) ; close(2) ; close(1) 
!
!
! --- Computing the eustatic Z array...
	ze(:,:) = 0.  		
	do k=0,nn	
		ze(:,k) = - rhoi_o_rhow*(iiii(1,k)/oc(1))*oc(:)
	enddo
!
!
! --- Computing the eustatic S array...
	se(:,:) = 0.
	se(1,:) = - rhoi_o_rhow*(iiii(1,:)/oc(1)) 
!
!
! --- No ice load for imode==5 
!
      if(imode==5) IIII(:,:)=0.0 
!
!
! --- Computing the A array...
	aaaa(:,:)=0.
	do j=1, jmax 
 	    do k=0,NN
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
! --- Computing the ocean average of A 
	do k=0, NN 
		aavv(k)=0.
		do j=1, jmax 
		aavv(k) = aavv(k) + & 
		          dm(j)*real(oc(j)*conjg(aaaa(j,k)))/oc(1)  
		enddo
	enddo 
!
!
! --- Computing the modified ocean average of A 
	aaaa_mod(:,:) = aaaa(:,:)
	aaaa_mod(1,:) = aaaa(1,:)-aavv(:) 
!
!
! --- Computing the R-array...
	do i=1, np 
!
!		   
	do k=0, NN
	    x(i,k)=0.			
	        do j=1, jmax  
        	x(i,k) = x(i,k) + & 
			 ALF(j,anc(i))*dm(j)*real(aaaa_mod(j,k)*long_table(mm(j),i)) 
		enddo   
	    enddo
	enddo 
!
	do k=0, NN
	do j=1, jmax 
	    hhhh(j,k)=0.  
	    do i=1, np 
	    if(wet(i)==1) hhhh(j,k) = & 
			  hhhh(j,k) + x(i,k)*alf(j,anc(i))*conjg(long_table(mm(j),i))  
	    enddo
	enddo
	enddo 
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
! --- Computing the 'B' array...
        bbbb(:,:)=0.
        do j=1, jmax 
	    do k=0,NN
	    bbbb(j,k) = ES(ll(j))*Z(j,k,is-1)        
	        do p=0, k
	        if(p==0) bbbb(j,k) = bbbb(j,k)-(Z(j,p,is-1)-0.           )*BETAS(ll(j),k-p)
	        if(p/=0) bbbb(j,k) = bbbb(j,k)-(Z(j,p,is-1)-Z(j,p-1,is-1))*BETAS(ll(j),k-p)		 
	        enddo 
            bbbb(j,k)= RHOW_o_RHOE_X3*bbbb(j,k)
        enddo
	enddo
!	
!
! --- Computing the ocean-average of array B array
	bbvv(:)=0.
	do k=0, NN 
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
! --- Computing array K...
	do i=1, np 
!	
	do k=0, NN
	    x(i,k)=0.			
	        do j=1, jmax  
        	x(i,k) = x(i,k) + & 
			 ALF(j,anc(i))*dm(j)*real(bbbb_mod(j,k)*long_table(mm(j),i)) 
		enddo   
	    enddo
	enddo 
!
	do k=0, NN
	do j=1, jmax 
	    kkkk(j,k)=0.  
	    do i=1, np 
	    if(wet(i)==1) kkkk(j,k) = & 
			  kkkk(j,k) + x(i,k)*alf(j,anc(i))*conjg(long_table(mm(j),i))  
	    enddo
	enddo
	enddo 
	kkkk(:,:)=kkkk(:,:)/float(np) 
!
!
! --- Solving for arrays 'Z' and 'S' 
	Z(:,:,is) = HHHH(:,:) + KKKK(:,:) 
	S(:,:)    = AAAA_MOD(:,:) + SE(:,:) + BBBB_MOD(:,:)

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
! --- Array "B" for vertical displacement 
	bbbb(:,:)=0.
 	do j=1, jmax 
		do k=0,NN
		bbbb(j,k) = EU(ll(j))*Z(j,k,SMAX)       
		do p=0, k
		if(p==0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-0.           )*BETAU(ll(j),k-p)
		if(p/=0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-Z(j,p-1,SMAX))*BETAU(ll(j),k-p)
		enddo 
        	bbbb(j,k)= RHOW_o_RHOE_X3*bbbb(j,k)
        	enddo
	enddo
!
! --- Array "A" for vertical displacement 
 	aaaa(:,:)=0.
 	do j=1, jmax 
	 	do k=0,NN       
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
! --- Array "B" for Geoid heigth   
	bbbb(:,:)=0.
 	do j=1, jmax 
		do k=0,NN
		bbbb(j,k) = EN(ll(j))*Z(j,k,SMAX)       
		do p=0, k
		if(p==0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-0.           )*BETAN(ll(j),k-p)
		if(p/=0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-Z(j,p-1,SMAX))*BETAN(ll(j),k-p)
		enddo 
        	bbbb(j,k)= RHOW_o_RHOE_X3*bbbb(j,k)
        	enddo
	enddo
!
! --- Array "A" for Geoid heigth  
 	aaaa(:,:)=0.
 	do j=1, jmax 
	 	do k=0,NN       
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
! --- Adding a constant to geoid undulations 
	N(1,:) = N(1,:) +  SE(1,:) - AAVV(:) - BBVV(:)  
!
! --- Geoid undulations (previous formulation) 
!       N(:,:) = S(:,:) + U(:,:)
!
!
! --- Array "B" for horizontal displacement 
	bbbb(:,:)=0.
 	do j=1, jmax 
		do k=0,NN
		bbbb(j,k) = EV(ll(j))*Z(j,k,SMAX)       
		do p=0, k
		if(p==0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-0.           )*BETAV(ll(j),k-p)
		if(p/=0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-Z(j,p-1,SMAX))*BETAV(ll(j),k-p)
		enddo 
        	bbbb(j,k)= RHOW_o_RHOE_X3*bbbb(j,k)
        	enddo
	enddo
!
! --- Array "A" for horizontal displacement 
 	aaaa(:,:)=0.
 	do j=1, jmax 
	 	do k=0,NN       
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
! --- Array "B" for FA Gravity 
	bbbb(:,:)=0.
 	do j=1, jmax 
		do k=0,NN
		bbbb(j,k) = EFA(ll(j))*Z(j,k,SMAX)       
		do p=0, k
		if(p==0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-0.           )*BETAFA(ll(j),k-p)
		if(p/=0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-Z(j,p-1,SMAX))*BETAFA(ll(j),k-p)
		enddo 
        	bbbb(j,k)= RHOW_o_RHOE_X3*100.*1000.*(GRA_REF/RAD_REF)*bbbb(j,k)      
        	enddo
	enddo
!
! --- Array "A" for FA Gravity 
 	aaaa(:,:)=0.
 	do j=1, jmax 
	 	do k=0,NN       
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
! --- Array "B" for SS Gravity 
	bbbb(:,:)=0.
 	do j=1, jmax 
		do k=0,NN
		bbbb(j,k) = ESS(ll(j))*Z(j,k,SMAX)       
		do p=0, k
		if(p==0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-0.           )*BETASS(ll(j),k-p)
		if(p/=0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-Z(j,p-1,SMAX))*BETASS(ll(j),k-p)
		enddo 
        	bbbb(j,k)= RHOW_o_RHOE_X3*100.*1000.*(GRA_REF/RAD_REF)*bbbb(j,k)      
        	enddo
	enddo
!
! --- Array "A" for SS Gravity 
 	aaaa(:,:)=0.
 	do j=1, jmax 
	 	do k=0,NN       
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   SS Gravity DRAFT 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   SS Gravity DRAFT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   SS Gravity DRAFT
!
! --- Array "B" for SS Gravity 
	bbbb(:,:)=0.
 	do j=1, jmax 
		do k=0,NN
		bbbb(j,k) = ET(ll(j))*Z(j,k,SMAX)       
		do p=0, k
		if(p==0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-0.           )*BETAT(ll(j),k-p)
		if(p/=0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-Z(j,p-1,SMAX))*BETAT(ll(j),k-p)
		enddo 
        	bbbb(j,k)= RHOW_o_RHOE_X3 *bbbb(j,k)      
        	enddo
	enddo
!
! --- Array "A" for SS Gravity 
 	aaaa(:,:)=0.
 	do j=1, jmax 
	 	do k=0,NN       
 	 	aaaa(j,k) = ET(ll(j))*IIII(j,k) 
 	 	do p=0, k
 	 	if(p==0) aaaa(j,k) = aaaa(j,k) - (IIII(j,p)-0.         )*BETAT(ll(j),k-p) 	 
 	 	if(p/=0) aaaa(j,k) = aaaa(j,k) - (IIII(j,p)-IIII(j,p-1))*BETAT(ll(j),k-p)
 		enddo 
 		aaaa(j,k)= RHOI_o_RHOE_X3 *aaaa(j,k)
 		enddo
	enddo
!
! --- SS Gravity gravity anomaly          (in units of "MILLIGAL")
	TT(:,:) = aaaa(:,:) + bbbb(:,:)
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
! --- Total (ice+ocean) load             			         <<<< new 
	Out_Filename=filelt 
 	open(3,file=Out_Filename,status='unknown',form='unformatted') 
        write(3) LO_T 
	close(3)
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Write(*,*) "    +++ SLE solved <3 +++" 
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
	End program SLE
!
!
!
!
!
!
