!
!*******************************
	MODULE MAIN_MODULE 
!*******************************
!
!
include 'alma.inc'
!
!
! This module provides a number of declarations 
!	
! Model parameters 
	TYPE(FM) R(0:NLA+2)    ! radii of the interfaces 
        TYPE(FM) RK(0:NLA+2)   ! radii of the interfaces (km) 
	TYPE(FM) RHO(0:NLA+1)  ! density of core (0), mantle, and litho (NLA+1) 
	TYPE(FM) RIG(0:NLA+1)  ! rigidity of core (0), mantle, and litho (NLA+1) 
	TYPE(FM) LAM(0:NLA+1)  ! 2nd Lame constant lambda 
	TYPE(FM) VIS(0:NLA+1)  ! viscosity of mantle layers 
	TYPE(FM) GRA(0:NLA+2)  ! gravity at the radii
	TYPE(FM) AC(0:NLA+1)   ! <<a>> constants 
	TYPE(FM) EMASS    ! Earth mass  
	TYPE(FM) LT       ! Thickness of the lithosphere 
	TYPE(FM) DELTA    ! Thickness of each mantle layer 
	TYPE(FM) NT 				! Newton constant 
        TYPE(FM) RA0, RHO0, RIG0, LAM0, T0	! normalization pivots 	
!
! gaver variables 
 	type(fm) gh(0:ng,1:2*ng)
 	type(fm) gl(0:ng,1:2*ng)
 	type(fm) gk(0:ng,1:2*ng)
!
	real ggh(l1:l2,P+1,0:ng), &
	     ggl(l1:l2,P+1,0:ng), &
	     ggk(l1:l2,P+1,0:ng)
!	
! salzer variables
 	integer, parameter :: ms_max=ng    ! Order of the Salzer recurrence  
	type(fm) FH(ms_max),    FL(ms_max),    FK(ms_max) 
	real	 FH_SP(l1:l2,P+1,ms_max), & 
	         FL_SP(l1:l2,P+1,ms_max), & 
		 FK_SP(l1:l2,P+1,ms_max) 
!
! time variable 
 	type(fm) time(p+1) 
!
! Heaviside variables 
 	type(fm) hd(p+1), ld(p+1), kd(p+1) 
	type(fm) rec_h(p+1), rec_l(p+1), rec_k(p+1) 	
!
! Love numbers 
 	type(fm) acca,   elle,   kapp
!
! Elastic and fluid Ln's
	real acca_e(l1:l2), elle_e(l1:l2), kapp_e(l1:l2)
	real acca_f(l1:l2), elle_f(l1:l2), kapp_f(l1:l2)	
!
! Other variables
        real ssm, ssc, time_sp     
!
! Earth radius and CMB radius of PREM 
 	real, parameter :: ea=6.371e6	       
 	real, parameter :: ec=3.480e6	
!
! Working directory 
	character*100 wdir 
!
	end MODULE MAIN_MODULE 
!
!
!
   PROGRAM ALMA
!
!******************
!  "Main program"
!******************        
!
   USE MAIN_MODULE
   IMPLICIT NONE	
!	
   integer i, j, l, n 
!
!
   Open (33,file='working-directory.txt',status='unknown')
   read (33,'(a100)')wdir 
   close(33)
!
   write(*,*) '    - alma.f90: Initializing the multiprecision routines' 
   CALL FM_SET(nsd) 
!	
   Write(*,*) '    - alma.f90: Building model: ', trim(adjustl(mod2))   
   call NORM 
!
   Write(*,*) '    - alma.f90: Building the time steps'  
   call time_steps    
!
   Write(*,*) '    - alma.f90: Setting the working directory'  


! # Loop on the harmonic degrees 
!
do 100 l=l1, l2, ls
!
!
   if(mod(l,10)==0.or.l==l1.or.l==l2) write(*,*) '    - alma.f90: Harmonic degree =', l, 'of', l2  
!
   if(ih==1.or.il==1.or.ik==1) then 
!
! # Counts time 
   call time_c(ssm)
!
! # Time do-loop 
!--------------- 
do 200 j=1, p+1 
!---------------
  if(mod(j,50)==0)write(*,*) ' time = ', j, 'of', p+1
!
! # Computes the Gaver functionals 
 call gaver_coeff(l,time(j))   	
!
! # Stores the functional 
 do n=1,ng 
	if(ih==1) ggh(l,j,n)=to_sp(gh(n,n))
	if(il==1) ggl(l,j,n)=to_sp(gl(n,n))
	if(ik==1) ggk(l,j,n)=to_sp(gk(n,n))
 enddo
!
if(isalz==1) then
!
! # "Salzer acceleration"  
 call salzer_coeff
!
! # Stores the accelerated Ln's 
 do n=1,ng 
	if(ih==1) fh_sp(l,j,n)=to_sp(fh(n))
	if(il==1) fl_sp(l,j,n)=to_sp(fl(n))
	if(ik==1) fk_sp(l,j,n)=to_sp(fk(n))
 enddo
!
	     endif
!
200 CONTINUE    
!
!
 call time_c(ssc)

! write(*,*) ' Elapsed time (s) =', ssc-ssm
! write(*,*) ' Time for each step (s) =', (ssc-ssm)/float(p+1)
!
	     endif	
!
100 continue
!
!
 Write(*,*) '    - alma.f90: Writing the Love numbers ...'
! 
! Opens the output files 
 Call oc (0) 
!
do 300 j=1, p+1 
!
time_sp=to_sp(time(j))
if(isalz==0) then
if(ih==1) write(71,'(1024(e19.8))') time_sp, (ggh(l,j,ng),l=l1,l2,ls)
if(il==1) write(72,'(1024(e19.8))') time_sp, (ggl(l,j,ng),l=l1,l2,ls)
if(ik==1) write(73,'(1024(e19.8))') time_sp, (ggk(l,j,ng),l=l1,l2,ls)
endif
if(isalz==1) then
if(ih==1) write(71,'(1024(e19.8))') time_sp, (FH_sp(l,j,ng),l=l1,l2,ls)
if(il==1) write(72,'(1024(e19.8))') time_sp, (FL_sp(l,j,ng),l=l1,l2,ls)
if(ik==1) write(73,'(1024(e19.8))') time_sp, (FK_sp(l,j,ng),l=l1,l2,ls)
endif    
!
300 continue
!
! Closes the output files 
 Call oc (1)
!
! 
 
 
 
 
! 
! 	
 END PROGRAM ALMA 
!
!
!
!
!
   SUBROUTINE OC(SWITCH) 
!
!-------------------------------------
! Opens and closes the output files... 
!-------------------------------------
!
   USE MAIN_MODULE 
   IMPLICIT NONE
   integer deg_min, deg_max, switch 
!
   deg_min=l1; deg_max=l2
!
!
! For switch=0, the files are open ----------
   if(switch==0) then 
! 
!write(*,'(a100)')wdir
!
	if(ih==1) open(71,file=trim(adjustl(wdir))//'/ALMA/'//'h.dat',status='unknown') 
	if(il==1) open(72,file=trim(adjustl(wdir))//'/ALMA/'//'l.dat',status='unknown') 
	if(ik==1) open(73,file=trim(adjustl(wdir))//'/ALMA/'//'k.dat',status='unknown')
! 
 	if(ih==1) Call header (deg_min,deg_max,71)
 	if(il==1) Call header (deg_min,deg_max,72)
 	if(ik==1) Call header (deg_min,deg_max,73)
!      
   endif
!
!
! For switch=1, the files are closed -----------
   if(switch==1) then 
!
   if(ih==1) write(71,*) ' # End of file'  
   if(il==1) write(72,*) ' # End of file' 
   if(ik==1) write(73,*) ' # End of file'      
!	
if(ih==1) close(71) 
if(il==1) close(72) 
if(ik==1) close(73)
!  
   endif
!
   END SUBROUTINE OC 
!   
!
!
!
!
 SUBROUTINE TIME_C(SECONDS)
!
!------------------------------------------- 
! Counts the seconds elapsed since midnight
!-------------------------------------------
!
 implicit NONE
 character*20 date, timc
 integer iu, i1, i2, i3, i4 
 real seconds 
!
 call date_and_time (date, timc)

! write(*,*) date(1:4),'.',date(5:6),'.',date(7:8),' ', & 
!            timc(1:2),'.',timc(3:4),'.',timc(5:6),'.',timc(8:10) 
!
 iu=97 
 open(iu,file='fake.dat',status='unknown') 
 write(iu,'(a2)') timc(1:2); write(iu,'(a2)') timc(3:4)
 write(iu,'(a2)') timc(5:6); write(iu,'(a3)') timc(8:10)
 close(iu)  
 open(iu,file='fake.dat',status='unknown') 
 read(iu,*) i1; read(iu,*) i2; read(iu,*) i3; read(iu,*) i4  
 close(iu) 
!
! # Seconds since midnight 
  seconds=float(i1)*3600.+float(i2)*60.+float(i3)+ float(i4)/1000.
! 
 END SUBROUTINE TIME_C
!
!
!
!
!
SUBROUTINE HEADER (Lmin,Lmax,IU)
!
!----------------------------------------------------
!  Writes an header on the Love numbers output files    
!----------------------------------------------------
!
	USE MAIN_MODULE 
	implicit none
	integer lmin, lmax, iu 
        character*20 date, timc
!
 	call date_and_time (date, timc)
!
		write(iu,*) ' # '  
if(           iload==1) write(iu,*) ' # -------------------' 
if(iu==71.and.iload==1) write(iu,*) ' # h load Love number '
if(iu==72.and.iload==1) write(iu,*) ' # l load Love number '
if(iu==73.and.iload==1) write(iu,*) ' # k load Love number '
if(           iload==1) write(iu,*) ' # -------------------' 
if(           iload==0) write(iu,*) ' # --------------------' 
if(iu==71.and.iload==0) write(iu,*) ' # h tidal Love number '
if(iu==72.and.iload==0) write(iu,*) ' # l tidal Love number '
if(iu==73.and.iload==0) write(iu,*) ' # k tidal Love number '
if(           iload==0) write(iu,*) ' # --------------------' 
		write(iu,*) ' # '  
!	   
		write(iu,*) ' # Done by alma.f90 (SELEN port) on: ', & 
 				date(1:4),'.',date(5:6),'.',date(7:8),' time: ',&
		        	timc(1:2),'.',timc(3:4),'.',timc(5:6)
if(th=='h') write(iu,*) ' # Heaviside loading'
if(th=='r') write(iu,*) ' # Ramp loading'
		write(iu,*) ' # Range of harmonic degrees:', lmin,'-',lmax 
		write(iu,*) ' # Number of digits:', nsd 
		write(iu,*) ' # Order of the Gaver sequence:', ng    
if(isalz==0)    write(iu,*) ' # The sequence is not accelerated'
if(isalz==1)    write(iu,*) ' # The sequence is accelerated' 
if(imode==1)    write(iu,*) ' # PREM volume-averaged model'
if(imode==2)    write(iu,*) ' # User-supplied stratification, model: ', mod2
if(imode==3)    write(iu,*) ' # User-supplied stratification, model: ', mod3
if(imode==3)	write(iu,*) ' # General viscoelatic rheology (***** nla=2 *****)'
if(imode==3) then 
if(ire==1)write(iu,*)' # Hooke (H)'
if(ire==2)write(iu,*)' # Newton (N)'
if(ire==3)write(iu,*)' # Maxwell (M)'
if(ire==4)write(iu,*)' # Kelvin-Voigt (KV)'
if(ire==5)write(iu,*)' # Standar Linear Solid, Maxwell type (SLS-M)'
if(ire==6)write(iu,*)' # Standar Linear Solid, Kelvin type (SLS-K)'
if(ire==7)write(iu,*)' # Burgers (B)'
if(ire==8)write(iu,*)' # Caputo (C)'
if(ire==9)write(iu,*)' # Weichert (W)'
          endif
                write(iu,*) ' # Number of mantle layers:', NLA
!
!
if(imode==1)    write(iu,*) ' # Kind of viscosity profile:', kv 
	        write(iu,*) ' # Time interval is from ',m1,'to ',m2, 'kyrs' 
	        write(iu,*) ' # Time points:', P+1 
		write(iu,*) ' # ' 
if(iu==71) 	write(iu,*) ' #    time (kyrs)      h(t) for degrees l=', lmin,'...',lmax
if(iu==72) 	write(iu,*) ' #    time (kyrs)      l(t) for degrees l=', lmin,'...',lmax
if(iu==73) 	write(iu,*) ' #    time (kyrs)      k(t) for degrees l=', lmin,'...',lmax
		write(iu,*) ' # '  
!
END SUBROUTINE HEADER 
!
!
!
!
!
!
SUBROUTINE LOVE_NUMBERS (L, X, FLAG) 
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Computes the Love numbers at degree 'l' at frequency 'x' 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
	use FMZM
!
	USE MAIN_MODULE 
	implicit NONE
!
 	INTEGER I, J, K, L, FLAG, NNN, MMM1, MMM2, MMM
        INTEGER, PARAMETER :: N=6 
	INTEGER INDX(N/2)  
	TYPE(FM) PI
 	TYPE(FM) NUL, ONE, X, D, MUSS, LASS 
 	TYPE(FM) PJ1(N/2,N), PJ2(N/2,N), CMB_A(N,N/2)	
	TYPE(FM) PJ1_WWW_CMB(N/2,N/2), PJ2_WWW_CMB(N/2,N/2)
 	TYPE(FM) IDE(N,N), BC(N/2,1), CC(N/2,1), XA(N/2,1)  
        TYPE(FM) WWW(N,N), DDDD(N,N), IIII(N,N), XXX(N,N)
!
	TYPE(FM) RIG2, VIS2, TAUC 
!
	TYPE(FM) MUW, EXPO, RRR(100), ETA(100)  
	
	TYPE(FM) ARGO, EXPP, EXPM
!
        TYPE(FM) ESSETAUL, TAUL, THF	
!	
!
!
        pi = 2.0*asin(1.0)
!
!
! p1 and p2 projectors 
! 
	nul=to_fm(0.); one=to_fm(1.)
!
	If     (l.ne.1) then 
		pj1(:,:)=nul; pj1(1,1)=one; pj1(2,2)=one; pj1(3,5)=one 	
		pj2(:,:)=nul; pj2(1,3)=one; pj2(2,4)=one; pj2(3,6)=one 
	Elseif (l.eq.1) then 
		pj1(:,:)=nul; pj1(1,1)=one; pj1(2,2)=one; pj1(3,5)=one 	   
		pj2(:,:)=nul; pj2(1,3)=one; pj2(2,4)=one; pj2(3,5)=one 
	Endif	
!
!
! Initializing the W propagator and the unit matrix 
!
	www(:,:)=nul 
	ide(:,:)=nul 
	do i=1,6
		www(i,i)=one
		ide(i,i)=one 
	enddo
!
! W propagator 
	do 100 j=nla+1,1,-1
!
! The lithosphere is always elastic 
	if(j==nla+1) muss=rig(j)
!
! Maxwell rheology 	
         muss =rig(j)*x/(x+rig(j)/vis(j))       ! Maxwell 
!
!
  If(imode==3) then 
!
!
! ********************************************************************
! This part is to test various rheologies - for one layer mantle
! If you want to use this part, just delete the comments and choose
! the rheological model setting the 'ire' switch. Here the rheologies
! are applied to layer J=1. This is the whole mantle in the case
! nla=1. 
! ******************************************************************** 
!
	IF(J==2) THEN     ! J=1 ---- > LOWER MANTLE ; J=2 ---- > UPPER MANTLE 
!		
	if(ire==1) muss= rig(j)				! Hooke 
!
	if(ire==2) muss= x*vis(j)			! Newton 
!	
	if(ire==3) muss= rig(j)*x/(x+rig(j)/vis(j))	! Maxwell 							         
!
	if(ire==4) then 				! Kelvin-Voigt   
	rig2=0.1*rig(j) 
	vis2=0.1*vis(j) 	
	muss= rig2 + x*vis2 
	endif				
!
	if(ire==5) then					! Standard Linear Solid (M)
	rig2=0.1*rig(j) 
	vis2=0.1*vis(j) 
	muss= (x*(rig(j)+rig2)+rig2*rig(j)/vis(j))/(x+rig(j)/vis(j)) 
	endif
!
	if(ire==6) then					! Standard Linear Solid (V)
        rig2=0.1*rig(j) 
	vis2=0.1*vis(j) 
	muss= (x*rig(j)+rig(j)*rig2/vis2)/(x+(rig(j)+rig2)/rig2)    
	endif
!
	if(ire==7) then					! Burgers    ! UPDATED on June 10, 2009.- 
	rig2=0.3*rig(j) 
	vis2=vis(j)  
 	muss=rig(j)*x*(x+rig2/vis2)/&
 		(x**2 + x*((rig(j)+rig2)/vis2 + rig(j)/vis(j)) + &
 		rig(j)*rig2/vis(j)/vis2 )	
	endif	
!
	if(ire==8) then					! Caputo
	tauc=6.3E9/60./24./365./1000.
	muss=rig(j)*(x*tauc)**0.85/(1.+(x*tauc)**0.85) 
	endif	
!	
!
	if(ire==9) then					! Weichert (generalized Maxwell)
!
!
! Number of Maxwell elements 
  	NNN=30
!	
! Wiechert shear modulus (elastic element)		
	MUW=0.0 
!
! Shear moduli of Maxwell elements	
	DO I=1,NNN
	   RRR(I)=RIG(J)	
	   ETA(I)=VIS(J)/FLOAT(I)/FLOAT(I)	
	ENDDO
!
! Option #1:Complex ridigity for a CONTINUOUS Wiechert model 
!
	ARGO=PI*SQRT(X*VIS(J)/RIG(J))
	EXPP=EXP(+ARGO)
	EXPM=EXP(-ARGO)
!		
	MUSS= MUW + 0.5*RIG(J)*(-1.0 + ARGO*(EXPP+EXPM)/(EXPP-EXPM))
!	
!
! Option #2: Complex ridigity for a DISCRETE Wiechert model of NNN elements 
!
!	MUSS=MUW
!	DO I=1, NNN
!		MUSS=MUSS + RRR(I)*X/(X+RRR(I)/ETA(I))	
!	ENDDO
!
	endif
!
!
!
	ENDIF   ! on J=!  
!
!
	ENDIF
!
!       End of the test with various rheologies
!===========================================================================================             	
!		
!
! Direct and Inverse propagators... 
!
		call direct (l,r(j+1),rho(j),ac(j),gra(j+1),muss,DDDD)
		call inverse(l,r(j  ),rho(j),ac(j),gra(j  ),muss,IIII)	
!
		xxx=matmul(www,matmul(DDDD,IIII))
		www=xxx
!
 100 Continue 
!
!
! CMB interface matrix, at the CMB   
	call cmb (l, r(1), rho(0), ac(0), cmb_a)

! Surface boundary conditions 
        call surf_bc (l, r(nla+2), nt, gra(nla+2), bc)
!
! Computing the Ln's 
        pj2_www_cmb = matmul(matmul(pj2,www),cmb_a)
!
        call ludcmp(pj2_www_cmb, N/2, N/2, indx, d)
        call lubksb(pj2_www_cmb, N/2, N/2, indx, bc)
!
        pj1_www_cmb = matmul(matmul(pj1,www),cmb_a)
!	
	xa = matmul(pj1_www_cmb,bc) 
!
!
! Love numbers's h, l, and k
!
	if(flag==1) then  
! 
        if      (th=='r') then                ! Ramp loading 
! 	
	TAUL=to_fm('1.0')                     ! loading period (kyrs) 
!
        THF=((ONE - EXP(-X*TAUL))/(X*TAUL))   ! TH factor 
!
	if(ih==1) acca = (xa(1,1)*emass/r(nla+2)/x)*THF		  
	if(il==1) elle = (xa(2,1)*emass/r(nla+2)/x)*THF			  
	if(ik==1) kapp = ((-one - (emass/r(nla+2)/gra(nla+2))*xa(3,1))/x)*THF	
!
        elseif  (th=='h') then  ! Heaviside loading 
!
 	if(ih==1) acca = (xa(1,1)*emass/r(nla+2)/x)		  
 	if(il==1) elle = (xa(2,1)*emass/r(nla+2)/x)			  
 	if(ik==1) kapp = ((-one - (emass/r(nla+2)/gra(nla+2))*xa(3,1))/x)
!
        endif
!
endif 
!
END SUBROUTINE love_numbers
!
!
!
!
!
!
SUBROUTINE SALZER_COEFF 
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Computes the Salzer accelerators (Valko & Abate, 2002) 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
USE FMZM
USE MAIN_MODULE 
implicit NONE
!
 integer kk, nn 
! 
 TYPE(FM) W1, W2, W3, W4, Weight, fkk, fms
!
 do nn=1, ms_max 
!
if(ih==1) FH(nn)=to_fm(0.0)
if(il==1) FL(nn)=to_fm(0.0)
if(ik==1) FK(nn)=to_fm(0.0)   
!
 do kk=1, NN  
!
 	CALL FMI2M (kk,Fkk)
 	CALL FMI2M (nn,fms) 
 	call FMIPWR (to_fm(-1.0), kk+nn, W1)    
 	CALL FMIPWR (fkk,nn,W2) 
 	CALL FMFACT (fms,W3)   
 	CALL FMCOMB (fms,fkk,W4)    
 	Weight = W1*W2*W4/W3
!
if(ih==1) FH(nn) = FH(nn) + Weight*gh(kk,kk)
if(il==1) FL(nn) = FL(nn) + Weight*gl(kk,kk) 
if(ik==1) FK(nn) = FK(nn) + Weight*gk(kk,kk) 
! 
 enddo 
 enddo
!
END SUBROUTINE SALZER_COEFF 
!
!
!
!
!
SUBROUTINE GAVER_COEFF (L, EPOCH) 
!
!*****************************************************
! Computes the Gaver functionals (Valko & Abate, 2002) 
!*****************************************************
!
use FMZM
USE MAIN_MODULE 
implicit NONE
integer kk, mm, j, l   
type(fm) log2, fmm, fkk, ratio, epoch, argument 
! 
!
 gh(:,:)=to_fm(0.0)
 gl(:,:)=to_fm(0.0)
 gk(:,:)=to_fm(0.0)  
!
kk=0
do mm=1,2*ng
!
	CALL FMI2M(MM,FMM)
	CALL FMLN(TO_FM(2.),LOG2) 
! 
 	ARGUMENT=FMM*LOG2/EPOCH 
!
 	call love_numbers(l, ARGUMENT, 1) 
!
	if(ih==1) gh(kk,mm) = argument*acca 
	if(il==1) gl(kk,mm) = argument*elle
	if(ik==1) gk(kk,mm) = argument*kapp
!	
enddo
! 
do kk=1, ng
! 
	do mm=1, 2*ng-kk  
!
	CALL FMI2M(MM,FMM)
	CALL FMI2M(KK,FKK)
!
	RATIO=FMM/FKK
!
	if(ih==1) gh(kk,mm)=(1.+RATIO)*gh(kk-1,mm) - RATIO*gh(kk-1,mm+1) 
	if(il==1) gl(kk,mm)=(1.+RATIO)*gl(kk-1,mm) - RATIO*gl(kk-1,mm+1) 
	if(ik==1) gk(kk,mm)=(1.+RATIO)*gk(kk-1,mm) - RATIO*gk(kk-1,mm+1) 
!
	enddo
!
enddo
!
END SUBROUTINE gaver_coeff 
!
!
!
!
!
SUBROUTINE INVERSE (N, R, RHO, AC, GRA, MU, INVE) 
! 
!------------------------------------- 
! # Inverse of the fundamental matrix
!------------------------------------- 
!
USE FMZM
implicit none
!
integer, parameter :: nn=6 
integer  i, j, n  
TYPE(FM) nul, gra, r, rho, ac, mu, fn 
TYPE(FM) h1, h2, h3, h4, h5
TYPE(FM) v(24),a(nn,nn),b(nn,nn),inve(nn,nn) 
!
!
!
nul=to_fm(0.0)
inve(:,:)=nul; a(:,:)=nul; b(:,:)=nul 
!
fn = float (n)  
!
h1 = fn+1.  
h2 = 1.+2.*fn  
h3 = 2.*fn-1.  
h4 = 3.+2.*fn  
h5 = fn-1.  
!
v(1) = (fn+1.)/(1.+2.*fn)  
v(2) = -2.*(fn+1.)*(fn+2.)/(1.+2.*fn)
v(3) = fn*(fn+1.)/2./(2.*fn-1.)/(1.+2.*fn)  
v(4) = - (fn*fn+3.*fn-1.)*fn/(1.+2.*fn)/(2.*fn-1.)  
v(5) = fn/(1.+2.*fn)  
v(6) = 2.*fn*(fn-1.)/(1.+2.*fn)
v(7) = fn*(fn+1.)/2./(1.+2.*fn)/(3.+2.*fn)  
v(8) = (fn+1.)*(fn*fn-fn-3.)/(1.+2.*fn )/(3.+2.*fn)  
v(9) = -2.*fn*(fn+1.)*(2.+fn)/(1.+2.*fn)  
v(10) = -(fn-1.)*fn*(fn+1.)*(fn+1.)/(2.*fn-1.)/ (1.+2.*fn)  
v(11) = -2.*(fn-1.)*fn*(fn+1.)/(1.+2.*fn)  
v(12) = -fn*fn*(fn+1.)*(2.+fn)/(1.+2.*fn )/(3.+2.*fn)  
!
v(13) = - v (1)  
v(14) = - v (3)  
v(15) = - v (5)  
v(16) = - v (7) 
!
v(17) = -fn*(fn+1.)/(1.+2.*fn)  
v(18) = (2.-fn)*(fn+1.)*fn/2./(2.*fn-1.)/(1.+2.*fn )  
v(19) = fn*(fn+1.)/(1.+2.*fn)  
v(20) = fn *(fn+1.)*(3.+fn)/2./(1.+2.*fn)/(3.+2.*fn) 
!
v(21) = - v (1)  
v(22) = - v (3)  
v(23) = - v (5)  
v(24) = - v (7)  
!
!
! 
! Array A - 1st column  
 call FMIPWR(r,-n-1,a(1,1)); 	a(1,1)=+a(1,1)*v(2) 			
 call FMIPWR(r,-n+1,a(2,1)); 	a(2,1)=-a(2,1)*v(4)  			
 call FMIPWR(r,-n+1,a(3,1)); 	a(3,1)=+a(3,1)*3.*ac/(2.*fn+1.)                
 call FMIPWR(r,+n  ,a(4,1)); 	a(4,1)=+a(4,1)*v(6) 			
 call FMIPWR(r,+n+2,a(5,1)); 	a(5,1)=-a(5,1)*v(8) 			
 call FMIPWR(r,+n+2,a(6,1)); 	a(6,1)=-a(6,1)*3.*ac/(2.*fn+1.) 	
!
!
! Array A - 2nd column  
 call FMIPWR(r,-n-1,a(1,2)); 	a(1,2)=-a(1,2)*v(9) 			
 call FMIPWR(r,-n+1,a(2,2)); 	a(2,2)=+a(2,2)*v(10) 			
 				a(3,2)=nul
 call FMIPWR(r,+n  ,a(4,2));    a(4,2)=-a(4,2)*v(11) 			
 call FMIPWR(r,+n+2,a(5,2));    a(5,2)=+a(5,2)*v(12) 			
				a(6,2)=nul
!
!
! Array A - 5th column
 do i=1,6
 	a(i,5)=to_fm(0.0) 	
 enddo
 call FMIPWR(r,+n+1,a(6,5));    a(6,5)=-a(6,5) 				
!
!
! Array A - 6th column  
 				a(1,6)=nul      
				a(2,6)=nul
 call FMIPWR(r,-n+1,a(3,6));    a(3,6)=-a(3,6)/(2.*fn+1.) 		
 				a(4,6)=nul       
				a(5,6)=nul 
 call FMIPWR(r,+n+2,a(6,6));    a(6,6)=+a(6,6)/(2.*fn+1.) 
!
!
! Array B - 1st column 
 call FMIPWR(r,-n  ,b(1,1));    b(1,1)=+b(1,1)*rho*gra*v(1)  		
 call FMIPWR(r,-n+2,b(2,1));	b(2,1)=-b(2,1)*rho*gra*v(3)		
				b(3,1)=nul
 call FMIPWR(r,+n+1,b(4,1));    b(4,1)=+b(4,1)*rho*gra*v(5)  		
 call FMIPWR(r,+n+3,b(5,1));	b(5,1)=-b(5,1)*rho*gra*v(7)		
				b(6,1)=nul
!
! Array B - 2nd column 
 do i=1,6
 	b(i,2)=to_fm(0.0) 	
 enddo
!
!
! Array B - 3rd column 
 call FMIPWR(r,-n  ,b(1,3));    b(1,3)=+b(1,3)*v(13)     		
 call FMIPWR(r,-n+2,b(2,3));	b(2,3)=-b(2,3)*v(14)  			
				b(3,3)=nul
 call FMIPWR(r,+n+1,b(4,3));    b(4,3)=+b(4,3)*v(15)  			
 call FMIPWR(r,+n+3,b(5,3));	b(5,3)=-b(5,3)*v(16)			
				b(6,3)=nul
!
!
! Array B - 4th column
 call FMIPWR(r,-n  ,b(1,4));    b(1,4)=-b(1,4)*v(17)  			
 call FMIPWR(r,-n+2,b(2,4));	b(2,4)=+b(2,4)*v(18)  			
				b(3,4)=nul
 call FMIPWR(r,+n+1,b(4,4));    b(4,4)=-b(4,4)*v(19) 			
 call FMIPWR(r,+n+3,b(5,4));	b(5,4)=+b(5,4)*v(20)			
				b(6,4)=nul		
!
!
! Array B - 5th column
 call FMIPWR(r,-n  ,b(1,5));    b(1,5)=-b(1,5)*rho*v(21)   		
 call FMIPWR(r,-n+2,b(2,5));	b(2,5)=+b(2,5)*rho*v(22)   		
				b(3,5)=nul
 call FMIPWR(r,+n+1,b(4,5));    b(4,5)=-b(4,5)*rho*v(23) 		
 call FMIPWR(r,+n+3,b(5,5));	b(5,5)=+b(5,5)*rho*v(24) 		
				b(6,5)=nul
!
do i=1, 6; do j=1, 6
	inve(i,j)  = a(i,j) + b(i,j)/mu  
enddo; enddo
!  
end subroutine inverse
!
!
!
!
!
SUBROUTINE DIRECT (N, R, RHO, AC, GRA, MU, DIRE)  
!
!-------------------------------
! # Direct fundamental matrix  
!-------------------------------
!
USE FMZM
implicit none
!
integer, parameter :: NN=6   
integer  i, j, n  
TYPE(FM) nul, gra, r, rho, ac, mu, fn  
TYPE(FM) a1, a2, a3, a4, b1, b2, b3, b4 
TYPE(FM) c1, c2, c3, c4, d1, d2, d3, d4
TYPE(FM) a(nn,nn), b(nn,nn), dire(nn,nn)
!
!
nul=to_fm(0.0)
dire(:,:)=nul; a(:,:)=nul; b(:,:)=nul 
!
!
fn = float (n)  
!
 a1 = fn/(2.*(2.*fn+3.))
 a2 = 1.  
 a3 = (fn+1.)/(2.*(2.*fn-1.))  
 a4=  1.
 b1 = (fn+3.)/(2.*(2.*fn+3.)*(1.*fn+1.))
 b2 = 1./fn  
 b3 = (-fn+2.)/(2.*fn*(2.*fn-1.))  
 b4 = - 1./(fn+1.)  
 c1 = (fn*fn-fn-3.)/(2.*fn+3.)  
 c2 = 2.*(fn-1.)  
 c3 = (-fn*fn-3.*fn+1.)/(2.*fn-1.)  
 c4 = - 2.*(fn+2.)  
 d1 = fn*(fn+2.)/((2.*fn+3.)*(fn+1.))  
 d2 = 2.*(fn-1.)/fn  
 d3 = (fn*fn-1.)/(fn*(2.*fn-1.))  
 d4 = 2.*(1.*fn+2.)/(fn+1.)  
!
! 
! Array A - 1st row 
 call FMIPWR(r,+n+1,a(1,1)); 	a(1,1)=a(1,1)*a1 
 call FMIPWR(r,+n-1,a(1,2)); 	a(1,2)=a(1,2)*a2 
 				a(1,3)=nul   
 call FMIPWR(r,-n,  a(1,4)); 	a(1,4)=a(1,4)*a3 
 call FMIPWR(r,-n-2,a(1,5)); 	a(1,5)=a(1,5)*a4 
 				a(1,6)=nul   
! 
! Array A - 2nd row 
 call FMIPWR(r,+n+1,a(2,1)); 	a(2,1)=a(2,1)*b1 
 call FMIPWR(r,+n-1,a(2,2)); 	a(2,2)=a(2,2)*b2  
				a(2,3)=nul   
 call FMIPWR(r,-n  ,a(2,4)); 	a(2,4)=a(2,4)*b3  
 call FMIPWR(r,-n-2,a(2,5)); 	a(2,5)=a(2,5)*b4 
 				a(2,6)=nul
!
!  Array A - 3rd row
 a(3,1)=a(1,1)*rho*gra
 a(3,2)=a(1,2)*rho*gra  
 call FMIPWR(r,+n  ,a(3,3)); 	a(3,3)=-rho*a(3,3) 
 a(3,4)=a(1,4)*rho*gra  
 a(3,5)=a(1,5)*rho*gra 
 call FMIPWR(r,-n-1,a(3,6)); 	a(3,6)=-rho*a(3,6) 
!
!  Array A - 4th row 
 do i=1,6
 	a(4,i)=nul	
 enddo
!
!  Array A - 5th row 
 do i=1,6
 	a(5,i)=nul	
 enddo
 a(5,3)=a(3,3)/rho 
 a(5,6)=a(3,6)/rho  
!
!  Array A - 6th row  
 a(6,1)=3.0*ac*a(1,1)  
 a(6,2)=3.0*ac*a(1,2)
 call FMIPWR(r,+n-1,a(6,3)); 	a(6,3)=-(2.*fn+1.)*a(6,3)
 a(6,4)=3.0*ac*a(1,4)  
 a(6,5)=3.0*ac*a(1,5) 
 a(6,6)=nul
!
!
!  Array B - 1st, 2nd, 5th and 6th rows 
 do i=1,6
 	b(1,i)=nul       
 	b(2,i)=nul 
 	b(5,i)=nul       
 	b(6,i)=nul 	      
 enddo
! 
!  Array B - 3rd row
 call FMIPWR(r,+n  ,b(3,1)); 	b(3,1)=c1*b(3,1) 
 call FMIPWR(r,+n-2,b(3,2)); 	b(3,2)=c2*b(3,2) 
 b(3,3)=nul	  
 call FMIPWR(r,-n-1,b(3,4)); 	b(3,4)=c3*b(3,4) 
 call FMIPWR(r,-n-3,b(3,5)); 	b(3,5)=c4*b(3,5) 
 b(3,6)=nul
! 
!  Array B - 4rd row
 call FMIPWR(r,+n  ,b(4,1)); 	b(4,1)=d1*b(4,1) 
 call FMIPWR(r,+n-2,b(4,2)); 	b(4,2)=d2*b(4,2) 
 b(4,3)=nul	  
 call FMIPWR(r,-n-1,b(4,4)); 	b(4,4)=d3*b(4,4) 
 call FMIPWR(r,-n-3,b(4,5)); 	b(4,5)=d4*b(4,5) 
 b(4,6)=nul
!
!
do i=1, 6; do j=1, 6
	dire(i,j) = a(i,j) + mu*b(i,j) 
enddo; enddo
!
! 
end subroutine direct
!
!
!
!
!
SUBROUTINE SURF_BC (L, RAD, GGG, GRAV, BC)  
!
!-----------------------------------------------
! # Surface boundary conditions (tidal/loading)
!-----------------------------------------------
!
use FMZM 
use MAIN_MODULE 
implicit none  
!
integer  l
TYPE(FM) one, pi, grav, rad, ggg, bc(3)
!
    one=to_fm(1.0)  
    pi = 2.0*asin(one)	
!
    If   (l.ne.1) then 
         bc(1) = -grav*(2.*float(l)+1.)/4./pi/rad**2  
         bc(2) = to_fm(0.)  
         bc(3) = -4.*pi*ggg*(2.*float(l)+1.)/4./pi/rad**2 
         if(iload==0) bc(1) = to_fm(0.)
    Elseif(l.eq.1) then 
         bc(1) = -grav*(2.*float(l)+1.)/4./pi/rad**2  
         bc(2) = to_fm(0.)  
         bc(3) = to_fm(0.)    
         if(iload==0) bc(1) = to_fm(0.)      
    Endif
!
END SUBROUTINE surf_bc
!
!
!
!
!
        SUBROUTINE NORM 
!
!------------------------------------------------------
! Defines and normalizes the model parameters profiles 
!------------------------------------------------------
!
	USE MAIN_MODULE 
        IMPLICIT NONE   
!
	integer i, j, k, na, imark(8)  	
  	type(fm) one, pi 
	real*8 avero, avemu, avela, avevi
        character*12 rb(1:nla+2), rhob(1:nla+2), mub(1:nla+2), visb(1:nla+2)
 	character*12 vc(4)
	character*96 riga  
        character*20 date, timc	 
!
!
	one=to_fm(1.) 
	pi   = 2.0*asin(one)			! pi
	nt   = to_fm('0.667e-10')		! Newton's constant 
!
	write(109,*) ' ' 
	if(imode==1) write(109,*) ' Option #1: Uniform PREM layering' 
	if(imode==2) write(109,*) ' Option #2: User-supplied model' 
	if(imode==3) write(109,*) ' Option #3: User-supplied model' 		
!
!
!
!+++++++++++++++++++++++++++++
     if(imode==1) then 
!+++++++++++++++++++++++++++++
!
	    LT        = lth 		! thickness of the lithosphere	 
	    r(nla+2)  = ea     	        ! radius of the Earth 
	    r(nla+1)  = r(nla+2)-LT     ! radius of the litho-mantle boundary 	
	    r(1)      = ec              ! radius of the core mantle boundary 
!
! Thickness of each mantle layer (m) 
	 delta=(r(nla+1)-r(1))/float(nla) 
!
! Radii of the mantle layers (m) 	 
	r(0)=to_fm('0.0') 
         do i=1, NLA+1
	     r(i) = r(1) + (float(i)-1.0)*delta
	 enddo
!	 
! Radii of the mantle layers (km) 
	do i=0, NLA+2
		rk(i)=r(i)/1000. 
	enddo
!
!
! Uniform PREM-averaged core 
!
! -------- Core parameters (uniform core) 
!
! Average core density 
	 call prem_average(to_dp(rk(0)), to_dp(rk(1)), avero, avemu, avela) 
	 rho(0)=to_fm(avero)
! Average core Lame' Lambda
	 call prem_average(to_dp(rk(0)), to_dp(rk(1)), avero, avemu, avela) 
	 lam(0)=to_fm(avela)	 
! Core rigidity 
	 rig(0)=to_fm('0.0')
!
!
! Rigidity and density of all mantle layers (litho is i=NLA+1) 	 
	 do i=1, NLA+1
	     	call prem_average (to_dp(rk(i)),to_dp(rk(i+1)),avero, avemu, avela) 
	     	rho(i)=to_fm(avero)
	     	rig(i)=to_fm(avemu) 
		lam(i)=to_fm(avela)  
	 enddo
!
! Viscosity of the Core (uniform fluid core) 
        vis(0)=to_fm('0e0')
!
!	 
!
! Viscosity of each mantle layer 
         do i=1,NLA    
	 	call visco_average(to_dp(rk(i)),to_dp(rk(i+1)),avevi)
  	 	vis(i)=to_fm(avevi)	 	 
	 enddo
!
! Viscosity of the "Lithosphere" 
  	 vis(nla+1)=to_fm('1.0e35')
! The viscosity of the lithosphere is effectively = 'oo' 
!
! Plots of the profiles 
	 call  earth_par 
!
!++++++++++++++++
      endif 
!++++++++++++++++
!
!
!+++++++++++++++++++++++++++++++++++
     if(imode==2.or.imode==3) then 
!+++++++++++++++++++++++++++++++++++
!
! Reads the user-supplied profile- 
!
!
! Counts the rows and checks that their number is consistent with `nla' - see alma.inc  
if(imode==2)open(1,file=mod2,status='unknown')
if(imode==3)open(1,file=mod3,status='unknown')
!
do i=1,4; read(1,'(a96)') riga; enddo
k=0
do i=1,100000
read(1,'(a96)',end=102) riga; k=k+1
enddo
102 continue
if(k-2.ne.nla) then
!write(*,*) k-2, nla
	write(*,*) ' The number of rows declared is not consistent with the input file...'
	write(*,*) ' JOB ABORTED'
	stop 
        endif
  close(1)
!
!
if(imode==2)open(1,file=mod2,status='unknown')
if(imode==3)open(1,file=mod3,status='unknown')
do i=1,4; read(1,'(a96)') riga; enddo
!
	do 101 k=1,nla+2	
        read(1,'(a96)') riga 
!
! Reads line by line  
!
  na=0
  do i = 1,96
      if(riga(i:i).eq."'") then 
      na=na+1; imark(na)=i; endif
  enddo
!
  do i=1,4 
      vc(i)=riga(imark(2*i-1)+1:imark(2*i)-1)
  enddo
!write(*,'(4(a12,1x))') (vc(j),j=1,4)
!
  rb(k)   = vc(1); rhob(k) = vc(2)
  mub(k)  = vc(3); visb(k) = vc(4)	
!
!write(*,'(4(a12,1x))') rb(k), rhob(k), mub(k), visb(k) 
!		
101 continue 
!read(*,*) 	
!		
! Radii according to our conventions (m)
! 
	r(0)=to_fm('0.e0')
	do i=1, nla+2   
        	r(i)=to_fm(rb(nla+3-i))   
        enddo 
	LT=r(nla+2)-r(nla+1)
!read(*,*)
!
! Densities according to our conventions (kg/m^3)
	do i=0, nla+1
		rho(i)=to_fm(rhob(nla+2-i)) 
        enddo 
!
! Rigidity according to our conventions (Pa) 	
	do i=0, nla+1
		rig(i)=to_fm(mub(nla+2-i))    
        enddo 
!
! Viscosity according to our conventions (Pa.s)
	vis(0)=to_fm('0e0')
	do i=1, nla
		vis(i)=to_fm(visb(nla+2-i))   
        enddo
! The viscosity of the lithosphere is effectively = 'oo' 
        vis(nla+1)=to_fm('1e35')
!
! Plots of the profiles 
	 call earth_par 
!
!++++++++++++++++
      endif 
!++++++++++++++++
!
!
!
! Mass of the Earth  
	emass=rho(0)*r(1)*r(1)*r(1)  
	do i=1, NLA+1
		emass=emass+rho(i)*(r(i+1)*r(i+1)*r(i+1)-r(i)*r(i)*r(i))
	enddo 
 	emass=emass*4.0*pi/3.0
!
! Gravity at the boundaries 
	gra(0)=to_fm('0.0')
	 do i=1, nla+2
	 	gra(i)=to_fm('0.0')
		do k=0, i-1 	
		gra(i)=gra(i)+rho(k)*(r(k+1)**3-r(k)**3)
		enddo
		gra(i)=gra(i)*4.*pi*nt/3./r(i)/r(i)
	 enddo
!
	open(109,file=trim(adjustl(wdir))//'/ALMA/'//'alma-logfile.dat',status='unknown')	 
!
 	call date_and_time (date, timc)
!
	write(109,*) ' '	
        write(109,*) ' # Done by alma.f90 (SELEN port) on: ', & 
 				date(1:4),'.',date(5:6),'.',date(7:8),' time: ',&
		        	timc(1:2),'.',timc(3:4),'.',timc(5:6)
	write(109,*) ' '
	write(109,*) ' # Model : ', mod2	
	write(109,*) ' '
!
! Logs the first stage (dimensional variables) 
	call LOG_1
!
!
!'********************' 	
!    NORMALIZATION'
!'********************' 	
!
! pivots 
         ra0  = r(nla+2)                        ! raggio di riferimento
         rho0 = rho(1)                          ! densit' di riferimento
         rig0 = rig(1)                          ! rigidita' di riferimento
         t0 = to_fm(1000.0*365.250*24.0*3600.0) ! tempo di riferimento (1 kyr)
!		 
! NORMALIZED radii  
         do i=1, NLA+2
	     r(i) = r(i)/ra0
	 enddo
!
! NORMALIZED thickness 
	  delta=delta/ra0  
!	 
! NORMALIZED rigidities from bottom to top' 
         do i=0, NLA+1
	     rig(i) = rig(i)/rig0
	 enddo
!
! NORMALIZED "Lambdas" from bottom to top' 
         do i=0, NLA+1
	     lam(i) = lam(i)/rig0
	 enddo
!
! NORMALIZED densities from bottom to top' 
         do i=0, NLA+1
	     rho(i) = rho(i)/rho0
	 enddo
!
! NORMALIZED viscosity profile' 
	do i=0, nla+1
	 	vis(i)=vis(i)/t0/rig0
	enddo 
!
! NORMALIZED Newton constant' 
	nt=nt*rho0*rho0*ra0*ra0/rig0 
!
! NORMALIZED mass of the Earth (adim)' 
	emass=rho(0)*r(1)*r(1)*r(1)  
	do i=1, NLA+1
		emass=emass+rho(i)*(r(i+1)*r(i+1)*r(i+1)-r(i)*r(i)*r(i))
	enddo 
 	emass=emass*4.0*pi/3.0
!
! NORMALIZED gravity at the boundaries, bottom to top' 
	 do i=1, nla+2
	 	gra(i)=to_fm('0.0')
		do k=0, i-1 	
		gra(i)=gra(i)+rho(k)*(r(k+1)**3-r(k)**3)
		enddo
		gra(i)=gra(i)*4.*pi*nt/3./r(i)/r(i)
	 enddo
!
! NORMALIZED <<a>> constants  
         do i=0, NLA+1
	     ac(i) = (4.*pi/3.)*rho(i)*nt  
	 enddo
!	 
! Logs the second stage (a-dimensional variables) 
	 call LOG_2
!
! Closes the Log file
	 close(109)   
!
         return
         end 
!
!
!
!
! 
SUBROUTINE  EARTH_PAR 
!
!---------------------------------------------------------
! Profiles of the Earth properties for plotting purposes
!---------------------------------------------------------
!
USE FMZM
USE MAIN_MODULE 
implicit NONE
!
integer, parameter :: npu=5000   ! number of radial points, minus one 
integer i,j 
real*8 heaviside, hdiff 
real*8 densita, rigidita, lambda, viscosita
real*8 rho_prem, mu_prem, lambda_prem, visco
real*8 aalfa, bbeta, raggio  
!
!	 
! Radii of the mantle layers (km) 
do i=0, NLA+2
   rk(i)=to_dp(r(i))/1000.d0
enddo
!
aalfa=(rk(nla+2)-rk(0))/dfloat(npu)
bbeta=rk(0)-aalfa
!
open(20,file=trim(adjustl(wdir))//'/ALMA/'//'rho-r.dat',status='unknown') 
open(30,file=trim(adjustl(wdir))//'/ALMA/'//'rig-r.dat',status='unknown')
open(40,file=trim(adjustl(wdir))//'/ALMA/'//'vis-r.dat',status='unknown')
!open(5,file='lam-r.dat',status='unknown')
!
!
do i=1,npu+1
!
	raggio=aalfa*dfloat(i)+bbeta
!
	densita   =0.d0
	rigidita  =0.d0 
	lambda    =0.d0 
	viscosita =0.d0 
!
	do j=0,NLA+1
	hdiff    = heaviside(raggio-rk(j+1))-heaviside(raggio-rk(j)) 
	densita  = densita     -to_dp(rho(j))*hdiff
	rigidita = rigidita    -to_dp(rig(j))*hdiff
	lambda   = lambda      -to_dp(lam(j))*hdiff
	viscosita= viscosita   -to_dp(vis(j))*hdiff		
	enddo
!		
write(20,*) raggio, densita,        rho_prem(raggio)
write(30,*) raggio, rigidita/1e11,  mu_prem(raggio)/1e11
write(40,*) raggio, viscosita,      visco(raggio)
!write(5,*) raggio, lambda/1e11,    lambda_prem(raggio)/1e11
!
enddo 
!
 close(20); close(30); close(40)
!
END SUBROUTINE EARTH_PAR 
!
!
!
!
!
SUBROUTINE TIME_STEPS
!
!-----------------------------------------------------------
! Generates P+1 points between t=10^{M1} and t=10^{M2} kyrs 
! according to linearly or logarithmically spaced scales... 
! Revised GS May 8, 2009. 
!-----------------------------------------------------------
!
USE FMZM
USE MAIN_MODULE 
implicit NONE
!
integer jj, k
TYPE(fm) A, B, PP, MM1, MM2, DM1, DM2 
!
! CALL FMI2M(M1,MM1)
! CALL FMI2M(M2,MM2)
 
!if(m1==0.)m1=1.e-6
 
 MM2=M2
 MM1=M1 
 	
 CALL FMI2M(P,PP)
!
!
 If(tscale=='log')     then 
!
 B=(MM2-MM1)/PP ; A=MM1-B 
!
 do jj=1, p+1 
     call FMPWR(to_fm(10.), A+float(jj)*B, time(jj))
 enddo 
!
 elseif(tscale=='lin') then 
!
	k=0
	do jj=0, int(m2) 
		k=k+1
		if(jj==0) then 
			   time(k)=1.E-6 
			  else 
			   time(k)=float(jj) 
		endif
	enddo
!		  
 endif

!
END SUBROUTINE TIME_STEPS
!
!
!
! 
!
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
!
!++++++++++++++++++++++++++++++++++      
! Adapted from Numerical recipes
!++++++++++++++++++++++++++++++++++
!
      USE FMZM
      IMPLICIT NONE 
!
      INTEGER NMAX 
      PARAMETER (NMAX=500)   
!      
      INTEGER N,NP,INDX(N)
      INTEGER I,J,K, IMAX
!
      TYPE(FM) D,A(NP,NP),TINNY
      TYPE(FM) AAMAX,DUM,SUM,VV(NMAX)
!       
      tinny=to_fm(1.e-20)      
      d=to_fm(1.)
!
      do 12 i=1,n
        aamax=0.0
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.0) pause 'singular matrix in ludcmp'
        vv(i)=1.0/aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.0
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINNY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software W^3.
!
!
!
!
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
!
!++++++++++++++++++++++++++++++++++      
! Adapted from Numerical recipes
!++++++++++++++++++++++++++++++++++
!
      USE FMZM
      IMPLICIT NONE 
! 
      INTEGER n,np,indx(n)
      TYPE(FM) a(np,np),b(n)
      INTEGER i,ii,j,ll
      TYPE(FM) sum
!
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software W^3.
!
!
!
!
!
  SUBROUTINE PREM_AVERAGE(R1, R2, AVERO, AVEMU, AVELA) 
!  
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Computes the volume-averaged parameters in the range [r1:r2] 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
  implicit none	
!
  real*8 r1, r2, avero, avemu, avela 
  real*8 rho, mu, lambda
  real*8 rhor2, mur2, lambdar2 
  external rhor2, mur2, lambdar2 
!
  if(r1.gt.r2) then 
  write(*,*) 'Error in sbr prem_average'		   
  stop
  endif	 
!
  call qtrap (rhor2,r1,r2,avero)
  avero=3.d0*avero/(r2**3-r1**3)
!	
  call qtrap (mur2,r1,r2,avemu)
  avemu=3.d0*avemu/(r2**3-r1**3)	
!	
  call qtrap (lambdar2,r1,r2,avela)
  avela=3.d0*avela/(r2**3-r1**3)
!		
  END SUBROUTINE PREM_AVERAGE 
!
!
!
!
!
	SUBROUTINE VISCO_AVERAGE(R1, R2, AVEVI) 
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Computes the volume-averaged viscosity in the range [r1:r2] 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
	implicit none	
!
	real*8 r1, r2, avevi
	real*8 visco
	real*8 viscor2
	external viscor2
!
	if(r1.gt.r2) then 
		     write(*,*) 'Error in sbr visco_average '
		     stop
		     endif	 
!
	call qtrap (viscor2,r1,r2,avevi)
	avevi=3.d0*avevi/(r2**3-r1**3)
!		
	END SUBROUTINE VISCO_AVERAGE 
!
!
!
!
!
      SUBROUTINE QTRAP(FUNC,A,B,S)
! Adapted from Numerical Recipes 
      INTEGER JMAX
      REAL*8 a,b,func,s,EPS
      EXTERNAL func
      PARAMETER (EPS=1.d-3, JMAX=20)
      INTEGER j
      REAL*8 olds
      olds=-1.d30
      do 11 j=1,JMAX
        call trapzd(func,a,b,s,j)
        if (abs(s-olds).lt.EPS*abs(olds)) return
        olds=s
11    continue
      pause 'too many steps in qtrap'
      END
!      
!
!
!
!     
      SUBROUTINE TRAPZD(FUNC,A,B,S,N)
! Adapted from Numerical Recipes 
      INTEGER n
      REAL*8 a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL*8 del,sum,tnm,x
      if (n.eq.1) then
        s=0.5d0*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5d0*del
        sum=0.d0
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5d0*(s+(b-a)*sum/tnm)
      endif
      return
      END
!
!
!
!
!
	FUNCTION VISCOR2(R) 
! Viscosity * r^2
	implicit NONE 
		real*8 r, visco, viscor2 
		viscor2=r**2*visco(r)  
	end function viscor2 
!
!
!
!
!
	FUNCTION MUR2(R) 
! Rigidity * r^2
	implicit NONE 
		real*8 r, mu_prem, mur2 
		mur2=r**2*mu_prem(r)  
	end function mur2 
!
!
!
!
!
	FUNCTION LAMBDAR2(R) 
! 	2nd Lame constant * r^2
	implicit NONE  
		real*8 r, lambda_prem, lambdar2 
		lambdar2=r**2*lambda_prem(r)  
	end function lambdar2 
!
!
!
!
!
	FUNCTION RHOR2(R) 
!	 Density * r^2 
	implicit NONE 
		real*8 r, rho_prem, rhor2 
		rhor2=r**2*rho_prem(r)  
	end function rhor2 
!
!
!
!
!
	FUNCTION MU_PREM(R)
! 	Rigidity (Pa) according to PREM  
	implicit NONE 
	real*8 r, mu_prem, rho_prem, vs_prem 
		mu_prem=rho_prem(r)*vs_prem(r)**2  
	end function mu_prem
!
!
!
!
!
	FUNCTION LAMBDA_PREM(R) 
! 	2nd Lame constant (Pa) according to PREM  
	implicit NONE 
	real*8 r, lambda_prem, rho_prem, vp_prem, vs_prem  	
		LAMBDA_prem=rho_prem(r)*(vp_prem(r)**2-2.*vs_prem(r)**2)  
	end function lambda_prem
!
!
!
!
!
	FUNCTION RHO_PREM (RAD) 
!	
!--------------------------------------------------------
! Returns the PREM density (kg/m^3) at radius 'rad' (km)  
!--------------------------------------------------------
!
	USE MAIN_MODULE 
	implicit NONE 
!
	real*8, parameter :: erad = ea/1.e3 
	real*8  z, rad, rho_prem
!	
	z=rad/erad 
!
! Outside the Earth  
	if(rad.lt.0.d0.or.rad.gt.6371.d0) rho_prem = 0.d0           	        
	
! Inner core  
	if(0.0d0.le.rad.and.rad.lt.1221.5d0) & 
	rho_prem = 13.0885d0 - 8.8381d0*(z)**2  
!                      
! Outer core  					 
	if(1221.5d0.le.rad.and.rad.lt.3480.0d0) & 
	rho_prem = 12.5815d0 -1.2638d0*(z) -3.6426d0*(z)**2  -5.5281d0*(z)**3  
!  
! Lower mantle  
	if(3480.0d0.le.rad.and.rad.lt.5701.0d0 ) & 
	rho_prem =  7.9565d0 -6.4761d0*(z) +5.5283d0*(z)**2  -3.0807d0*(z)**3 
!   	
! Transition zone (I)  
	if(5701.0d0.le.rad.and.rad.lt.5771.0d0 ) & 
	rho_prem = 5.3197d0  -1.4836d0*(z)
!             
! Transition zone (II)  
	if(5771.0d0.le.rad.and.rad.lt.5971.0d0 ) & 
	rho_prem = 11.2494d0 -8.0298d0*(z)  
!          					
! Transition zone (III) 
	if(5971.0d0.le.rad.and.rad.lt.6151.0d0) &
	rho_prem = 7.1089d0  -3.8045d0*(z)          									
!
! LVZ & LID  
	if(6151.0d0.le.rad.and.rad.lt.6346.6d0) & 
	rho_prem = 2.6910d0  +0.6924d0*(z)             				
!
! Lower Crust    
	if(6346.6d0.le.rad.and.rad.lt.6356.0d0) rho_prem = 2.900d0                  
!	
! Continental Upper Crust  
	if(ioc.eq.0.and.6356.0d0 .le.rad.and.rad.le.6371.0d0) rho_prem = 2.900d0   
!
! Oceanic Upper Crust  
	if(ioc.eq.1.and.6356.0d0 .le.rad.and.rad.le.6368.0d0) rho_prem = 2.900d0    
!
! Ocean layer 
	if(ioc.eq.1.and.6368.0d0 .le.rad.and.rad.le.6371.0d0) rho_prem = 1.020d0    
!
!	
	RHO_PREM=RHO_PREM*1000.0d0     ! Density is now in units of kg/m**3 
!					
	END FUNCTION RHO_PREM 
!
!
!
!
!
	FUNCTION VS_PREM (RAD) 
!
!---------------------------------------------------------
! Returns the S waves velocity (m/s) at radius 'rad' (km)  
!---------------------------------------------------------
!
	USE MAIN_MODULE 
	implicit NONE 
!
	real*8, parameter :: erad = ea/1.e3 
	real*8 z, rad, vs_prem  	
!
	z=rad/erad 
!
! Outside the Earth
	if(rad.lt.0.0d0.or.rad.gt.6371.0d0) vs_prem = 0.0d0           	
!	
! Inner core 
	if(0.0d0.le.rad.and.rad.lt.1221.5d0)  & 
	vs_prem = 3.6678d0  - 4.4475d0*(z)**2        
!	
! Outer core 						 
	if(1221.5d0.le.rad.and.rad.lt.3480.0d0) vs_prem = 0.0d0	               

! Lower mantle (I) 
	if(3480.0d0.le.rad.and.rad.lt.3630.0d0) & 
	vs_prem = 6.9254d0  +1.4672d0*(z) -2.0834d0*(z)**2  +0.9783d0*(z)**3    	
!
! Lower mantle (II) 
	if(3630.0d0.le.rad.and.rad.lt.5600.0d0) & 
	vs_prem = 11.1671d0 -13.7818d0*(z) +17.4575d0*(z)**2 -9.2777d0*(z)**3    	
!
! Lower mantle (III) 
	if(5600.0d0.le.rad.and.rad.lt.5701.0d0) & 
	vs_prem = 22.3459d0 -17.2473d0*(z) -2.0834d0*(z)**2 +0.9783d0*(z)**2   	     
!
! Transition zone (I)   
	if(5701.0d0.le.rad.and.rad.lt.5771.0d0) & 
	vs_prem = 9.9839d0 -4.9324d0*(z)             
!
! Transition zone (II)   
	if(5771.0d0.le.rad.and.rad.lt.5971.0d0) & 
	vs_prem = 22.3512d0 -18.5856d0*(z)           							
!
! Transition zone (III)   
	if(5971.0d0.le.rad.and.rad.lt.6151.0d0) & 
	vs_prem = 8.9496d0 -4.4597d0*(z)           									
!
 ! LVZ & LID  
	if(6151.0d0.le.rad.and.rad.lt.6346.6d0) & 
	vs_prem = 2.1519d0 +2.3481d0*(z)             				
!
! Lower Crust    
	if(6346.6d0.le.rad.and.rad.lt.6356.0d0) vs_prem = 3.900d0                
!	
! Continental Upper Crust  
	if(ioc.eq.0.and.6356.0d0.le.rad.and.rad.le.6371.0d0) vs_prem = 3.200d0    
!
! Oceanic Upper Crust  
	if(ioc.eq.1.and.6356.0d0.le.rad.and.rad.le.6368.0d0) vs_prem = 3.200d0     
!
! Ocean layer 
	if(ioc.eq.1.and.6368.0d0.le.rad.and.rad.le.6371.0d0) vs_prem = 0.0d0     
!			
	vs_prem = vs_prem*1000.d0     ! vs_prem is now in m/s 
!		
	END FUNCTION vs_prem
!
!
!
!
!
	FUNCTION VP_PREM (RAD) 
!
!---------------------------------------------------------
! Returns the P waves velocity (m/s) at radius 'rad' (km)  
!---------------------------------------------------------
!
	USE MAIN_MODULE 
	implicit NONE 
!
	real*8, parameter :: erad = ea/1.e3 
	real*8  z, rad, vp_prem  
!	
	z=rad/erad 	
!	
! Outside the Earth
	if(rad.lt.0.0d0.or.rad.gt.6371.0d0) vp_prem = 0.0d0           	        
	
! Inner core 
	if(0.0d0.le.rad.and.rad.lt.1221.5d0) & 
	vp_prem = 11.2622d0 - 6.3640d0*(z)**2        
!				
! Outer core 			 
	if(1221.5d0.le.rad.and.rad.lt.3480.0d0) & 
	vp_prem = 11.0487d0 -4.0362d0*(z) +4.8023d0*(z)**2 -13.5732d0*(z)**3 
!
! Lower mantle (I) 
	if(3480.0d0.le.rad.and.rad.lt.3630.d0) & 
	vp_prem  = 15.3891d0 -5.3181d0*(z) +5.5242d0*(z)**2 -2.5114d0*(z)**3    	
!
! Lower mantle (II) 
	if(3630.0d0.le.rad.and.rad.lt.5600.0d0) & 
	vp_prem  = 24.9520d0 -40.4673d0*(z) +51.4823d0*(z)**2 -26.6419d0*(z)**3    	
!
! Lower mantle (III) 
	if(5600.0d0.le.rad.and.rad.lt.5701.0d0) & 
	vp_prem  = 29.2766d0 -23.6027d0*(z) +5.5242d0*(z)**2 -2.5514d0*(z)**3  
!	        
 ! Transition zone (I)   
	if(5701.0d0.le.rad.and.rad.lt.5771.0d0) & 
	vp_prem = 19.0957d0 -9.8672d0*(z)             
!
! Transition zone (II) 
	if(5771.0d0.le.rad.and.rad.lt.5971.0d0) & 
	vp_prem = 39.7027d0 -32.6166d0*(z)           							
!
! Transition zone (III)  
	if(5971.0d0.le.rad.and.rad.lt.6151.0d0) & 
	vp_prem = 20.3926d0  -12.2596d0*(z)           									
!
! LVZ & LID  
	if(6151.0d0.le.rad.and.rad.lt.6346.6d0) & 
	vp_prem = 4.1875d0 +3.9382d0*(z)             				
!
! Lower Crust    
	if(6346.6d0.le.rad.and.rad.lt.6356.0d0) vp_prem = 6.800d0               
!	
! Continental Upper Crust  
	if(ioc.eq.0.and.6356.0d0.le.rad.and.rad.le.6371.0d0) vp_prem = 5.800d0     
!
! Oceanic Upper Crust  
	if(ioc.eq.1.and.6356.0d0.le.rad.and.rad.le.6368.0d0) vp_prem = 5.800d0   
!
! Ocean layer 
	if(ioc.eq.1.and.6368.0d0.le.rad.and.rad.le.6371.0d0) vp_prem = 1.450d0    
!			
	vp_prem=vp_prem*1000.d0     ! vp_prem is now in m/s 
!		
	END FUNCTION vp_prem
!
!
!
!
!
FUNCTION HEAVISIDE (X)
!
! Heaviside function 
!

   USE MAIN_MODULE 

  implicit none
  TYPE(FM) x
  real*8 xx
  real*8 heaviside 
!
  xx=to_dp(x)
!
  heaviside=0.d0 
  if(x.ge.0.d0) heaviside=1.d0 
!
END FUNCTION HEAVISIDE
!
!
!
!
!
!
FUNCTION VISCO (RAD) 
!
!++++++++++++++++++++++++++++++++++++++++++++++ 
! Defines the viscosity profile for imode = 1 
!++++++++++++++++++++++++++++++++++++++++++++++
!
USE MAIN_MODULE
implicit none 
real*8 rad, z, visco 
! 
! General settings 
real*8, parameter :: ERAD=6371.d0     ! Radius of the Earth, km 
real*8, parameter :: EAK=EA/1d3       ! Radius of the Earth, km 
real*8, parameter :: LITK=LTH/1d3      ! Thickness of the lithosphere, km  
real*8, parameter :: ECK=EC/1d3       ! Radius of the CMB, km 
real*8, parameter :: ELK=EAK-LITK     ! Radius of the base of the lithosphere, km 
!
!
! -----------------------
!  Settings for kv=1 
! -----------------------
!
! Radius of the viscosity minimum in the asthenosphere, km   
real*8, parameter :: ELVZK=6000d0  
!
! Viscosity just at the base of the lithosphere, Pa.s 
real*8, parameter :: ETA1=1d21
!
! Viscosity at the CMB and at radius ELVZK, Pa.s
real*8, parameter :: ETA0=1d20
!
!
! -----------------------  
!  Settings for kv=2 
! ----------------------- 
!
! A scaling factor (~ the viscosity at the base of the lithosphere), Pa.s
real*8, parameter :: ETA_SCALE=1d21
!
!
! -----------------------  
!  Settings for kv=3 
! ----------------------- 
!
! Thickness of the LVZ, km 
 real*8, parameter :: LVZTK=200d0
!
! Viscosity of the LVZ, Pa.s 
 real*8, parameter :: VISCO_LVZ=1d19
!
!
! -----------------------  
!  Settings for kv=4 
! ----------------------- 
!
! Radii of the interfaces "670" and "420"
 real*8, parameter :: E670K=EAK-670d0
 real*8, parameter :: E420K=EAK-420d0 
!
! Viscosities of the "lower mantle, "TZ", and "Shallow upper mantle"
 real*8, parameter :: visco_lm=2d21
 real*8, parameter :: visco_tz=1d22
 real*8, parameter :: visco_um=1d21  
!
!
!
!
! **************************************
! Definition of the viscosity profiles
! **************************************
  visco=0d0
!
!
!
!
! kv=1
! ========================================= 
!    ''Continuous viscosity profile'' 
! Spada & Boschi, GJI 2006, Figure 9(c)
! ========================================= 
!
if(kv==1) then 
!
! Defines the viscosity between the CMB and the base of the litho 
if(ECK.le.rad.and.rad.le.ELK) &
       visco = eta0 + (eta1-eta0)*((rad-ELVZK)/(ELK-ELVZK))**2*((rad-ECK)/(ELK-ECK))			
endif 			
!
!
!
! kv=2
! =======================================
!      ''Convex viscosity profile'' 
!       L. Hanik - Thesis, page 59. 
! =======================================
!
if(kv==2) then 
!
! Normalized depth 
z=(EAK-rad)/(ELK-ECK)
! 
! Defines the viscosity between the CMB and the base of the litho 
if(ECK.le.rad.and.rad.le.ELK) & 
       visco=eta_scale*(1d0 + 214.3d0*z*dexp(-16.07d0*(0.7d0-z)**2))   
endif 	
!
!
!
! kv=3
! ====================================================
!       ''Convex viscosity profile with a LVZ'' 
!            L. Hanik - Thesis, page 90. 
! ====================================================
!
if(kv==3) then 
!
! Within the LVZ
if(EAK-LITK-LVZTK.le.rad.and.rad.le.ELK) visco=visco_lvz
!
! Normalized depth 
z=(EAK-rad)/(ELK-ECK)
!
! Between the CMB and the base of the LVZ 
if(ECK.le.rad.and.rad.lt.EAK-LITK-LVZTK) & 
       visco=1d21*(1d0 + 214.3d0*z*dexp(-16.07d0*(0.7d0-z)**2))   
endif 	
!
!
!
! kv=4
! =================================================
!        ''Three-layers viscosity profile'' 
!  Lower mantle, transition zone, & upper mantle 
! =================================================
!
if(kv==4) then 
!
! "Lower mantle" viscosity 
if      (ECK.le.rad.and.rad.le.E670K) visco=visco_lm
!
! "Transition zone viscosity"
if    (E670K.lt.rad.and.rad.le.E420K) visco=visco_tz 
!
! "Upper mantle viscosity"
if (E420K.lt.rad.and.rad.le.EAK-LITK) visco=visco_um 
!
endif
!
END FUNCTION VISCO 
!
!
!
!
!
	SUBROUTINE LOG_1
!
!------------------------------------------------------------
! Prints the un-normalized model parameters on 'alma-logfile.dat'
!------------------------------------------------------------
!
	USE MAIN_MODULE 
        IMPLICIT NONE
	integer i, onefound, densiv 
	type(fm)average_density, pi 
!
if(imode==1)write(109,*) ' Thickness of each mantle layer (km) = ', to_sp(delta)/1000.0       
	   write(109,*) ' Litho thickness (km) = ', to_sp(LT)/1000.0 
	   write(109,*) ' Earth radius (km) = ', to_sp(r(nla+2))/1000.0 
	   write(109,*) ' CMB radius (km) = ',   to_sp(r(1))/1000.0 
           write(109,*) ' Number of mantle layers = ', nla
           write(109,*) ' Total number of  layers = ', nla+2	
!	 
write(109,*) ' '
write(109,*) ' Radii of the interfaces from bottom to top, (km)'
            i=0; write(109,*) ' ', i, ':', to_sp(rk(i)), '(Earth center)'
            i=1; write(109,*) ' ', i, ':', to_sp(rk(i)), '(CMB)'	
	do i=2, NLA+1
        	 write(109,*) ' ', i, ':', to_sp(rk(i))
	enddo
        i=NLA+2; write(109,*) ' ', i, ':', to_sp(rk(i)), '(External surface)' 
!
write(109,*) ' '
write(109,*) ' Rigidity of the layers from bottom to top, (*10^11 Pa)' 
             i=0; write(109,*) ' ', i, ':', to_sp(rig(i))/1.e11, '(Core)'
	 do i=1, NLA
	          write(109,*) ' ', i, ':', to_sp(rig(i))/1.e11   
	 enddo
	 i=NLA+1; write(109,*) ' ', i, ':', to_sp(rig(i))/1.e11, '(Lithosphere)'	
!
!     
!if(compr==1) then 
!write(109,*) ' '
!write(109,*) ' "Lambda" of the layers from bottom to top, (*10^11 Pa)' 
!             i=0; write(109,*) ' ', i, ':', to_sp(lam(i))/1.e11, '(Core)'
!	 do i=1, NLA
!	          write(109,*) ' ', i, ':', to_sp(lam(i))/1.e11   
!	 enddo
!	 i=NLA+1; write(109,*) ' ', i, ':', to_sp(lam(i))/1.e11, '(Lithosphere)'	
!endif
!
!
write(109,*) ' '
write(109,*) ' Density of the layers from bottom to top, (kg/m^3)' 
             i=0; write(109,*) ' ', i, ':', to_sp(rho(i)), '(Core)'
	 do i=1, NLA
	          write(109,*) ' ', i, ':', to_sp(rho(i))    
	 enddo	
	 i=NLA+1; write(109,*) ' ', i, ':', to_sp(rho(i)), '(Lithosphere)'
!
write(109,*) ' '
	densiv=0
	onefound=0
	do i=0,nla
		if(rho(i).lt.rho(i+1)) then 
		onefound=1 
		densiv=1
		write(109,*) ' Density inversion for layers ', i, i+1
		densiv=0
		endif		
	enddo	 
!
	if(onefound==0) write(109,*) ' NO density inversions found!'
!
!
write(109,*) ' '
write(109,*) ' Viscosity of the layers from bottom to top, (Pa.s)'	
             i=0; write(109,*) ' ', i, ':', to_dp(vis(i)), '(Core)'
	 do i=1, NLA
	          write(109,*) ' ', i, ':', to_dp(vis(i))    
	 enddo	
	 i=NLA+1; write(109,*) ' ', i, ':', ' Infinity', ' (Lithosphere)'
!
write(109,*) ' '
write(109,*) ' The *density*   profile is reported on file "rho-r.dat" '
write(109,*) ' The *rigidity*  profile is reported on file "rig-r.dat" '
write(109,*) ' The *viscosity* profile is reported on file "vis-r.dat" '  
!
write(109,*) ' '
write(109,*) ' Mass of the Earth (kg) = ', to_dp(emass) 
!
!
pi = 2.0*asin(1.0)
average_density= emass/((4./3.)*pi*r(nla+2)**3) 
write(109,*) ' '
write(109,'(A5,A40,E14.8)') '>>>>>', 'Average density of the Earth model is: ', & 
			  	      to_dp(average_density) 
!
write(109,*) ' '
write(109,*) ' Gravity at the interfaces from bottom to top, (m/s^2)'
            i=0; write(109,*) ' ', i, ':', to_sp(rk(i)), to_sp(gra(i)), '(Earth center)'
            i=1; write(109,*) ' ', i, ':', to_sp(rk(i)), to_sp(gra(i)), '(CMB)'	
            do i=2, NLA+1
	         write(109,*) ' ', i, ':', to_sp(rk(i)), to_sp(gra(i))
            enddo
            i=NLA+2; write(109,*) ' ', i, ':', to_sp(rk(i)), to_sp(gra(i)), '(External surface)'	 
!
end subroutine LOG_1 
!
!
!
!
!
	SUBROUTINE LOG_2
!
!---------------------------------------------------------
! Prints the normalized model parameters on 'alma-logfile.dat'
!---------------------------------------------------------
!
	USE MAIN_MODULE 
        IMPLICIT NONE
	integer i 	
!
write(109,*) ' '
	write(109,*) ' Reference radius (m) =', to_sp(ra0) 
	write(109,*) ' Reference density (kg/m^3) =', to_sp(rho0) 
	write(109,*) ' Reference rigidity (Pa) =', to_sp(rig0) 
	write(109,*) ' Reference time (s) =', to_sp(t0) 
!
!
write(109,*) ' '
write(109,*) ' NORMALIZED Radii of the interfaces from bottom to top'
            i=0; write(109,*) ' ', i, ':', to_sp(r(i)), '(Earth center)'
            i=1; write(109,*) ' ', i, ':', to_sp(r(i)), '(CMB)'	
	do i=2, NLA+1
        	 write(109,*) ' ', i, ':', to_sp(r(i))
	enddo
        i=NLA+2; write(109,*) ' ', i, ':', to_sp(r(i)), '(External surface)' 
!
write(109,*) ' '
write(109,*) ' NORMALIZED Rigidity of the layers from bottom to top' 
             i=0; write(109,*) ' ', i, ':', to_sp(rig(i)), '(Core)'
	 do i=1, NLA
	          write(109,*) ' ', i, ':', to_sp(rig(i))  
	 enddo
	 i=NLA+1; write(109,*) ' ', i, ':', to_sp(rig(i)), '(Lithosphere)'	
!
!if(compr==1)then
!write(109,*) ' '
!write(109,*) ' NORMALIZED "Lambda" of the layers from bottom to top' 
!             i=0; write(109,*) ' ', i, ':', to_sp(lam(i)), '(Core)'
!	 do i=1, NLA
!	          write(109,*) ' ', i, ':', to_sp(lam(i))  
!	 enddo
!	 i=NLA+1; write(109,*) ' ', i, ':', to_sp(lam(i)), '(Lithosphere)'
!endif	
!
!
write(109,*) ' '
write(109,*) ' NORMALIZED Density of the layers from bottom to top' 
             i=0; write(109,*) ' ', i, ':', to_sp(rho(i)), '(Core)'
	 do i=1, NLA
	          write(109,*) ' ', i, ':', to_sp(rho(i))    
	 enddo	
	 i=NLA+1; write(109,*) ' ', i, ':', to_sp(rho(i)), '(Lithosphere)'
!
write(109,*) ' '
write(109,*) ' NORMALIZED Viscosity of the layers from bottom to top, (Pa.s)'	
             i=0; write(109,*) ' ', i, ':', to_dp(vis(i)), '(Core)'
	 do i=1, NLA
	          write(109,*) ' ', i, ':', to_dp(vis(i))    
	 enddo	
	 i=NLA+1; write(109,*) ' ', i, ':', ' Infinity ', '(Lithosphere)'
!
	 write(109,*) ' '
	 write(109,*) ' NORMALIZED Mass of the Earth = ', to_dp(emass) 
	 write(109,*) ' '
	 write(109,*) ' NORMALIZED Newton constant = ', to_dp(nt) 
!
write(109,*) ' '
write(109,*) ' NORMALIZED Gravity at the interfaces from bottom to top'
            i=0; write(109,*) ' ', i, ':', to_sp(rk(i)), to_sp(gra(i)), '(Earth center)'
            i=1; write(109,*) ' ', i, ':', to_sp(rk(i)), to_sp(gra(i)), '(CMB)'	
            do i=2, NLA+1
	         write(109,*) ' ', i, ':', to_sp(rk(i)), to_sp(gra(i))
            enddo
            i=NLA+2; write(109,*) ' ', i, ':', to_sp(rk(i)), to_sp(gra(i)), '(External surface)'	 
!
!write(109,*) ' '
!write(109,*) ' NORMALIZED "a" constants from bottom to top' 
!             i=0; write(109,*) ' ', i, ':', to_sp(ac(i)), '(Core)'
!	 do i=1, NLA
!	          write(109,*) ' ', i, ':', to_sp(ac(i))    
!	 enddo	
!	 i=NLA+1; write(109,*) ' ', i, ':', to_sp(ac(i)), '(Lithosphere)'
!
end subroutine LOG_2 
!
!
!
!
!
!
!
!
SUBROUTINE CMB (N, R, RHO, AC, CMB_A)
!
!------------------------------------
! # Boundary conditions at the CMB
!------------------------------------
!
USE FMZM
implicit NONE 
integer  n
type(FM) r, rho, ac, fn, cmb_a(6,3)    
!
!
fn=float(n) 
!
 call FMIPWR(r,+n-1,cmb_a(1,1));cmb_a(1,1)=-cmb_a(1,1)/ac  	 
 				cmb_a(1,2)=to_fm(0.0) 
				cmb_a(1,3)=to_fm(1.0) 	 
! 
				cmb_a(2,1)=to_fm(0.0) 
				cmb_a(2,2)=to_fm(1.0)  
				cmb_a(2,3)=to_fm(0.0)
!				 
				cmb_a(3,1)=to_fm(0.0)
				cmb_a(3,2)=to_fm(0.0)
				cmb_a(3,3)=rho*ac*r   
!				 
				cmb_a(4,1)=to_fm(0.0)
				cmb_a(4,2)=to_fm(0.0)
				cmb_a(4,3)=to_fm(0.0)
!
 call FMIPWR(r,+n  ,cmb_a(5,1));cmb_a(5,1)=cmb_a(5,1)  
				cmb_a(5,2)=to_fm(0.0)
				cmb_a(5,3)=to_fm(0.0)
!
 call FMIPWR(r,+n-1,cmb_a(6,1));cmb_a(6,1)=cmb_a(6,1)*2.*(fn-1.)	 
   				cmb_a(6,2)=to_fm(0.0)
				cmb_a(6,3)=3.*ac 		 
!
end subroutine CMB 
!
