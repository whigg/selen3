!
! originally written in ??/??/2009 - Retouched for the sea level benchmark 
! by GS on august 2010. 
!
!
!
 INCLUDE "harmonics.f90"
 PROGRAM R
 IMPLICIT NONE
 INCLUDE "data.inc"
 INTEGER I, J, N, K, L, DOM
!
!
! Input binary files... 
 CHARACTER*20, PARAMETER :: SH_FILE = "sh-r44-l128_bc00.bin"
 CHARACTER*20, PARAMETER :: SH_S_FILE = "shs_bc00.bin"
 CHARACTER*20, PARAMETER :: SH_U_FILE = "shu_bc00.bin"
 CHARACTER*20, PARAMETER :: SH_N_FILE = "shn_bc00.bin"
 CHARACTER*20, PARAMETER :: SH_V_FILE = "shv_bc00.bin"
!
! Output files for S, U, N, and V
 CHARACTER*20, PARAMETER :: OUT_S_FILE = "s.dat"
 CHARACTER*20, PARAMETER :: OUT_U_FILE = "u.dat"
 CHARACTER*20, PARAMETER :: OUT_N_FILE = "n.dat"
 CHARACTER*20, PARAMETER :: OUT_V_FILE = "v.dat"
!
 INTEGER, PARAMETER :: NG=401 
 REAL*4, PARAMETER :: ALFA=180./FLOAT(NG-1), BETA=-90.-ALFA 


 COMPLEX*16 ARMOY (JMAX),   & 
            DARMOY(JMAX),   &
            CS(JMAX,0:NN), & 
	    CU(JMAX,0:NN), & 
	    CN(JMAX,0:NN), &
	    CV(JMAX,0:NN)
!	    
 REAL*8 REC_S(0:NN,NG), REC_U(0:NN,NG), REC_N(0:NN,NG) 
 REAL*8 LAT(NG) 
!
!
!
!
 open (11,file=SH_S_FILE,status=un,form='UNFORMATTED') 
 open (12,file=SH_U_FILE,status=un,form='UNFORMATTED') 
 open (13,file=SH_N_FILE,status=un,form='UNFORMATTED') 
 open (14,file=SH_V_FILE,status=un,form='UNFORMATTED')
!
!
!
 open (101,file=SH_S_FILE,status=un,form='UNFORMATTED') 
 open (102,file=SH_U_FILE,status=un,form='UNFORMATTED') 
 open (103,file=SH_N_FILE,status=un,form='UNFORMATTED') 
 open (104,file=SH_V_FILE,status=un,form='UNFORMATTED')
!
 read (101) CS
 read (102) CU
 read (103) CN
 read (104) CV 
! 
 close(101) ; close(102) ; close(103) ; close(104) 

! --- Harmonics 

	do 10 k=0, NN 
	write(*,*) "time: ",k 


	do 10 i=1, NG 
!	
	lat(i)=alfa*float(i)+beta 
!	
		if(lat(i)>+90.) lat(i)=+90.  
		if(lat(i)<-90.) lat(i)=-90.  
!
	REC_S(k,i)=0d0
	REC_U(k,i)=0d0 
	REC_N(k,i)=0d0 
	REC_U(k,i)=0d0 
!	
		call HARMO           (LMAX, 90d0, LAT(i),  ARMOY)
                call GRAD_THETA_HARMO(LMAX, 90d0, LAT(i), DARMOY)
!		
        do 10 j=1, jmax 
!
			REC_S(k,i)=REC_S(k,i)+(2.-DOM(J))*DBLE(CS(j,k)* ARMOY(j))
			REC_U(k,i)=REC_U(k,i)+(2.-DOM(J))*DBLE(CU(j,k)* ARMOY(j))
			REC_N(k,i)=REC_N(k,i)+(2.-DOM(J))*DBLE(CN(j,k)* ARMOY(j))
			REC_V(k,i)=REC_V(k,i)+(2.-DOM(J))*DBLE(CV(j,k)*DARMOY(j))
!		
10	continue
!
	do i=1, NG 
!	
		write(11,'(f14.6,30(1x,f14.6))') 90.-lat(i), (rec_s(k,i),k=0,nn)
		write(12,'(f14.6,30(1x,f14.6))') 90.-lat(i), (rec_u(k,i),k=0,nn)
		write(13,'(f14.6,30(1x,f14.6))') 90.-lat(i), (rec_n(k,i),k=0,nn)
		write(14,'(f14.6,30(1x,f14.6))') 90.-lat(i), (rec_v(k,i),k=0,nn)
!
	enddo 
!
	close(14) ; close(13) ; close(12) ; close(11) 
!
	stop
	end
!
!
!
