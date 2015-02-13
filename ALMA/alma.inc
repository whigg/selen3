 ! - - - - - - - - - - - - - - - - - - - - - - - -
 ! This is file 'alma.inc' - SELEN port July 2009 
 ! - - - - - - - - - - - - - - - - - - - - - - - -
 USE FMZM
 integer, parameter :: imode=2   !  Mode 2: Multi-layer user-supplied model
 ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
 ! # Parameters for imode=2  
 ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
 character*100, parameter :: mod2 = '/Users/giorgiospada/LAVORI/WDS/SELEN32/VSC/vsca_SANT.dat'
 ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ! # Parameters for imode=1 & imode=3 (declared but not used) 
 ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 integer, parameter ::      kv=0       ! DUMMY PARAMETER Type of viscosity profile
 integer, parameter ::      ioc=0      ! DUMMY PARAMETER PREM switch (1 oceanic, 0 continental)
 real, parameter::          lth=0      ! DUMMY PARAMETER Lithospheric thickness (m)
 character*30, parameter :: mod3='0'   ! DUMMY PARAMETER nla=2, user supplied model
 integer, parameter ::      ire=0      ! DUMMY PARAMETER Mantle rheology (Table 7)
 !~~~~~~~~~~~~~~~~~~~~~
 !# General parameters 
 !~~~~~~~~~~~~~~~~~~~~~
 integer, parameter :: nla=     3      ! Number of mantle layers
 character*1, parameter :: th= 'h'     ! Stands for 'Heaviside' Love numbers 
 integer, parameter :: iload=  1     ! Loading (1) or tidal (0) Love numbers 
 integer, parameter ::    L1=  1     ! Minimum degree
 integer, parameter ::    L2=  94      ! Maximum degree
 integer, parameter ::    LS=  1     ! Step degree
 integer, parameter ::    IH=  1     ! h Love number (0/1)
 integer, parameter ::    IL=  1     ! l Love number (0/1)
 integer, parameter ::    IK=  1     ! k Love number (0/1)
 integer, parameter ::   NSD=  64     ! Significant digits
 integer, parameter ::    NG=  8      ! Order of the Gaver method
 integer, parameter :: ISALZ=  1      !Salzer acceleration
 character*3, parameter :: tscale= 'lin'   ! toggle between lin/log time scale
 integer, parameter :: p=  21  ! Number of time points minus one
 real, parameter :: m1=0, m2=21  ! tmin and tmax
 ! 
 ! ********* File closed by config.f90 *********
 ! 
