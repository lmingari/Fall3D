!***************************************************************
!*
!*		Module for master operations
!*
!***************************************************************
MODULE Master
  use KindType
  IMPLICIT NONE
  SAVE
  !
  !***  Version of input file
  !
  real(rp),parameter :: version = 7.0
  !
  !***  Constants
  !
  real(rp), parameter :: pi     = 3.14159265358979323846_rp
  real(rp), parameter :: g      = 9.81_rp
  real(rp), parameter :: Rearth = 6356d3                  ! Earth's radius
  !
  !***  Master control
  !
  logical :: terrain_following = .true.        ! Source written in terrain following
  logical :: correct_vent      = .true.        ! Correct the dbs vent height
  logical :: standard_atmosphere = .false.     ! Use air properties form dbs or standard
  !
  character(len=s_work) :: use_mesh,model_name,type_source
  !
  !***  Time related variables
  !
  integer(ip) :: iebeg,ieend,irend,iebeg1,ieend1
  integer(ip) :: ndt,idt
  integer(ip) :: idbsbeg,idbsend,time,time2,idbsincr,idbsyr,idbsmo,idbsday
  !
  !***  Granulometry
  !
  integer(ip) :: npart,ngas,nc
  real   (rp) :: rhomean
  !
  !***  Aggregation (one class only)
  !
  character(len=s_work) :: type_aggr
  logical               :: aggregation = .false.
  real(rp)              :: phi_aggr,diam_aggr,rho_aggr,frac_aggr
  real(rp)              :: Df
  !
  !***  Volatiles
  !
  logical  :: SO2 = .false.
  logical  :: H2O = .false.
  real(rp) :: SO2_percentage
  real(rp) :: H2O_percentage
  !
  !***  Plume variables
  !
  logical               :: MER_wind_coupling
  character(len=s_work) :: MER_vs_h
  !
  character(len=10) :: solve_plume_for
  integer(ip)       :: np,ns,modv,MFR_iter
  real   (rp)       :: Hplume,Asuzu,Lsuzu,n_MFR(2)
  real   (rp)       :: erumass = 0.0
  !
  !***  Mesh related variables
  !
  character(len=10) :: coord_sys
  integer(ip) :: nx,ny,nz
  real   (rp) :: lonmin,latmin,lonmax,latmax
  real   (rp) :: xmin  ,ymin  ,xmax  ,ymax
  real   (rp) :: xorigr,yorigr,dx,dy,dX1,dX2
  !
  !***  Data at vent (BC)
  !
  character(len=3) :: UTM_zone
  real(rp) :: M0,u0,T0,W0,X0,Y0,Z0,PSIa,Tair0,Rair0,dZ0,R0,rho0,rhomax
  real(rp) :: X0_UTM,Y0_UTM,X0_LL,Y0_LL,LonPlume,LatPlume
  !
  !***  Resuspension related variables
  !
  logical           :: moisture_correction
  character(len=15) :: emission_scheme
  integer(ip)       :: nz_emission
  real   (rp)       :: diam_max,tust_cte,w_threshold,emission_factor,z_emission
  !
  !***  Arrays
  !
  integer(ip),allocatable :: idt_src(:)     ! (ndt)
  !
  real(rp), allocatable :: M0_dt    (:)
  real(rp), allocatable :: HPlume_dt(:)
  real(rp), allocatable :: Asuzu_dt (:)
  real(rp), allocatable :: Lsuzu_dt (:)
  real(rp), allocatable :: u0_dt    (:)
  real(rp), allocatable :: T0_dt    (:)
  real(rp), allocatable :: w0_dt    (:)
  !
  real(rp), allocatable :: fc  (:)         ! (nc)
  real(rp), allocatable :: rhop(:)
  real(rp), allocatable :: diam(:)
  real(rp), allocatable :: sphe(:)
  real(rp), allocatable :: psi (:)
  character(len=25), allocatable :: classname(:)
  !
  real(rp), allocatable :: Xplum(:)       ! (ns)
  real(rp), allocatable :: Yplum(:)
  real(rp), allocatable :: Zplum(:)
  real(rp), allocatable :: Lplum(:)
  real(rp), allocatable :: Hplum(:)
  real(rp), allocatable :: Uplum(:)
  real(rp), allocatable :: Tplum(:)
  real(rp), allocatable :: Dplum(:)
  real(rp), allocatable :: Rplum(:)
  real(rp), allocatable :: Mair (:)
  real(rp), allocatable :: Mplum(:,:)
  !
  real(rp), allocatable :: topg(:,:)
  real(rp), allocatable :: zmodel(:)
  real(rp), allocatable :: src(:,:,:,:)
  real(rp), allocatable :: Vair(:)        ! Velocity
  real(rp), allocatable :: Tair(:)        ! Temperature
  real(rp), allocatable :: Aair(:)        ! Azimuth
  real(rp), allocatable :: Rair(:)        ! Density
  real(rp), allocatable :: Nair(:)        ! N2 buoyancy frequency
  !
  real(rp), allocatable :: lon(:)           ! nx
  real(rp), allocatable :: lat(:)           ! ny
  real(rp), allocatable :: load(:,:)        ! (nx,ny)
  real(rp), allocatable :: ust (:,:)        ! (nx,ny)
  real(rp), allocatable :: smoi(:,:)        ! (nx,ny)
  real(rp), allocatable :: densi(:,:)       ! (nx,ny)         Air density
  real(rp), allocatable :: Hm1 (:,:)        ! (nx,ny)
  real(rp), allocatable :: emis(:,:,:)      ! (nx,ny,nc    )  Emission (source term)
  real(rp), allocatable :: emis_mas(:,:)    ! (nx,ny       )  Total Emited mass (source term)
  real(rp), allocatable :: clas(:,:,:)      ! (nx,ny,nc_dep)  Class load
  real(rp), allocatable :: tust(:,:,:)      ! (nx,ny,nc_dep)  Threshold friction velocity

  !
  !
CONTAINS
  !
  !
  !
  subroutine upcase(word)
    !***********************************************************************
    !*
    !*    This routine converts word to upper case
    !*
    !***********************************************************************
    implicit none
    character(len=*) :: word
    integer(ip)      :: iposi,ioctv,item1,item2,item3
    integer(ip)      :: ilen
    !
    item1 = int(o'141')
    item2 = int(o'172')
    item3 = int(o'40')
    !
    ilen=LEN_TRIM(word)
    !
    do iposi=1,ilen                                      ! process all positions
       ioctv=ichar(word(iposi:iposi))                    ! octal value
       if(item1.le.ioctv.and.item2.ge.ioctv) then        ! it is a lower case
          ioctv=ioctv-item3                              ! equivalent upper case
          word(iposi:iposi)=char(ioctv)                  ! convert it to upcase
       end if
    end do ! iposi=1,ilen
    !
  end subroutine upcase
  !
  !
  !
  real(rp)  function erff(x)
    !************************************************************************
    !*
    !*    Computes the error function (area below the typifiyed normal
    !*    distribution)
    !*
    !************************************************************************
    use KindType
    implicit none
    real(rp) ::  x
    !
    logical  :: go_on
    real(rp) :: t1,t2,dt
    real(rp), parameter ::  twopi = 6.283185307179586
    !
    !***  Computes integral between t=0 and t=abs(x)
    !
    if(x.eq.0) then
       erff=0.5
       return
    endif
    erff = 0.0
    dt  = abs(x)*1d-5
    !
    t1 = 0.0
    t2 = dt
    do while(go_on)
       erff = erff + 0.5*(exp(-t1*t1/2.0)+exp(-t2*t2/2.0))*dt
       t1 = t2
       t2 = t2 + dt
       if(t2.ge.abs(x)) go_on = .false.
    end do
    erff = erff/sqrt(twopi)

    if(x.ge.0.0) then
       erff = 0.5 + erff
    else
       erff = 0.5 - erff
    end if

    return
  end function erff

  !
END MODULE Master
