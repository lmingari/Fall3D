!***************************************************************
!*
!*		Module for master operations
!*
!***************************************************************
MODULE Master
  use KindType
  implicit none
  save
  !
  !***  Code version
  !
  real(rp), parameter :: version = 7.0
  !
  !***  Number of CPU groups
  !
  character(len=2) :: c_cpu_groups              ! A max of 99 groups allowed
  integer  (ip)    :: number_of_cpu_groups
  !
  !***  Time related variables
  !
  logical :: v_drydep   = .true.
  logical :: meteotime  = .true.
  logical :: sourcetime = .true.
  logical :: restart    = .false.
  !
  character(len=3), dimension(12) :: month= (/'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC' /)
  integer(ip) :: ibyr,ibmo,ibdy
  integer(ip) :: iiter
  integer(ip) :: ibdbs,iedbs,idbsyr,idbsmo,idbsdy
  !
  real   (rp) :: time,beg_time,end_time,print_time,meteo_time,source_time,time_lag
  !
  real   (rp) :: safety_factor,dt,dtc
  real   (rp) :: erumass,flowrate,flowrate_part,flowrate_gas,outmass
  real   (rp) :: massvol,massdep,massdry,masswet
  real   (rp) :: cpu_time_begin,cpu_time_end
  !
  !***  Particle properties
  !
  integer(ip), parameter              :: ncmax = 50
  integer(ip)                         :: nc,np,ng,nsou
  real   (rp)      , dimension(ncmax) :: rhop,diam,sphe,psi,fc
  character(len=25), dimension(ncmax) :: classname
  !
  !*** Wet deposition
  !
  logical :: wetdeposition = .false.
  !
  !***  Aggregation
  !
  logical     :: aggregation     = .false.  ! Switch-on for Costa et al. (2010) aggregation model
  logical     :: aggregationtime = .false.
  integer     :: ndtagr
  integer(ip) :: iaggr
  real   (rp) :: Df,vset_factor
  real   (rp) :: Alagrmax    ,Alagrmin    ,Alagrave, &
       Alagrmax_ice,Alagrmin_ice,          &
       Alagrmax_wat,Alagrmin_wat,          &
       Abagrmax    ,Abagrmin    ,Abagrave, &
       Asagrmax    ,Asagrmin    ,Asagrave, &
       Adagrmax    ,Adagrmin    ,Adagrave
  real   (rp) :: Kbagrmax,Kbagrmin, &
       Ksagrmax,Ksagrmin, &
       Kdagrmax,Kdagrmin
  !
  !*** Gravity current model
  !
  logical          :: gravity_current  = .false.
  character(len=3) :: UTM_zone
  real(rp)         :: xvent_ll,yvent_ll,xvent_utm,yvent_utm
  !
  real(rp) :: c_flow_rate                       ! Constant                (Read from file.inp)
  real(rp) :: lambda_grav                       ! Constant                (Read from file.inp)
  real(rp) :: k_entrain                         ! Entrainment coefficient (Read from file.inp)
  real(rp) :: brunt_vaisala                     ! Brunt-Vaisala frequency (Read from file.inp)
  real(rp) :: radial_wind_lon,radial_wind_lat   ! (Same as vent coordinates)
  real(rp) :: radial_wind_x  ,radial_wind_y     ! (Same as vent coordinates)
  !
  !***  Memory storage
  !
  integer(ip) :: memo
  !
  !***  Integer arrays
  !
  integer(ip), allocatable :: isrc(:)              ! isrc(nsou)
  integer(ip), allocatable :: jsrc(:)              ! jsrc(nsou)
  integer(ip), allocatable :: ksrc(:)              ! ksrc(nsou)
  !
  !***  Real arrays
  !
  real(rp), allocatable :: Hm1   (:,:)     ! Hm1   (nx,ny)  Scaling factor
  real(rp), allocatable :: zlayer(:)       ! zlayer(nz)
  real(rp), allocatable :: dz    (:)       ! dz    (0:nz)
  real(rp), allocatable :: topg  (:,:)     ! topg  (nx,ny)
  real(rp), allocatable :: LDU   (:,:)     ! ldu(nx,ny) land use
  real(rp), allocatable :: vx    (:,:,:)   ! vx    (nx,ny,nz)
  real(rp), allocatable :: vy    (:,:,:)   ! vy    (nx,ny,nz)
  real(rp), allocatable :: vz    (:,:,:)   ! vz    (nx,ny,nz)
  real(rp), allocatable :: divu  (:,:,:)   ! divu  (nx,ny,nz)
  real(rp), allocatable :: tempe (:,:,:)   ! tempe (nx,ny,nz)
  real(rp), allocatable :: vptemp(:,:,:)   ! vptemp(nx,ny,nz)
  real(rp), allocatable :: p     (:,:,:)   ! p     (nx,ny,nz)
  real(rp), allocatable :: qv    (:,:,:)   ! qv    (nx,ny,nz)
  real(rp), allocatable :: ustar (:,:)     ! ustar (nx,ny)
  real(rp), allocatable :: zi    (:,:)     ! zi    (nx,ny)
  real(rp), allocatable :: rl    (:,:)     ! rl    (nx,ny)
  real(rp), allocatable :: znt   (:,:)     ! znt   (nx,ny)
  real(rp), allocatable :: prate (:,:)     ! prate (nx,ny)
  !
  real(rp), allocatable :: rkhor (:,:,:)   ! rkhor (nx,ny,nz)
  real(rp), allocatable :: rkh1  (:,:,:)   ! rkh1  (nx,ny,nz) (scaled)
  real(rp), allocatable :: rkver (:,:,:)   ! rkver (nx,ny,nz)
  real(rp), allocatable :: rho   (:,:,:)   ! rho   (nx,ny,nz)
  real(rp), allocatable :: vset  (:,:,:,:) ! vset  (nx,ny,nz,nc)
  real(rp), allocatable :: vdep  (:,:,:)   ! vdep  (nx,ny,nc)
  real(rp), allocatable :: c     (:,:,:,:) ! c     (0:nx+1,0:ny+1,0:nz+1,nc)
  real(rp), allocatable :: cwet  (:,:,:)   ! cwet  (nx,ny,0:nc)
  real(rp), allocatable :: cdry  (:,:,:)   ! cdry  (nx,ny,0:nc)
  real(rp), allocatable :: cload (:,:,:)   ! cload (nx,ny,0:nc)
  real(rp), allocatable :: work  (:,:)     ! work  (max(nx+1,ny+1),2)
  !
  real(rp), allocatable :: tarat (:,:,:,:) ! tarat (nx,ny,nz,nc)   Aggregation rate
  real(rp), allocatable :: Alagr (:,:,:)   ! Alagr (nx,ny,nz)      alfa mean  aggregation
  real(rp), allocatable :: Abagr (:,:,:)   ! Abagr (nx,ny,nz)      A Brownian aggregation
  real(rp), allocatable :: Asagr (:,:,:)   ! Asagr (nx,ny,nz)      A Shear    aggregation
  real(rp), allocatable :: Adagr (:,:,:)   ! Asagr (nx,ny,nz)      A Differential settling aggregation
  real(rp), allocatable :: Maggr (:,:,:)   ! Maggr (nx,ny,nz)      Cumulative mass of aggregates
  real(rp), allocatable :: tplume(:,:,:)   ! tplume(nx,ny,nz)      Plume temperature (for aggregation)
  real(rp), allocatable :: wplume(:,:,:)   ! wplume(nx,ny,nz)      Plume water rate  (for aggregation)
  real(rp), allocatable :: mclass(:,:)     ! mclass(nc,2)
  !
  real(rp), allocatable :: xs    (:)       ! xs    (nsou)
  real(rp), allocatable :: ys    (:)       ! ys    (nsou)
  real(rp), allocatable :: zs    (:)       ! zs    (nsou)
  real(rp), allocatable :: tmrat (:,:)     ! tmrat (nsou,nc)
  !
CONTAINS
  !
  subroutine getmem1
    !**********************************************************************
    !*
    !*  Allocates memory for arrays related to nx,ny,nz,nc
    !*  except for zlayer and dz which are allocated during the lecture
    !*  of the input file
    !*
    !**********************************************************************
    use KindType
    use Numeric
    implicit none
    !
    integer :: rsize
    if(rp == 4) then
       rsize = 4      ! Assume that sizeof(real(4)) = 4 bytes
    else
       rsize = 8      ! Assume that sizeof(real(8)) = 8 bytes
    end if
    !
    allocate(Hm1(nx,ny))
    Hm1 = 0.0_rp
    memo = memo + rsize*nx*ny

    allocate(zlayer(nz))
    zlayer = 0.0_rp
    memo = memo + rsize*nz

    allocate(dz(0:nz))
    dz = 0.0_rp
    memo = memo + rsize*(nz+1)

    allocate(topg(nx,ny))
    topg = 0.0_rp
    memo = memo + rsize*nx*ny

    allocate(LDU(nx,ny))
    LDU = -99.
    memo = memo + rsize*nx*ny

    allocate(vx(nx,ny,nz))
    vx = 0.0_rp
    memo = memo + rsize*nx*ny*nz

    allocate(vy(nx,ny,nz))
    vy = 0.0_rp
    memo = memo + rsize*nx*ny*nz

    allocate(vz(nx,ny,nz))
    vz = 0.0_rp
    memo = memo + rsize*nx*ny*nz

    allocate(divu(nx,ny,nz))
    divu = 0.0_rp
    memo = memo + rsize*nx*ny*nz

    allocate(tempe(nx,ny,nz))
    tempe = 0.0_rp
    memo = memo + rsize*nx*ny*nz

    allocate(vptemp(nx,ny,nz))
    vptemp = 0.0_rp
    memo = memo + rsize*nx*ny*nz

    allocate(p(nx,ny,nz))
    p = 0.0_rp
    memo = memo + rsize*nx*ny*nz

    allocate(qv(nx,ny,nz))
    qv = 0.0_rp
    memo = memo + rsize*nx*ny*nz

    allocate(ustar(nx,ny))
    ustar = 0.0_rp
    memo = memo + rsize*nx*ny

    allocate(zi(nx,ny))
    zi = 0.0_rp
    memo = memo + rsize*nx*ny

    allocate(rl(nx,ny))
    rl = 0.0_rp
    memo = memo + rsize*nx*ny

    allocate(znt(nx,ny))
    znt = 0.0_rp
    memo = memo + rsize*nx*ny

    allocate(rkhor(nx,ny,nz))
    rkhor = 0.0_rp
    memo = memo + rsize*nx*ny*nz

    allocate(rkh1(nx,ny,nz))
    rkh1 = 0.0_rp
    memo = memo + rsize*nx*ny*nz

    allocate(rkver(nx,ny,nz))
    rkver = 0.0_rp
    memo = memo + rsize*nx*ny*nz

    allocate(rho(nx,ny,nz))
    rho = 0.0_rp
    memo = memo + rsize*nx*ny*nz

    allocate(vset(nx,ny,nz,nc))
    vset = 0.0_rp
    memo = memo + rsize*nx*ny*nz*nc

    allocate(vdep(nx,ny,nc))
    vdep = 0.0_rp
    memo = memo + rsize*nx*ny*nc

    allocate(c(0:nx+1,0:ny+1,0:nz+1,nc))
    c = 0.0_rp
    memo = memo + rsize*(nx+2)*(ny+2)*(nz+2)*nc

    allocate(cwet(nx,ny,0:nc))
    cwet = 0.0_rp
    memo = memo + rsize*nx*ny*(nc+1)

    allocate(cdry(nx,ny,0:nc))
    cdry = 0.0_rp
    memo = memo + rsize*nx*ny*(nc+1)

    allocate(cload(nx,ny,0:nc))
    cload = 0.0_rp
    memo = memo + rsize*nx*ny*(nc+1)

    allocate(work(max(nx+1,ny+1),2))
    work = 0.0_rp
    memo = memo + rsize*max(nx+1,ny+1)*2

    allocate(prate(nx,ny))
    prate = 0.0_rp
    memo  = memo + rsize*nx*ny
    !
    !** Aggregation (np>1)
    !
    if(aggregation) then
       allocate(tarat(nx,ny,nz,nc))
       tarat = 0.0_rp
       memo = memo + rsize*nx*ny*nz*np

       allocate(Alagr(nx,ny,nz))
       Alagr = 0.0_rp
       memo = memo + rsize*nx*ny*nz

       allocate(Abagr(nx,ny,nz))
       Abagr = 0.0_rp
       memo = memo + rsize*nx*ny*nz

       allocate(Asagr(nx,ny,nz))
       Asagr = 0.0_rp
       memo = memo + rsize*nx*ny*nz

       allocate(Adagr(nx,ny,nz))
       Adagr = 0.0_rp
       memo = memo + rsize*nx*ny*nz

       allocate(Maggr(nx,ny,nz))
       Maggr = 0.0_rp
       memo = memo + rsize*nx*ny*nz

       allocate(tplume(nx,ny,nz))
       tplume = 0.0_rp
       memo = memo + 8*nx*ny*nz   

       allocate(wplume(nx,ny,nz))
       wplume = 0.0_rp
       memo = memo + 8*nx*ny*nz 

       allocate(mclass(nc,2))
       mclass = 0.0_rp
       memo = memo + rsize*2*np

    end if
    !
  end subroutine getmem1
  !
  subroutine getmem2
    !**********************************************************************
    !*
    !*  Allocates memory for arrays related to nsou
    !*
    !**********************************************************************
    implicit none
    !
    allocate(isrc(nsou))
    isrc = 0

    allocate(jsrc(nsou))
    jsrc = 0

    allocate(ksrc(nsou))
    ksrc = 0

    allocate(xs(nsou))
    xs = 0.0_rp

    allocate(ys(nsou))
    ys = 0.0_rp

    allocate(zs(nsou))
    zs = 0.0_rp

    allocate(tmrat(nsou,nc))
    tmrat = 0.0_rp
    !
  end subroutine getmem2
  !
  subroutine outmem2
    !**********************************************************************
    !*
    !*  Deallocates memory for arrays related to nsou
    !*
    !**********************************************************************
    implicit none
    integer(ip) :: ierror
    !
    deallocate(isrc,STAT=ierror)
    if(ierror /= 0) call runend('Error deallocating memory for isrc')
    !
    deallocate(jsrc,STAT=ierror)
    if(ierror /= 0) call runend('Error deallocating memory for jsrc')
    !
    deallocate(ksrc,STAT=ierror)
    if(ierror /= 0) call runend('Error deallocating memory for ksrc')
    !
    deallocate(xs,STAT=ierror)
    if(ierror /= 0) call runend('Error deallocating memory for xs')
    !
    deallocate(ys,STAT=ierror)
    if(ierror /= 0) call runend('Error deallocating memory for ys')
    !
    deallocate(zs,STAT=ierror)
    if(ierror /= 0) call runend('Error deallocating memory for zs')
    !
    deallocate(tmrat,STAT=ierror)
    if(ierror /= 0) call runend('Error deallocating memory for tmrat')
    !
  end subroutine outmem2
  !
END MODULE Master
