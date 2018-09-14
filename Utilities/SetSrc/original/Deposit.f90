!***************************************************************
!*
!*		Module for deposit
!*
!***************************************************************
MODULE Deposit
  use KindType
  IMPLICIT NONE
  SAVE
  !
  !***  Despoit mesh related variables
  !
  character(len=10) :: coord_sys_dep
  integer(ip)       :: nx_dep,ny_dep,nt_dep
  real(rp)          :: xorigr_dep,yorigr_dep,dx_dep,dy_dep
  !
  !***  Granulometry
  !
  integer(ip), parameter :: ncmax = 50
  integer(ip) :: nc_dep
  !
  !***  netCDF
  !
  integer(ip)  :: ncID,wordID
  !
  character(len=25 ) :: attr_units
  character(len=100) :: attr_title,attr_desc
  !
  character(len=25) :: RES_name
  character(len=25) :: lon_RES_name    = 'lon'
  character(len=25) :: lat_RES_name    = 'lat'
  character(len=25) :: nc_dep_RES_name = 'nc_dep'
  character(len=25) :: nc_RES_name     = 'nc'
  character(len=25) :: time_RES_name   = 'time'
  !
  character(len=25)  :: load_RES_name        = 'LOAD'
  character(len=25)  :: clas_RES_name(ncmax) = 'LOAD_CLASS'   ! LOAD_CLASSxx
  character(len=25)  :: ust_RES_name         = 'UST'
  character(len=25)  :: smoi_RES_name        = 'SMOI'
  character(len=25)  :: densi_RES_name       = 'RHO'
  character(len=25)  :: tust_RES_name        = 'TUST'
  character(len=25)  :: emis_RES_name        = 'FLUX'
  !
  integer(ip)  :: nx_RES_ID,ny_RES_ID,nt_RES_ID,nc_RES_ID,nc_dep_RES_ID
  integer(ip)  :: lon_RES_ID
  integer(ip)  :: lat_RES_ID
  integer(ip)  :: time_RES_ID
  integer(ip)  :: load_RES_ID
  integer(ip)  :: clas_RES_ID(ncmax)
  integer(ip)  :: ust_RES_ID
  integer(ip)  :: smoi_RES_ID
  integer(ip)  :: densi_RES_ID
  integer(ip)  :: tust_RES_ID
  integer(ip)  :: emis_RES_ID
  !
  !***  Arrays
  !
  real(rp), allocatable :: rhop_dep(:)       ! (nc_dep)
  real(rp), allocatable :: diam_dep(:)
  real(rp), allocatable :: sphe_dep(:)
  !
  real(rp), allocatable :: lon_dep(:)        ! (nx_dep)
  real(rp), allocatable :: lat_dep(:)        ! (ny_dep)
  !
  real(rp), allocatable :: work_dep(:,:)     ! (nx_dep,ny_dep)
  real(rp), allocatable :: load_dep(:,:)     ! (nx_dep,ny_dep)
  real(rp), allocatable :: clas_dep(:,:,:)   ! (nx_dep,ny_dep,nc_dep)
  !
END MODULE Deposit
