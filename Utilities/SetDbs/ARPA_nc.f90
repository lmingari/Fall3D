!***************************************************************
!*
!*    Module for ARPA netcdf operations
!*
!***************************************************************
MODULE ARPA_nc
  use KindType, ONLY :  ip,rp
  implicit none
  save
  !
  !*** Database
  !
  integer(ip) :: ncID ,nDim, nVar, nAttr, nDimMax
  integer(ip) ::       iDim, iVar, iAttr
  !
  !*** Database arrays
  !
  character(len=50), allocatable :: dimName(:),varName(:)
  integer(ip)      , allocatable :: dimlen (:),varDims(:),varType(:),varDimID(:,:)
  integer(ip)      , allocatable :: ivoid  (:)
  !
  !*** ARPA dimension names
  !
  character(len=25) :: nt_ARPA_name = 'time'
  character(len=25) :: nx_ARPA_name = 'lon'
  character(len=25) :: ny_ARPA_name = 'lat'
  character(len=25) :: nz_ARPA_name = 'pres'
  !
  integer :: nx_ARPA_dimid   ! Dimension IDs
  integer :: ny_ARPA_dimid
  integer :: nz_ARPA_dimid
  integer :: nt_ARPA_dimid
  !
  integer(ip) :: nt_ARPA,it_ARPA
  integer(ip) :: nx_ARPA
  integer(ip) :: ny_ARPA
  integer(ip) :: nz_ARPA
  !
  !*** ARPA variable names
  !
  character(len=25) :: time_ARPA_name  = 'time'
  character(len=25) :: p_ARPA_name     = 'pres'
  character(len=25) :: lat_ARPA_name   = 'RLAT'
  character(len=25) :: lon_ARPA_name   = 'RLON'
  character(len=25) :: fi_ARPA_name    = 'GP'
  character(len=25) :: fis_ARPA_name   = 'GPsfc'
  character(len=25) :: T_ARPA_name     = 'T'
  character(len=25) :: Q_ARPA_name     = 'QV'
  character(len=25) :: u_ARPA_name     = 'U'
  character(len=25) :: v_ARPA_name     = 'V'
  character(len=25) :: w_ARPA_name     = 'OMEGA'
  !
  !*** ARPA attribute names
  !
  character(len=25) :: ibyr_ARPA_name = 'YEAR'
  character(len=25) :: ibmo_ARPA_name = 'MONTH'
  character(len=25) :: ibdy_ARPA_name = 'DAY'
  character(len=25) :: ibhr_ARPA_name = 'HOUR'
  character(len=25) :: dt_ARPA_name   = 'TIME_INCR'
  character(len=25) :: ARPA_cen_lon_name  = 'CEN_LON'
  character(len=25) :: ARPA_cen_lat_name  = 'CEN_LAT'
  !
  !*** ARPA variables ID
  !
  integer(ip) :: fiID,fisID,tID,qID,uID,vID,wID
  !
  !*** ARPA attribute values
  !
  integer(ip) :: ibyr_ARPA,ibmo_ARPA,ibdy_ARPA,ibhr_ARPA,ibmi_ARPA,ibse_ARPA
  integer(ip) :: dt_ARPA
  real   (rp) :: ARPA_cen_lon,ARPA_cen_lat
  !
  !*** ARPA variable arrays
  !
  real(rp), allocatable :: time_ARPA(:)
  real(rp), allocatable :: timesec_ARPA(:)
  !
  real(rp), allocatable :: p_ARPA(:)
  real(rp), allocatable :: lat_ARPA(:,:)
  real(rp), allocatable :: lon_ARPA(:,:)
  real(rp), allocatable :: topg_ARPA(:,:)
  real(rp), allocatable :: z_ARPA    (:,:,:)
  !
  !**  Working arrays
  !
  logical :: do_transpose  ! Flag: transpose matrices
  real(rp), allocatable :: work_ARPA1(:,:,:),work_ARPA2(:,:,:)
  !
  !*** Interpolation stuff
  !
  integer(ip) :: nelem_ARPA
  integer(ip) :: npoin_DAT
  !
  integer(ip), allocatable :: lnods_ARPA(:,:)
  integer(ip), allocatable :: lelpo_ARPA(:)
  !
  real   (rp), allocatable :: s_po_ARPA(:)
  real   (rp), allocatable :: t_po_ARPA(:)
  !
  !
  !
CONTAINS
  !
  subroutine ARPA_getmem01
    !***************************************************************
    !*
    !*   Gets memory for time-independent variables
    !*
    !***************************************************************
    implicit none
    !
    allocate( dimName (nDim) )
    allocate( varName (nVar) )
    !
    allocate( dimLen  (nDim) )
    allocate( varDims (nVar) )
    allocate( varType (nVar) )
    allocate( varDimID(nVar,nDimMax) )
    allocate( ivoid   (nDimMax) )
    !
  end subroutine ARPA_getmem01
  !
  !
  !
  subroutine ARPA_getmem02
    !***************************************************************
    !*
    !*   Gets memory for time-independent variables
    !*
    !***************************************************************
    implicit none
    !
    allocate(time_ARPA(nt_ARPA))
    allocate(timesec_ARPA(nt_ARPA))
    !
    !*** Variables
    !
    allocate(p_ARPA   (nz_ARPA))
    allocate(lat_ARPA (nx_ARPA,ny_ARPA))
    allocate(lon_ARPA (nx_ARPA,ny_ARPA))
    allocate(topg_ARPA(nx_ARPA,ny_ARPA))
    allocate(z_ARPA   (nx_ARPA,ny_ARPA,nz_ARPA))
    !
    !*** Interpolation
    !
    allocate(lnods_ARPA(4,nelem_ARPA))
    allocate(lelpo_ARPA(npoin_DAT))
    allocate(s_po_ARPA(npoin_DAT))
    allocate(t_po_ARPA(npoin_DAT))
    !
  end subroutine ARPA_getmem02
  !
END MODULE ARPA_nc
