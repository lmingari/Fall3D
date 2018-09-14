!***************************************************************
!*
!*    Module for ETA netcdf operations
!*
!***************************************************************
MODULE ETA_nc
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
  !*** ETA dimension names
  !
  character(len=25) :: nt_ETA_name = 'time'
  character(len=25) :: nx_ETA_name = 'lon'
  character(len=25) :: ny_ETA_name = 'lat'
  character(len=25) :: nz_ETA_name = 'pres'
  !
  integer(ip) :: nt_ETA,it_ETA
  integer(ip) :: nx_ETA
  integer(ip) :: ny_ETA
  integer(ip) :: nz_ETA
  !
  !*** ETA variable names
  !
  character(len=25) :: time_ETA_name  = 'time'
  character(len=25) :: p_ETA_name     = 'pres'
  character(len=25) :: fis_ETA_name   = 'HGT:sfc'
  character(len=25) :: lmask_ETA_name = 'LAND'
  character(len=25) :: fi_ETA_name    = 'HGT'
  character(len=25) :: T_ETA_name     = 'TMP'
  character(len=25) :: Q_ETA_name     = 'RH'
  character(len=25) :: u_ETA_name     = 'UGRD'
  character(len=25) :: v_ETA_name     = 'VGRD'
  character(len=25) :: w_ETA_name     = 'VVEL'
  !
  !*** ETA attribute names
  !
  character(len=25) :: ibyr_ETA_name = 'YEAR'
  character(len=25) :: ibmo_ETA_name = 'MONTH'
  character(len=25) :: ibdy_ETA_name = 'DAY'
  character(len=25) :: ibhr_ETA_name = 'HOUR'
  character(len=25) :: dt_ETA_name   = 'TIME_INCR'
  character(len=25) :: lonmin_ETA_name = 'LONMIN'
  character(len=25) :: lonmax_ETA_name = 'LONMAX'
  character(len=25) :: latmin_ETA_name = 'LATMIN'
  character(len=25) :: latmax_ETA_name = 'LATMAX'
  !
  !*** ETA variables ID
  !
  integer(ip) :: fiID,fisID,tID,qID,uID,vID,wID,lmaskID,pblhID
  !
  !*** ETA attribute values
  !
  integer(ip) :: ibyr_ETA,ibmo_ETA,ibdy_ETA,ibhr_ETA,ibmi_ETA,ibse_ETA
  integer(ip) :: dt_ETA
  real   (rp) :: lonmin_ETA,lonmax_ETA,latmin_ETA,latmax_ETA
  !
  !*** ETA variable arrays
  !
  real(rp), allocatable :: time_ETA(:)
  real(rp), allocatable :: timesec_ETA(:)
  !
  real(rp), allocatable :: p_ETA(:)
  real(rp), allocatable :: lat_ETA (:,:)
  real(rp), allocatable :: lon_ETA (:,:)
  real(rp), allocatable :: topg_ETA(:,:)
  real(rp), allocatable :: z_ETA   (:,:,:)
  !
  !**  Working arrays
  !
  real(rp), allocatable :: work_ETA1(:,:,:),work_ETA2(:,:,:)
  !
  !*** Interpolation stuff
  !
  integer(ip) :: nelem_ETA
  integer(ip) :: npoin_DAT
  !
  integer(ip), allocatable :: lnods_ETA(:,:)
  integer(ip), allocatable :: lelpo_ETA(:)
  !
  real   (rp), allocatable :: s_po_ETA(:)
  real   (rp), allocatable :: t_po_ETA(:)
  !
  !
  !
CONTAINS
  !
  subroutine ETA_getmem01
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
  end subroutine ETA_getmem01
  !
  !
  !
  subroutine ETA_getmem02
    !***************************************************************
    !*
    !*   Gets memory for time-independent variables
    !*
    !***************************************************************
    implicit none
    !
    allocate(time_ETA(nt_ETA))
    allocate(timesec_ETA(nt_ETA))
    !
    !*** Variables
    !
    allocate(p_ETA   (nz_ETA))
    allocate(lat_ETA (nx_ETA,ny_ETA))
    allocate(lon_ETA (nx_ETA,ny_ETA))
    allocate(topg_ETA(nx_ETA,ny_ETA))
    allocate(z_ETA   (nx_ETA,ny_ETA,nz_ETA))
    !
    !*** Interpolation
    !
    allocate(lnods_ETA(4,nelem_ETA))
    allocate(lelpo_ETA(npoin_DAT))
    allocate(s_po_ETA(npoin_DAT))
    allocate(t_po_ETA(npoin_DAT))
    !
  end subroutine ETA_getmem02
  !
END MODULE ETA_nc
