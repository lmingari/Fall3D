!***************************************************************
!*
!*    Module for GFS netcdf operations
!*
!***************************************************************
MODULE GFS_nc
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
  !*** GFS dimension names
  !
  character(len=25) :: nt_GFS_name = 'time'
  character(len=25) :: nx_GFS_name = 'lon'
  character(len=25) :: ny_GFS_name = 'lat'
  character(len=25) :: nz_GFS_name = 'pres'
  !
  integer(ip) :: nt_GFS,it_GFS
  integer(ip) :: nx_GFS
  integer(ip) :: ny_GFS
  integer(ip) :: nz_GFS
  !
  !*** GFS variable names
  !
  character(len=25) :: time_GFS_name  = 'time'
  character(len=25) :: p_GFS_name     = 'pres'
  character(len=25) :: fis_GFS_name   = 'HGT:surface'
  character(len=25) :: lmask_GFS_name = 'LAND'
  character(len=25) :: pblh_GFS_name  = 'HPBL'
  character(len=25) :: fi_GFS_name    = 'HGT'
  character(len=25) :: T_GFS_name     = 'TMP'
  character(len=25) :: Q_GFS_name     = 'RH'
  character(len=25) :: u_GFS_name     = 'UGRD'
  character(len=25) :: v_GFS_name     = 'VGRD'
  character(len=25) :: w_GFS_name     = 'VVEL'
  !
  !*** GFS attribute names
  !
  character(len=25) :: ibyr_GFS_name = 'YEAR'
  character(len=25) :: ibmo_GFS_name = 'MONTH'
  character(len=25) :: ibdy_GFS_name = 'DAY'
  character(len=25) :: ibhr_GFS_name = 'HOUR'
  character(len=25) :: dt_GFS_name   = 'TIME_INCR'
  character(len=25) :: lonmin_GFS_name = 'LONMIN'
  character(len=25) :: lonmax_GFS_name = 'LONMAX'
  character(len=25) :: latmin_GFS_name = 'LATMIN'
  character(len=25) :: latmax_GFS_name = 'LATMAX'
  !
  !*** GFS variables ID
  !
  integer(ip) :: fiID,fisID,tID,qID,uID,vID,wID,lmaskID,pblhID
  !
  !*** GFS attribute values
  !
  integer(ip) :: ibyr_GFS,ibmo_GFS,ibdy_GFS,ibhr_GFS,ibmi_GFS,ibse_GFS
  integer(ip) :: dt_GFS
  real   (rp) :: lonmin_GFS,lonmax_GFS,latmin_GFS,latmax_GFS
  !
  !*** GFS variable arrays
  !
  real(rp), allocatable :: time_GFS(:)
  real(rp), allocatable :: timesec_GFS(:)
  !
  real(rp), allocatable :: p_GFS(:)
  real(rp), allocatable :: lat_GFS (:,:)
  real(rp), allocatable :: lon_GFS (:,:)
  real(rp), allocatable :: topg_GFS(:,:)
  real(rp), allocatable :: z_GFS   (:,:,:)
  !
  !**  Working arrays
  !
  real(rp), allocatable :: work_GFS1(:,:,:),work_GFS2(:,:,:)
  !
  !*** Interpolation stuff
  !
  integer(ip) :: nelem_GFS
  integer(ip) :: npoin_DAT
  !
  integer(ip), allocatable :: lnods_GFS(:,:)
  integer(ip), allocatable :: lelpo_GFS(:)
  !
  real   (rp), allocatable :: s_po_GFS(:)
  real   (rp), allocatable :: t_po_GFS(:)
  !
  !
  !
CONTAINS
  !
  subroutine GFS_getmem01
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
  end subroutine GFS_getmem01
  !
  !
  !
  subroutine GFS_getmem02
    !***************************************************************
    !*
    !*   Gets memory for time-independent variables
    !*
    !***************************************************************
    implicit none
    !
    allocate(time_GFS(nt_GFS))
    allocate(timesec_GFS(nt_GFS))
    !
    !*** Variables
    !
    allocate(p_GFS   (nz_GFS))
    allocate(lat_GFS (nx_GFS,ny_GFS))
    allocate(lon_GFS (nx_GFS,ny_GFS))
    allocate(topg_GFS(nx_GFS,ny_GFS))
    allocate(z_GFS   (nx_GFS,ny_GFS,nz_GFS))
    !
    !*** Interpolation
    !
    allocate(lnods_GFS(4,nelem_GFS))
    allocate(lelpo_GFS(npoin_DAT))
    allocate(s_po_GFS(npoin_DAT))
    allocate(t_po_GFS(npoin_DAT))
    !
  end subroutine GFS_getmem02
  !
END MODULE GFS_nc
