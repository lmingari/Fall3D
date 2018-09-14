!***************************************************************
!*
!*    Module for ECMWF netcdf operations
!*
!***************************************************************
MODULE ECMWF_nc
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
  !*** ECMWF dimension names
  !
  character(len=25) :: nt_ECMWF_name = 'time'
  character(len=25) :: nx_ECMWF_name = 'lon'
  character(len=25) :: ny_ECMWF_name = 'lat'
  character(len=25) :: nz_ECMWF_name = 'pres'
  !
  integer(ip) :: nt_ECMWF,it_ECMWF
  integer(ip) :: nx_ECMWF
  integer(ip) :: ny_ECMWF
  integer(ip) :: nz_ECMWF
  !
  !*** ECMWF variable names
  !
  character(len=25) :: time_ECMWF_name  = 'time'
  character(len=25) :: p_ECMWF_name     = 'pres'
  character(len=25) :: fis_ECMWF_name   = 'Z:sfc'
  character(len=25) :: us_ECMWF_name    = 'U10:sfc'
  character(len=25) :: vs_ECMWF_name    = 'V10:sfc'
  character(len=25) :: Ts_ECMWF_name    = 'T2:sfc'
  character(len=25) :: lmask_ECMWF_name = 'LSM:sfc'
  character(len=25) :: pblh_ECMWF_name  = 'BLH:sfc'
  !
  character(len=25) :: fi_ECMWF_name    = 'Z'
  character(len=25) :: T_ECMWF_name     = 'T'
  character(len=25) :: Q_ECMWF_name     = 'R'
  character(len=25) :: u_ECMWF_name     = 'U'
  character(len=25) :: v_ECMWF_name     = 'V'
  character(len=25) :: w_ECMWF_name     = 'W'
  !
  !*** ECMWF attribute names
  !
  character(len=25) :: ibyr_ECMWF_name = 'YEAR'
  character(len=25) :: ibmo_ECMWF_name = 'MONTH'
  character(len=25) :: ibdy_ECMWF_name = 'DAY'
  character(len=25) :: ibhr_ECMWF_name = 'HOUR'
  character(len=25) :: dt_ECMWF_name   = 'TIME_INCR'
  character(len=25) :: lonmin_ECMWF_name = 'LONMIN'
  character(len=25) :: lonmax_ECMWF_name = 'LONMAX'
  character(len=25) :: latmin_ECMWF_name = 'LATMIN'
  character(len=25) :: latmax_ECMWF_name = 'LATMAX'
  !
  !*** ECMWF variables ID
  !
  integer(ip) :: fiID,fisID,tID,qID,uID,vID,wID,lmaskID,pblhID,usID,vsID,TsID
  !
  !*** ECMWF attribute values
  !
  integer(ip) :: ibyr_ECMWF,ibmo_ECMWF,ibdy_ECMWF,ibhr_ECMWF,ibmi_ECMWF,ibse_ECMWF
  integer(ip) :: dt_ECMWF
  real   (rp) :: lonmin_ECMWF,lonmax_ECMWF,latmin_ECMWF,latmax_ECMWF
  !
  !*** ECMWF variable arrays
  !
  real(rp), allocatable :: time_ECMWF(:)
  real(rp), allocatable :: timesec_ECMWF(:)
  !
  real(rp), allocatable :: p_ECMWF(:)
  real(rp), allocatable :: lat_ECMWF (:,:)
  real(rp), allocatable :: lon_ECMWF (:,:)
  real(rp), allocatable :: topg_ECMWF(:,:)
  real(rp), allocatable :: z_ECMWF   (:,:,:)
  !
  !**  Working arrays
  !
  real(rp), allocatable :: work_ECMWF1(:,:,:),work_ECMWF2(:,:,:)
  !
  !*** Interpolation stuff
  !
  integer(ip) :: nelem_ECMWF
  integer(ip) :: npoin_DAT
  !
  integer(ip), allocatable :: lnods_ECMWF(:,:)
  integer(ip), allocatable :: lelpo_ECMWF(:)
  !
  real   (rp), allocatable :: s_po_ECMWF(:)
  real   (rp), allocatable :: t_po_ECMWF(:)
  !
  !
  !
CONTAINS
  !
  subroutine ECMWF_getmem01
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
  end subroutine ECMWF_getmem01
  !
  !
  !
  subroutine ECMWF_getmem02
    !***************************************************************
    !*
    !*   Gets memory for time-independent variables
    !*
    !***************************************************************
    implicit none
    !
    allocate(time_ECMWF(nt_ECMWF))
    allocate(timesec_ECMWF(nt_ECMWF))
    !
    !*** Variables
    !
    allocate(p_ECMWF   (nz_ECMWF))
    allocate(lat_ECMWF (nx_ECMWF,ny_ECMWF))
    allocate(lon_ECMWF (nx_ECMWF,ny_ECMWF))
    allocate(topg_ECMWF(nx_ECMWF,ny_ECMWF))
    allocate(z_ECMWF   (nx_ECMWF,ny_ECMWF,nz_ECMWF))
    !
    !*** Interpolation
    !
    allocate(lnods_ECMWF(4,nelem_ECMWF))
    allocate(lelpo_ECMWF(npoin_DAT))
    allocate(s_po_ECMWF(npoin_DAT))
    allocate(t_po_ECMWF(npoin_DAT))
    !
  end subroutine ECMWF_getmem02
  !
END MODULE ECMWF_nc
