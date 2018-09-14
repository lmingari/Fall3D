!***************************************************************
!*
!*    Module for WRF netcdf operations
!*
!***************************************************************
MODULE WRF_nc
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
  !*** WRF dimension names
  !
  character(len=25) :: nt_WRF_name     = 'Time'
  character(len=25) :: nx_WRF_name     = 'west_east'
  character(len=25) :: nxstag_WRF_name = 'west_east_stag'
  character(len=25) :: ny_WRF_name     = 'south_north'
  character(len=25) :: nystag_WRF_name = 'south_north_stag'
  character(len=25) :: nz_WRF_name     = 'bottom_top'
  character(len=25) :: nzstag_WRF_name = 'bottom_top_stag'
  character(len=25) :: nsoil_WRF_name  = 'soil_layers_stag'
  !
  integer(ip) :: nt_WRF,it_WRF
  integer(ip) :: nx_WRF,nxstag_WRF
  integer(ip) :: ny_WRF,nystag_WRF
  integer(ip) :: nz_WRF,nzstag_WRF
  integer(ip) :: nsoil_WRF
  !
  !*** WRF variable names
  !
  character(len=25) :: time_WRF_name  = 'Times'
  character(len=25) :: lat_WRF_name   = 'XLAT'
  character(len=25) :: lon_WRF_name   = 'XLONG'
  character(len=25) :: hgt_WRF_name   = 'HGT'
  character(len=25) :: lmask_WRF_name = 'LANDMASK'
  character(len=25) :: phb_WRF_name   = 'PHB'
  character(len=25) :: ph_WRF_name    = 'PH'
  character(len=25) :: u_WRF_name     = 'U'
  character(len=25) :: v_WRF_name     = 'V'
  character(len=25) :: W_WRF_name     = 'W'
  character(len=25) :: p_WRF_name     = 'P'
  character(len=25) :: pb_WRF_name    = 'PB'
  character(len=25) :: t_WRF_name     = 'T'
  character(len=25) :: qv_WRF_name    = 'QVAPOR'
  character(len=25) :: pblh_WRF_name  = 'PBLH'
  character(len=25) :: ust_WRF_name   = 'UST'
  character(len=25) :: rmol_WRF_name  = 'RMOL'
  character(len=25) :: znt_WRF_name   = 'ZNT'
  character(len=25) :: v10_WRF_name   = 'V10'
  character(len=25) :: u10_WRF_name   = 'U10'
  character(len=25) :: hfx_WRF_name   = 'HFX'
  character(len=25) :: smoi_WRF_name  = 'SMOIS'
  character(len=25) :: ldu_WRF_name   = 'LU_INDEX'
  character(len=25) :: soil_WRF_name  = 'ISLTYP'
  character(len=25) :: vegfra_WRF_name= 'VEGFRA'
  character(len=25) :: prat_WRF_name  = 'PRATEC'
  !
  !*** WRF attribute names
  !
  character(len=25) :: WRF_map_proj_name = 'MAP_PROJ'
  character(len=25) :: WRF_cen_lon_name  = 'CEN_LON'
  character(len=25) :: WRF_truelat1_name = 'TRUELAT1'
  character(len=25) :: WRF_truelat2_name = 'TRUELAT2'
  !
  !*** Variables ID
  !
  integer(ip) :: phbID,phID,hgtID,lmaskID,uID,vID,wID,pID,pbID,tID,qvID,pblhID,ustID, &
                 rmolID,zntID,u10ID,v10ID,hfxID,lduID,soilID,vegfraID,smoiID,pratID
  !
  !*** WRF variable arrays
  !
  character(len=19), allocatable :: time_WRF_string(:)
  real(rp)         , allocatable :: time_WRF(:)
  real(rp)         , allocatable :: timesec_WRF(:)
  !
  real(rp)         , allocatable :: lat_WRF (:,:)
  real(rp)         , allocatable :: lon_WRF (:,:)
  real(rp)         , allocatable :: topg_WRF(:,:)
  real(rp)         , allocatable :: phb_WRF  (:,:,:)
  real(rp)         , allocatable :: z_WRF    (:,:,:)
  real(rp)         , allocatable :: zstag_WRF(:,:,:)
  !
  !*** WRF attribute values
  !
  integer(ip) :: WRF_map_proj
  real   (rp) :: WRF_cen_lon,WRF_truelat1,WRF_truelat2
  !
  !*** WRF time
  !
  integer(ip)   :: ibyr_WRF,ibmo_WRF,ibdy_WRF,ibhr_WRF,ibmi_WRF,ibse_WRF
  real   (rp)   :: dt_WRF
  !
  !**  Working arrays
  !
  real(rp), allocatable :: work_WRF(:,:,:),work_WRF1(:,:,:),work_WRF2(:,:,:)
  !
  !*** Interpolation stuff
  !
  integer(ip) :: nelem_WRF
  integer(ip) :: npoin_DAT
  !
  integer(ip), allocatable :: lnods_WRF(:,:)
  integer(ip), allocatable :: lelpo_WRF(:)
  !
  real   (rp), allocatable :: s_po_WRF(:)
  real   (rp), allocatable :: t_po_WRF(:)
  !
CONTAINS
  !
  !
  !
  subroutine WRF_getmem01
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
  end subroutine WRF_getmem01
  !
  !
  !
  subroutine WRF_getmem02
    !***************************************************************
    !*
    !*   Gets memory for time-independent variables
    !*
    !***************************************************************
    implicit none
    !
    allocate(time_WRF_string(nt_WRF))
    allocate(time_WRF(nt_WRF))
    allocate(timesec_WRF(nt_WRF))
    !
    !*** Variables
    !
    allocate(lat_WRF  (nx_WRF,ny_WRF))
    allocate(lon_WRF  (nx_WRF,ny_WRF))
    allocate(topg_WRF (nx_WRF,ny_WRF))
    allocate(phb_WRF  (nx_WRF,ny_WRF,nzstag_WRF))
    allocate(z_WRF    (nx_WRF,ny_WRF,nz_WRF    ))
    allocate(zstag_WRF(nx_WRF,ny_WRF,nzstag_WRF))
    !
    !*** Interpolation
    !
    allocate(lnods_WRF(4,nelem_WRF))
    allocate(lelpo_WRF(npoin_DAT))
    allocate(s_po_WRF(npoin_DAT))
    allocate(t_po_WRF(npoin_DAT))
    !
  end subroutine WRF_getmem02
  !
  !
  !
END MODULE WRF_nc
