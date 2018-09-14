!***************************************************************
!*
!*    Module for NCEP2 netcdf operations
!*
!***************************************************************
     MODULE NCEP2_nc
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
!*** NCEP2 dimension names
!
     character(len=25) :: nt_NCEP2_name = 'time'
     character(len=25) :: nx_NCEP2_name = 'lon'
     character(len=25) :: ny_NCEP2_name = 'lat'
     character(len=25) :: nz_NCEP2_name = 'pres'
!
     integer(ip) :: nt_NCEP2,it_NCEP2
     integer(ip) :: nx_NCEP2
     integer(ip) :: ny_NCEP2
     integer(ip) :: nz_NCEP2
!
!*** NCEP2 variable names
!
     character(len=25) :: time_NCEP2_name  = 'time'
     character(len=25) :: p_NCEP2_name     = 'pres'
     character(len=25) :: fis_NCEP2_name   = 'HGT:sfc'
     character(len=25) :: fi_NCEP2_name    = 'HGT'
     character(len=25) :: T_NCEP2_name     = 'TMP'
     character(len=25) :: Q_NCEP2_name     = 'RH'
     character(len=25) :: u_NCEP2_name     = 'UGRD'
     character(len=25) :: v_NCEP2_name     = 'VGRD'
     character(len=25) :: w_NCEP2_name     = 'VVEL'
!
!*** NCEP2 attribute names
!
     character(len=25) :: ibyr_NCEP2_name = 'YEAR'
     character(len=25) :: ibmo_NCEP2_name = 'MONTH'
     character(len=25) :: ibdy_NCEP2_name = 'DAY'
     character(len=25) :: ibhr_NCEP2_name = 'HOUR'
     character(len=25) :: dt_NCEP2_name   = 'TIME_INCR'
     character(len=25) :: lonmin_NCEP2_name = 'LONMIN'
     character(len=25) :: lonmax_NCEP2_name = 'LONMAX'
     character(len=25) :: latmin_NCEP2_name = 'LATMIN'
     character(len=25) :: latmax_NCEP2_name = 'LATMAX'
!
!*** NCEP2 variables ID
!
     integer(ip) :: fiID,fisID,tID,qID,uID,vID,wID,lmaskID,pblhID
!
!*** NCEP2 attribute values
!
     integer(ip) :: ibyr_NCEP2,ibmo_NCEP2,ibdy_NCEP2,ibhr_NCEP2,ibmi_NCEP2,ibse_NCEP2
     integer(ip) :: dt_NCEP2
     real   (rp) :: lonmin_NCEP2,lonmax_NCEP2,latmin_NCEP2,latmax_NCEP2
!
!*** NCEP2 variable arrays
!
     real(rp), allocatable :: time_NCEP2(:)
     real(rp), allocatable :: timesec_NCEP2(:)
!
     real(rp), allocatable :: p_NCEP2(:)
     real(rp), allocatable :: lat_NCEP2 (:,:)
     real(rp), allocatable :: lon_NCEP2 (:,:)
     real(rp), allocatable :: topg_NCEP2(:,:)
     real(rp), allocatable :: z_NCEP2   (:,:,:)
!
!**  Working arrays
!
     real(rp), allocatable :: work_NCEP21(:,:,:),work_NCEP22(:,:,:)
!
!*** Interpolation stuff
!
     integer(ip) :: nelem_NCEP2
     integer(ip) :: npoin_DAT
!
     integer(ip), allocatable :: lnods_NCEP2(:,:)
     integer(ip), allocatable :: lelpo_NCEP2(:)
!
     real   (rp), allocatable :: s_po_NCEP2(:)
     real   (rp), allocatable :: t_po_NCEP2(:)
!
!
!
     CONTAINS
!
     subroutine NCEP2_getmem01
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
     end subroutine NCEP2_getmem01
!
!
!
     subroutine NCEP2_getmem02
!***************************************************************
!*
!*   Gets memory for time-independent variables
!*
!***************************************************************
     implicit none
!
     allocate(time_NCEP2(nt_NCEP2))
     allocate(timesec_NCEP2(nt_NCEP2))
!
!*** Variables
!
     allocate(p_NCEP2   (nz_NCEP2))
     allocate(lat_NCEP2 (nx_NCEP2,ny_NCEP2))
     allocate(lon_NCEP2 (nx_NCEP2,ny_NCEP2))
     allocate(topg_NCEP2(nx_NCEP2,ny_NCEP2))
     allocate(z_NCEP2   (nx_NCEP2,ny_NCEP2,nz_NCEP2))
!
!*** Interpolation
!
     allocate(lnods_NCEP2(4,nelem_NCEP2))
     allocate(lelpo_NCEP2(npoin_DAT))
     allocate(s_po_NCEP2(npoin_DAT))
     allocate(t_po_NCEP2(npoin_DAT))
!
     end subroutine NCEP2_getmem02
!
     END MODULE NCEP2_nc
