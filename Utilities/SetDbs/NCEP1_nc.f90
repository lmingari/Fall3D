!***************************************************************
!*
!*    Module for NCEP1 netcdf operations
!*
!***************************************************************
     MODULE NCEP1_nc
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
!*** NCEP1 dimension names
!
     character(len=25) :: nt_NCEP1_name = 'time'
     character(len=25) :: nx_NCEP1_name = 'lon'
     character(len=25) :: ny_NCEP1_name = 'lat'
     character(len=25) :: nz_NCEP1_name = 'pres'
!
     integer(ip) :: nt_NCEP1,it_NCEP1
     integer(ip) :: nx_NCEP1
     integer(ip) :: ny_NCEP1
     integer(ip) :: nz_NCEP1
!
!*** NCEP1 variable names
!
     character(len=25) :: time_NCEP1_name  = 'time'
     character(len=25) :: p_NCEP1_name     = 'pres'
     character(len=25) :: fis_NCEP1_name   = 'HGT.sfc'
     character(len=25) :: fi_NCEP1_name    = 'HGT'
     character(len=25) :: Ps_NCEP1_name    = 'P.sfc'
     character(len=25) :: Ts_NCEP1_name    = 'TMP.sfc'
     character(len=25) :: T_NCEP1_name     = 'TMP'
     character(len=25) :: Q_NCEP1_name     = 'RH'
     character(len=25) :: us_NCEP1_name    = 'U.sfc'
     character(len=25) :: u_NCEP1_name     = 'UGRD'
     character(len=25) :: vs_NCEP1_name    = 'V.sfc'
     character(len=25) :: v_NCEP1_name     = 'VGRD'
     character(len=25) :: w_NCEP1_name     = 'VVEL'
!
!*** NCEP1 attribute names
!
     character(len=25) :: ibyr_NCEP1_name = 'YEAR'
     character(len=25) :: ibmo_NCEP1_name = 'MONTH'
     character(len=25) :: ibdy_NCEP1_name = 'DAY'
     character(len=25) :: ibhr_NCEP1_name = 'HOUR'
     character(len=25) :: dt_NCEP1_name   = 'TIME_INCR'
     character(len=25) :: lonmin_NCEP1_name = 'LONMIN'
     character(len=25) :: lonmax_NCEP1_name = 'LONMAX'
     character(len=25) :: latmin_NCEP1_name = 'LATMIN'
     character(len=25) :: latmax_NCEP1_name = 'LATMAX'
!
!*** NCEP1 variables ID
!
     integer(ip) :: fiID,fisID,psID,tsID,tID,qID,usID,uID,vsID,vID,wID,lmaskID,pblhID
!
!*** NCEP1 attribute values
!
     integer(ip) :: ibyr_NCEP1,ibmo_NCEP1,ibdy_NCEP1,ibhr_NCEP1,ibmi_NCEP1,ibse_NCEP1
     integer(ip) :: dt_NCEP1
     real   (rp) :: lonmin_NCEP1,lonmax_NCEP1,latmin_NCEP1,latmax_NCEP1
!
!*** NCEP1 variable arrays
!
     real(rp), allocatable :: time_NCEP1(:)
     real(rp), allocatable :: timesec_NCEP1(:)
!
     real(rp), allocatable :: p_NCEP1(:)
     real(rp), allocatable :: lat_NCEP1 (:,:)
     real(rp), allocatable :: lon_NCEP1 (:,:)
     real(rp), allocatable :: topg_NCEP1(:,:)
     real(rp), allocatable :: z_NCEP1   (:,:,:)
!
!**  Working arrays
!
     real(rp), allocatable :: work_NCEP11(:,:,:),work_NCEP12(:,:,:)
!
!*** Interpolation stuff
!
     integer(ip) :: nelem_NCEP1
     integer(ip) :: npoin_DAT
!
     integer(ip), allocatable :: lnods_NCEP1(:,:)
     integer(ip), allocatable :: lelpo_NCEP1(:)
!
     real   (rp), allocatable :: s_po_NCEP1(:)
     real   (rp), allocatable :: t_po_NCEP1(:)
!
!
!
     CONTAINS
!
     subroutine NCEP1_getmem01
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
     end subroutine NCEP1_getmem01
!
!
!
     subroutine NCEP1_getmem02
!***************************************************************
!*
!*   Gets memory for time-independent variables
!*
!***************************************************************
     implicit none
!
     allocate(time_NCEP1(nt_NCEP1))
     allocate(timesec_NCEP1(nt_NCEP1))
!
!*** Variables
!
     allocate(p_NCEP1   (nz_NCEP1))
     allocate(lat_NCEP1 (nx_NCEP1,ny_NCEP1))
     allocate(lon_NCEP1 (nx_NCEP1,ny_NCEP1))
     allocate(topg_NCEP1(nx_NCEP1,ny_NCEP1))
     allocate(z_NCEP1   (nx_NCEP1,ny_NCEP1,nz_NCEP1))
!
!*** Interpolation
!
     allocate(lnods_NCEP1(4,nelem_NCEP1))
     allocate(lelpo_NCEP1(npoin_DAT))
     allocate(s_po_NCEP1(npoin_DAT))
     allocate(t_po_NCEP1(npoin_DAT))
!
     end subroutine NCEP1_getmem02
!
     END MODULE NCEP1_nc
