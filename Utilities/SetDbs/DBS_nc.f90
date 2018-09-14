!***************************************************************
!*
!*    Module for WRF dbs operations
!*
!***************************************************************
MODULE DBS_nc
  use KindType, ONLY :  ip,rp
  implicit none
  save
  !
  !*** Database
  !
  integer(ip) :: ncID
  !
  !*** Attributes
  !
  character(len=25 ) :: attr_units
  character(len=100) :: attr_title,attr_desc
  !
  !*** DBS dimension (=coordinate variables) names and ID
  !
  character(len=25) :: lon_DBS_name   = 'lon'
  character(len=25) :: lat_DBS_name   = 'lat'
  character(len=25) :: x_DBS_name     = 'x'
  character(len=25) :: y_DBS_name     = 'y'
  character(len=25) :: alt_DBS_name   = 'alt'
  character(len=25) :: time_DBS_name  = 'time'
  !
  integer(ip)  :: nx_DBS_ID,ny_DBS_ID,nz_DBS_ID,nt_DBS_ID
  integer(ip)  :: lon_DBS_ID
  integer(ip)  :: lat_DBS_ID
  integer(ip)  :: x_DBS_ID
  integer(ip)  :: y_DBS_ID
  integer(ip)  :: alt_DBS_ID
  integer(ip)  :: time_DBS_ID
  !
  !*** DBS variable names and ID
  !
  character(len=25) :: xlon_DBS_name = 'XLON'
  character(len=25) :: xlat_DBS_name = 'XLAT'
  character(len=25) :: xutm_DBS_name = 'UTM_X'
  character(len=25) :: yutm_DBS_name = 'UTM_Y'
  character(len=25) :: topg_DBS_name = 'TOPG'
  character(len=25) :: u_DBS_name    = 'U'
  character(len=25) :: v_DBS_name    = 'V'
  character(len=25) :: w_DBS_name    = 'W'
  character(len=25) :: T_DBS_name    = 'T'
  character(len=25) :: Tp_DBS_name   = 'TP'
  character(len=25) :: p_DBS_name    = 'P'
  character(len=25) :: qv_DBS_name   = 'QV'
  character(len=25) :: ro_DBS_name   = 'RHO'
  character(len=25) :: pblh_DBS_name = 'PBLH'
  character(len=25) :: ust_DBS_name  = 'UST'
  character(len=25) :: smoi_DBS_name = 'SMOI'
  character(len=25) :: rmol_DBS_name = 'RMOL'
  character(len=25) :: znt_DBS_name  = 'ZNT'
  character(len=25) :: spd10_DBS_name= 'SPD10'
  character(len=25) :: L_DBS_name    = 'L'
  character(len=25) :: LDU_DBS_name  = 'LDU'
  character(len=25) :: SOIL_DBS_name = 'SOIL'
  character(len=25) :: VEGFRA_DBS_name = 'VEGFRA'
  character(len=25) :: prat_DBS_name = 'PRATE'
  !
  integer(ip)       :: xlon_DBS_ID
  integer(ip)       :: xlat_DBS_ID
  integer(ip)       :: xutm_DBS_ID
  integer(ip)       :: yutm_DBS_ID
  integer(ip)       :: topg_DBS_ID
  integer(ip)       :: u_DBS_ID
  integer(ip)       :: v_DBS_ID
  integer(ip)       :: w_DBS_ID
  integer(ip)       :: T_DBS_ID
  integer(ip)       :: Tp_DBS_ID
  integer(ip)       :: p_DBS_ID
  integer(ip)       :: qv_DBS_ID
  integer(ip)       :: ro_DBS_ID
  integer(ip)       :: pblh_DBS_ID
  integer(ip)       :: ust_DBS_ID
  integer(ip)       :: smoi_DBS_ID
  integer(ip)       :: rmol_DBS_ID
  integer(ip)       :: znt_DBS_ID
  integer(ip)       :: spd10_DBS_ID
  integer(ip)       :: L_DBS_ID
  integer(ip)       :: LDU_DBS_ID
  integer(ip)       :: SOIL_DBS_ID
  integer(ip)       :: VEGFRA_DBS_ID
  integer(ip)       :: prat_DBS_ID
  !
  !*** work array
  !
  real(rp), allocatable :: work_DBS(:)
  !
END MODULE DBS_nc
