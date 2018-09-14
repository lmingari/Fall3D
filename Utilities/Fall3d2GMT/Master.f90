!***************************************************************
!*
!*		Module for master operations
!*
!***************************************************************
MODULE Master
  use KindType
  IMPLICIT NONE
  SAVE
  !
  !***  Version of input file
  !
  real(rp),parameter :: version = 7.0
  !
  !***  Plot definition
  !
  real(rp)    :: lonmin,lonmax,latmin,latmax
  !
  logical     :: pp_topog = .false.
  logical     :: pp_load0 = .false.
  logical     :: pp_loadw = .false.
  logical     :: pp_loadc = .false.
  logical     :: pp_wetc  = .false.
  logical     :: pp_thick = .false.
  logical     :: pp_concg = .false.
  logical     :: pp_PMxxg = .false.
  logical     :: pp_cumul = .false.
  logical     :: pp_PMxxc = .false.
  logical     :: pp_fl    = .false.
  logical     :: pp_aod05 = .false.
  !
  !***  Plot related variables
  !
  character(len=s_name) :: varname,classname(100),classnamew(100)
  integer(ip) :: nclass,nclassw
  !
  integer(ip), parameter  :: mxnval = 100  ! max number of contour values
  integer(ip) :: &
       nval, &
       nval_topog, &
       nval_load0, &
       nval_loadw, &
       nval_loadc, &
       nval_wetc,  &
       nval_thick, &
       nval_concg, &
       nval_PMxxg, &
       nval_cumul, &
       nval_PMxxc, &
       nval_fl,    &
       nval_aod05
  !
  character(len=12) :: &
       fact, &
       fact_topog, &
       fact_load0, &
       fact_loadw, &
       fact_loadc, &
       fact_wetc,  &
       fact_thick, &
       fact_concg, &
       fact_PMxxg, &
       fact_cumul, &
       fact_PMxxc, &
       fact_fl,    &
       fact_aod05
  !
  real(rp),dimension(mxnval) :: &
       cval, &
       cval_topog, &
       cval_load0, &
       cval_loadw, &
       cval_loadc, &
       cval_wetc,  &
       cval_thick, &
       cval_concg, &
       cval_PMxxg, &
       cval_cumul, &
       cval_PMxxc, &
       cval_fl,    &
       cval_aod05
  !
  character(len=s_name) :: &
       unit, &
       unit_topog, &
       unit_load0, &
       unit_loadw, &
       unit_loadc, &
       unit_wetc,  &
       unit_thick, &
       unit_concg, &
       unit_PMxxg, &
       unit_cumul, &
       unit_PMxxc, &
       unit_fl,    &
       unit_aod05
  !
  !***  Time
  !
  integer(ip)             :: ifl,iyr,imo,idy,ihr,imi,irunb,irune
  integer(ip)             :: nt
  real   (rp)      ,allocatable :: times(:)
  character(len=18),allocatable :: dates(:)
  !
  !*** Attribute names
  !
  character(len=s_name) :: nc_title_name  = 'TITLE'
  character(len=s_name) :: nc_coord_name  = 'COORDINATES'
  character(len=s_name) :: nc_lonmin_name = 'LONMIN'
  character(len=s_name) :: nc_lonmax_name = 'LONMAX'
  character(len=s_name) :: nc_latmin_name = 'LATMIN'
  character(len=s_name) :: nc_latmax_name = 'LATMAX'
  character(len=s_name) :: nc_xmin_name   = 'XMIN'
  character(len=s_name) :: nc_xmax_name   = 'XMAX'
  character(len=s_name) :: nc_ymin_name   = 'YMIN'
  character(len=s_name) :: nc_ymax_name   = 'YMAX'
  character(len=s_name) :: nc_iyr_name    = 'YEAR'
  character(len=s_name) :: nc_imo_name    = 'MONTH'
  character(len=s_name) :: nc_idy_name    = 'DAY'
  character(len=s_name) :: nc_irunb_name  = 'RUN_START'
  character(len=s_name) :: nc_irune_name  = 'RUN_END'
  !
END MODULE Master
