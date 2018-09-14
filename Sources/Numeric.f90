!***************************************************************
!*
!*		Module for numeric operations
!*
!***************************************************************
MODULE Numeric
  use KindType
  IMPLICIT NONE
  SAVE
  !
  !***  Constants
  !
  real(rp) :: pi     = 3.14159265358979323846_rp   ! pi
  real(rp) :: Rearth = 6356d3                      ! Earth's radius
  !
  !***  Integer flags
  !
  integer(ip) :: modkv,modkh,modv,lmeth
  !
  !***  Atmophere
  !
  real(rp) :: rkh0   = 0.0_rp
  real(rp) :: rkv0   = 0.0_rp
  real(rp) :: rhoa0  = 1.225_rp    ! Atmospheric density at sea level
  real(rp) :: ztropo = 15d3        ! Standard value for tropopause
  real(rp) :: Csh    = 0.2275_rp   ! Parameter in RAMS formula
  !
  !***  Mesh related variables
  !
  integer(ip) :: nx,ny,nz
  real   (rp) :: lonmin,lonmax,latmin,latmax
  real   (rp) :: xmin  ,xmax  ,ymin  ,ymax
  real   (rp) :: xorigr,yorigr,dx,dy,tpgmax,tpgmin,ztop
  !
  !***  System of coordinates and transformations
  !
  character(len=10)  :: coord_sys
  real(rp)           :: dX1,dX2
  !
END MODULE Numeric
