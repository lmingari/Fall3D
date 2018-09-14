!***************************************************************
!*
!*    Module for Plume equations parameters
!*
!***************************************************************
MODULE Plume
  use KindType
  implicit none
  save
  !
  real(rp), parameter :: Cv = 1996.0   ! vapour
  real(rp), parameter :: Cp = 1250.0   ! particles
  real(rp), parameter :: Ca = 1005.0   ! air
  !
  real(rp), parameter :: zmin_wind = 100.0    ! Ignore wind entrainment below this zvalue (low jet region)
  real(rp), parameter :: zjet  = 500.0        ! Height of the jet region (above vent)
  real(rp), parameter :: c_umbrella = 1.320  ! Thichness of umbrella relative to Hb (>1)
  !
  real(rp)  :: xi    = 0.1      ! Bursik factor.
  real(rp)  :: Cbo   = Cp          ! bulk at z=0
  real(rp)  :: Cb                  ! bulk

  !
END MODULE Plume
