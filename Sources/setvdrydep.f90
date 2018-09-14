  subroutine setvdrydep
  !***********************************************************************
  !*
  !*    Computes the dry deposition velocity according to Feng:
  !*    J.Feng, A size-resolved model and a four-mode parameterization
  !*    of dry deposition of atmospheric aerosols, JGR, v113, 2008.
  !*
  !*    Computation is done only for particles smaller than 100 mic
  !*    (the limit of the aerosol Giant mode) and for the second z-layer
  !*    (first is ground, where w=0)
  !*
  !*    The Feng model discriminates between 6 land-use categories:
  !*
  !*    1 'Urban'
  !*    2 'Remote continental'
  !*    3 'Desert'
  !*    4 'Polar'
  !*    5 'Marine'
  !*    6 'Rural'
  !*
  !*    whereas the DBS is based on the USGS 24-categories. We assume
  !*    the following simplifications:
  !*
  !*               USGS                                   Equivalent
  !*    ------------------------------------------------------------
  !*    1 'Urban and Built-Up Land'                         1
  !*    2 'Dryland Cropland and Pasture'                    6
  !*    3 'Irrigated Cropland and Pasture'                  6
  !*    4 'Mixed Dryland/Irrigated Cropland and Pasture'    6
  !*    5 'Cropland/Grassland Mosaic'                       6
  !*    6 'Cropland/Woodland Mosaic'                        6
  !*    7 'Grassland'                                       6
  !*    8 'Shrubland'                                       6
  !*    9 'Mixed Shrubland/Grassland'                       6
  !*   10 'Savanna'                                         3
  !*   11 'Deciduous Broadleaf Forest'                      6
  !*   12 'Deciduous Needleleaf Forest'                     6
  !*   13 'Evergreen Broadleaf Forest'                      6
  !*   14 'Evergreen Needleleaf Forest'                     6
  !*   15 'Mixed Forest'                                    6
  !*   16 'Water Bodies'                                    5
  !*   17 'Herbaceous Wetland'                              6
  !*   18 'Wooded Wetland'                                  6
  !*   19 'Barren or Sparsely Vegetated'                    3
  !*   20 'Herbaceous Tundra'                               2
  !*   21 'Wooded Tundra'                                   2
  !*   22 'Mixed Tundra'                                    2
  !*   23 'Bare Ground Tundra'                              2
  !*   24 'Snow or Ice'                                     4

  !***********************************************************************
  use KindType
  use Master
  use Numeric
  implicit none
  !
  integer(ip), parameter :: nland = 6
  integer(ip), parameter :: nmode = 4
  real(rp)   , parameter :: vkarm =  0.4_rp        ! von karman constant
  real(rp)   , parameter :: visc0 = 1.827e-5_rp    ! reference viscosity
  real(rp)   , parameter :: Tair0 = 291.15_rp      ! reference temperature
  !
  logical    , save      :: have_to_work
  integer(ip), save      :: ipass = 0
  integer(ip), save      :: imode(ncmax)
  real(rp)   , save      :: a(nland,nmode),b(nland,nmode)
  !
  integer(ip) :: ix,iy,ic,iland
  real   (rp) :: T, ra, visc, Re, Vd1, Vd
  real   (rp) :: zeta, zeta0, eta, eta0
  !
  !***  First time initialize values
  !
  if(ipass == 0) then
     ipass = 1
     !
     !***     a(iland,imode) and b(iland,imode)
     !
     a(1,1) = 0.0048_rp      ! Mode 1
     a(2,1) = 0.0037_rp
     a(3,1) = 0.0042_rp
     a(4,1) = 0.0032_rp
     a(5,1) = 0.0043_rp
     a(6,1) = 0.0045_rp
     b(1,1) = 1.0_rp
     b(2,1) = 1.0_rp
     b(3,1) = 1.0_rp
     b(4,1) = 1.0_rp
     b(5,1) = 1.0_rp
     b(6,1) = 1.0_rp
     !
     a(1,2) = 0.0315_rp     ! Mode 2
     a(2,2) = 0.0120_rp
     a(3,2) = 0.2928_rp
     a(4,2) = 0.1201_rp
     a(5,2) = 0.1337_rp
     a(6,2) = 0.0925_rp
     b(1,2) = 2.7925_rp
     b(2,2) = 2.2413_rp
     b(3,2) = 3.8581_rp
     b(4,2) = 3.4407_rp
     b(5,2) = 3.5456_rp
     b(6,2) = 3.2920_rp
     !
     a(1,3) = 1.2891_rp     ! Mode 3
     a(2,3) = 1.3977_rp
     a(3,3) = 1.3970_rp
     a(4,3) = 1.1838_rp
     a(5,3) = 1.2834_rp
     a(6,3) = 1.2654_rp
     b(1,3) = 2.6878_rp
     b(2,3) = 2.5838_rp
     b(3,3) = 2.5580_rp
     b(4,3) = 2.8033_rp
     b(5,3) = 2.7157_rp
     b(6,3) = 2.7227_rp
     !
     a(1,4) = 1.0338_rp     ! Mode 4
     a(2,4) = 1.0707_rp
     a(3,4) = 0.9155_rp
     a(4,4) = 1.0096_rp
     a(5,4) = 1.1595_rp
     a(6,4) = 1.0891_rp
     b(1,4) = 1.2644_rp
     b(2,4) = 1.3247_rp
     b(3,4) = 1.0364_rp
     b(4,4) = 1.2069_rp
     b(5,4) = 1.4863_rp
     b(6,4) = 1.4240_rp
     !
     !***     Loop over classes to averiguate the mode (including particles)
     !
     have_to_work = .false.
     do ic = 1,nc
        if(diam(ic) > 100e-6_rp) then  ! not considered
           imode(ic) = 0
        else if(diam(ic) > 10e-6_rp) then    ! > 10  mic   giant mode
           imode(ic) = 4
           have_to_work = .true.
        else if(diam(ic) > 2.5e-6_rp) then   ! > 2.5 mic   coarse mode
           imode(ic) = 3
           have_to_work = .true.
        else if(diam(ic) > 0.1e-6_rp) then   ! > 0.1 mic   accumulation mode
           imode(ic) = 2
           have_to_work = .true.
        else                               ! nucleii mode
           imode(ic) = 1
           have_to_work = .true.
        end if
     end do
     !
  end if   ! ipass
  !
  !***   have_to_work ?
  !
  if(.not.have_to_work) return
  !
  !***   Loop over points
  !
  do iy = 1,ny
     do ix = 1,nx
        !
        !***      Get the land coefficient
        !
        if(INT(LDU(ix,iy)) == 1) then
           iland = 1
        else if(INT(LDU(ix,iy)) == 2) then
           iland = 6
        else if(INT(LDU(ix,iy)) == 3) then
           iland = 6
        else if(INT(LDU(ix,iy)) == 4) then
           iland = 6
        else if(INT(LDU(ix,iy)) == 5) then
           iland = 6
        else if(INT(LDU(ix,iy)) == 6) then
           iland = 6
        else if(INT(LDU(ix,iy)) == 7) then
           iland = 6
        else if(INT(LDU(ix,iy)) == 8) then
           iland = 6
        else if(INT(LDU(ix,iy)) == 9) then
           iland = 6
        else if(INT(LDU(ix,iy)) == 10) then
           iland = 3
        else if(INT(LDU(ix,iy)) == 11) then
           iland = 6
        else if(INT(LDU(ix,iy)) == 12) then
           iland = 6
        else if(INT(LDU(ix,iy)) == 13) then
           iland = 6
        else if(INT(LDU(ix,iy)) == 14) then
           iland = 6
        else if(INT(LDU(ix,iy)) == 15) then
           iland = 6
        else if(INT(LDU(ix,iy)) == 16) then
           iland = 5
        else if(INT(LDU(ix,iy)) == 17) then
           iland = 6
        else if(INT(LDU(ix,iy)) == 18) then
           iland = 6
        else if(INT(LDU(ix,iy)) == 19) then
           iland = 3
        else if(INT(LDU(ix,iy)) == 20) then
           iland = 2
        else if(INT(LDU(ix,iy)) == 21) then
           iland = 2
        else if(INT(LDU(ix,iy)) == 22) then
           iland = 2
        else if(INT(LDU(ix,iy)) == 23) then
           iland = 2
        else if(INT(LDU(ix,iy)) == 24) then
           iland = 4
        else
           iland = 6
        end if
        !
        zeta  = rl(ix,iy)*zlayer(1)
        zeta0 = rl(ix,iy)*znt(ix,iy)
        !
        if(ustar(ix,iy) > 0.0_rp) then
           !neutral
           ra = LOG(zlayer(1)/znt(ix,iy)) 
           if(zeta>0) then
             !stable
             ra = ra + 4.7_rp * (zeta-zeta0)
           else if(zeta<0) then
             !unstable
             eta  = (1.0_rp-15.0_rp*zeta )**0.25
             eta0 = (1.0_rp-15.0_rp*zeta0)**0.25
             ra = ra + LOG((eta0**2+1)*(eta0+1)**2 / ((eta**2+1)*(eta+1)**2))
             ra = ra + 2*(ATAN(eta)-ATAN(eta0))
           end if
           ra = ra / (ustar(ix,iy)*vkarm)
        else
           ra = 1e20_rp   ! prevent overflow
        end if
        !
        !***      Air viscosity (Sutherland's law)
        !
        T    = tempe(ix,iy,1)
        visc = visc0*((Tair0+120.0_rp)/(T+120.0_rp))*((T/Tair0)**1.5_rp)
        !
        !***      Re*
        !
        Re = ustar(ix,iy)*znt(ix,iy)*rho(ix,iy,1)/visc
        if(Re>40300.0_rp) Re = 40300.0_rp
        !
        !***      Fist term of deposition velocity
        !
        Vd1 = -0.5_rp*((Re-40300.0_rp)/15330.0_rp)*((Re-40300.0_rp)/15330.0_rp)
        Vd1 = 0.0226_rp*ustar(ix,iy)*exp(Vd1)
        !
        !***      Loop over classes
        !
        do ic = 1,nc
           if(imode(ic) > 0) then
              Vd = Vd1 + a(iland,imode(ic))*(ustar(ix,iy)**b(iland,imode(ic)))
              Vd = ra + 1.0_rp/Vd
              vdep(ix,iy,ic) = -1.0_rp/Vd   ! sign criteria
           end if
        end do
        !
     end do   !  ix = 1,nx
  end do      !  iy = 1,ny
  !
  return
end subroutine setvdrydep
