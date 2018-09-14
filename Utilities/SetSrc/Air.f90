!***************************************************************
!*
!*		Module for Air operations
!*
!***************************************************************
MODULE Air
  use KindType
  IMPLICIT NONE
  SAVE
  !
  real(rp), parameter :: ztropo  = 1d4          ! Tropoause height
  real(rp), parameter :: rhoair0 = 1.2255    ! standard at sea level
  !
CONTAINS
  !
  !
  !
  real(rp) function Ta(z)
    !**************************************************************
    !*
    !*   Gets air temperature (K) at z (a.s.l.)
    !*
    !**************************************************************
    use Master
    implicit none
    real(rp) :: z
    !
    logical     :: go_on
    integer(ip) :: iz1,iz2
    real(rp)    :: f1,f2
    !
    if(standard_atmosphere) then
       !
       !***  1.Standard atmosphere
       !
       Ta  = Tair(1) - (z-Z0-dZ0-zmodel(1))*6.5e-3
       !
    else
       !
       !***  2.Values interpolated from Tair(nz), extracted from dbs
       !
       !***    Value of z is below the first layer
       !
       if(z.lt.(Z0+dZ0+zmodel(1))) then
          Ta = Tair(1)
          return
       end if
       !
       !***    Gets the position indexes iz1,iz2  and the weights
       !
       iz1 = 0
       go_on = .true.
       do while(go_on)
          iz1 = iz1 + 1
          iz2 = iz1 + 1
          if(z.ge.(Z0+dZ0+zmodel(iz1)).and.z.lt.(Z0+dZ0+zmodel(iz2))) then
             go_on = .false.
             f1 = 1.0-((z-Z0-dz0-zmodel(iz1))/(zmodel(iz2)-zmodel(iz1)))
             f2 = 1.0-f1
          else if(iz2.eq.nz) then
             call runend('Source position not found in function Ta')
          end if
       end do
       !
       Ta = f1*Tair(iz1)+f2*Tair(iz2)
       !
    end if
    return
  end function Ta
  !
  !
  !
  real(rp) function Va(z)
    !**************************************************************
    !*
    !*    Interpolates air velocity (m/s) at z a.s.l.
    !*    from Vair
    !*
    !**************************************************************
    use Master
    use Plume
    implicit none
    real(rp) :: z
    !
    logical     :: go_on
    integer(ip) :: iz1,iz2
    real(rp)    :: f1,f2
    !
    !***  Value of z is below the first layer
    !
    if(z.lt.(Z0+dZ0+zmodel(1))) then
       Va = Vair(1)
       return
    end if
    !
    !***  Gets the position indexes iz1,iz2  and the weights
    !
    iz1 = 0
    go_on = .true.
    do while(go_on)
       iz1 = iz1 + 1
       iz2 = iz1 + 1
       if(z.ge.(Z0+dZ0+zmodel(iz1)).and.z.lt.(Z0+dZ0+zmodel(iz2))) then
          go_on = .false.
          f1 = 1.0-((z-Z0-dZ0-zmodel(iz1))/(zmodel(iz2)-zmodel(iz1)))
          f2 = 1.0-f1
       else if(iz2.eq.nz) then
          call runend('Function Va : Source position not found')
       end if
    end do
    !
    Va = f1*Vair(iz1)+f2*Vair(iz2)
    !
    !***  Sets wind velocity to zero below zmin_wind. This is to avoid
    !***  entrainment at the very low jet region
    !
    if(z.lt.zmin_wind) Va = 0.0
    !
    return
  end function Va
  !
  !
  !
  real(rp) function dVadz(z)
    !**************************************************************
    !*
    !*    Interpolates air velocity gradient
    !*    from Vair
    !*
    !**************************************************************
    use Master
    use Plume
    implicit none
    real(rp) :: z
    !
    logical     :: go_on
    integer(ip) :: iz1,iz2
    !
    !***  Value of z is below the first layer
    !
    if(z.lt.(Z0+dZ0+zmodel(1))) then
       dVadz = Vair(1)/zmodel(1)
       return
    end if
    !
    !***  Gets the position indexes iz1,iz2  and the weights
    !
    iz1 = 0
    go_on = .true.
    do while(go_on)
       iz1 = iz1 + 1
       iz2 = iz1 + 1
       if(z.ge.(Z0+dZ0+zmodel(iz1)).and.z.lt.(Z0+dZ0+zmodel(iz2))) then
          go_on = .false.
          dVadz = (Vair(iz2)-Vair(iz1))/(zmodel(iz2)-zmodel(iz1))
       else if(iz2.eq.nz) then
          call runend('Function dVadz : Source position not found')
       end if
    end do
    !
    return
  end function dVadz
  !
  !
  !
  real(rp) function rhoa(z)
    !****************************************************************
    !*
    !*    Computes air density at z a.s.l.
    !*
    !****************************************************************
    use Master
    use Plume
    implicit none
    logical     :: go_on
    integer(ip) :: iz1,iz2
    real(rp)    :: f1,f2,z
    !
    if(standard_atmosphere) then
       !
       !***  1.Standard atmosphere
       !
       if(z.lt.ztropo) then
          rhoa = rhoair0*((1.0-2.26e-5*z     )**(4.255))
       else
          rhoa = rhoair0*((1.0-2.26e-5*ztropo)**(4.255))
          rhoa = rhoa*exp( (ztropo-z)/(29.2745*Ta(z)) )
       end if
       !
    else
       !
       !***  2.Values interpolated from Rair(nz), extracted from dbs
       !
       !***    Value of z is below the first layer
       !
       if(z.lt.(Z0+dZ0+zmodel(1))) then
          rhoa = Rair(1)
          return
       end if
       !
       !***    Gets the position indexes iz1,iz2  and the weights
       !
       iz1 = 0
       go_on = .true.
       do while(go_on)
          iz1 = iz1 + 1
          iz2 = iz1 + 1
          if(z.ge.(Z0+dZ0+zmodel(iz1)).and.z.lt.(Z0+dZ0+zmodel(iz2))) then
             go_on = .false.
             f1 = 1.0-((z-Z0-dz0-zmodel(iz1))/(zmodel(iz2)-zmodel(iz1)))
             f2 = 1.0-f1
          else if(iz2.eq.nz) then
             call runend('Source position not found in function rhoa')
          end if
       end do
       !
       rhoa = f1*Rair(iz1)+f2*Rair(iz2)
       !
    endif
    !
    return
  end function rhoa
  !
  !
  !
  real(rp) function mua(z)
    !************************************************************
    !*
    !*    Computes air viscosity at z for the standard atmosphere
    !*
    !************************************************************
    implicit none
    real(rp) :: z,mua0,Tair0
    !
    mua0  = 1.827d-5  ! reference viscosity
    Tair0 = 291.15    ! reference temperature
    !
    mua = mua0*((Tair0+120.0_rp)/(Ta(z)+120.0_rp))*((Ta(z)/Tair0)**1.5_rp)   ! Sutherland's law
    !
    return
  end function mua
  !
  !
  !
END MODULE Air
