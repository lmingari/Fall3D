subroutine setcut_FL(nx,ny,nz,zlayer,var2d,var3d,topg,p,fl)
  !**************************************************************
  !*
  !*    Does the cut at a given Flight Level.
  !*    A Flight Level (FL) is a standard nominal altitude of an aircraft,
  !*    in hundreds of feet. This altitude is calculated from a world-wide
  !*    fixed pressure datum of 1013.25 hPa, the average sea-level pressure,
  !*    and therefore is not necessarily the same as the aircraft's true altitude
  !*    either above mean sea level or above ground level.
  !*
  !**************************************************************
  use KindType
  implicit none
  integer(ip) :: nx,ny,nz
  real   (rp) :: fl
  real   (rp) :: zlayer(nz),var2d(nx,ny),topg(nx,ny), &
       var3d(nx,ny,nz),p(nx,ny,nz)
  !
  logical     :: go_on
  integer(ip) :: ix,iy,iz
  real   (rp) :: s,pnom,Tnom
  !
  !***  Get the nominal pressure and temperature at height fl (in m)
  !
  call standard_atm(pnom,Tnom,fl)
  !
  !***  Loop over points
  !
  do iy = 1,ny
     do ix = 1,nx
        !
        !***        Gets the layer pressures
        !
        if(pnom > p(ix,iy,1)) then    ! First level
           var2d(ix,iy) = 0.0_rp
           !
        else
           go_on = .true.
           iz = 0
           do while(go_on)
              iz = iz + 1
              if( (pnom <= p(ix,iy,iz  )) .and. &
                   (pnom > p(ix,iy,iz+1)) ) then
                 go_on = .false.
                 s = (pnom-p(ix,iy,iz))/(p(ix,iy,iz+1)-p(ix,iy,iz))
                 var2d(ix,iy) = var3d(ix,iy,iz)*(1.0_rp-s) + var3d(ix,iy,iz+1)*s
              end if
              if( (iz+1) == nz .and. go_on ) then
                 var2d(ix,iy) = 0.0_rp
                 go_on = .false.
              end if
           end do
        end if
        !
     end do
  end do
  !
  return
end subroutine setcut_FL
!
!
!
subroutine standard_atm(p,T,z)
  !*****************************************************************************
  !*
  !*    Returns p and T at height z (a.s.l.) for the ICAO standard atmosphere
  !*    Only three layers are used (valid up to 32km)
  !*
  !*****************************************************************************
  use KindType
  implicit none
  !
  real(rp) :: p,T,z
  !
  real(rp) :: goo = 9.80665_rp    ! standard gravity
  real(rp) :: Roo = 8.31432_rp    ! universal gas constant Nm / (mol K)
  real(rp) :: Moo = 0.0289644_rp  ! molar mass of air   (kg/mol)
  !
  real(rp) :: zoo =  0.0_rp       ! First layer  (Troposphere)
  real(rp) :: Too = 288.15_rp
  real(rp) :: poo = 101325_rp
  real(rp) :: Loo = -0.0065_rp
  !
  real(rp) :: zo1 = 11000.0_rp     ! Second layer (Stratosphere)
  real(rp) :: To1 = 216.65_rp
  real(rp) :: po1 = 22632_rp
  !
  real(rp) :: zo2 = 20000.0_rp     ! Third layer (Stratosphere)
  real(rp) :: To2 = 216.65_rp
  real(rp) :: po2 = 5474_rp
  real(rp) :: Lo2 = 0.001_rp
  !
  !***  computes
  !
  if(z < zo1) then
     T = Too + Loo*(z-zoo)
     P = poo*( (Too/T)**(goo*Moo/(Roo*Loo)) )
  else if(z < zo2) then
     T = To1
     P = Po1*exp(-goo*Moo*(z-zo1)/(Roo*To1))
  else
     T = To2 + Lo2*(z-zo2)
     P = po2*( (To2/T)**(goo*Moo/(Roo*Lo2)) )
  end if
  !
  return
end  subroutine standard_atm
