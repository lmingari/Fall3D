subroutine get_par_ABL(z0,z1,Tz0,Tz1,Pz0,Pz1,u,v,ustar,rmonin)
  !*****************************************************************************
  !*
  !*     This routine estimates some parameters of the ABL
  !*
  !*     INPUTS:
  !*      z0       Reference level 1
  !*      z1       Reference level 2 (z1 > z0)
  !*      Tz0      Temperature (K) at z0
  !*      Tz1      Temperature (K) at z1
  !*      Pz0      Pressure    (Pa) at z0
  !*      Pz1      Pressure    (Pa) at z1
  !*      u        u-velocity  (at z1)
  !*      v        v-velocity  (at z1)
  !*     OUTPUTS:
  !*      ustar    Friction velocity
  !*      rmonin   Monin-Obukhov lenght
  !*
  !*****************************************************************************
  use KindType
  implicit none
  !
  real(rp) :: z0,z1,Tz0,Tz1,Pz0,Pz1,u,v
  real(rp) :: ustar,rmonin
  !
  real(rp), parameter :: g = 9.81_rp
  real(rp), parameter :: k = 0.4_rp
  real(rp)            :: Thz0,Thz1,umod,Ri,Gm,Gh,thstar
  !
  !***   Potential temperature Theta=T*(1bar/P))**(R/cp)
  !***   R  =  287  J/kg K   Air specific gas  constant
  !***   cp = 1006  J/kg K   Air specific heat capacity
  !
  Thz0 = Tz0*(1.01d5/Pz0)**(0.285)
  Thz1 = Tz1*(1.01d5/Pz1)**(0.285)
  !
  !***   Velocity
  !
  umod = max(sqrt(u*u+v*v),1e-4_rp)
  !
  !***   Bulk Richardson number
  !
  Ri = g*(z1-z0)*(Thz1-Thz0)/(umod*umod*0.5_rp*(Tz0+Tz1))
  !
  !***   Stable/Unstable ABL
  !
  if(Ri.ge.0.0_rp) then
     !                             Stable
     Gm = 1.0_rp + 4.7_rp*Ri
     Gm = 1.0_rp/(Gm*Gm)
     Gh = Gm
  else
     !                             Unstable
     Gm = log(z1/z0)
     Gm = 1.0_rp + (70.0_rp*k*k*sqrt(abs(Ri)*z1/z0))/(Gm*Gm)
     Gh = 1.0_rp + (50.0_rp*k*k*sqrt(abs(Ri)*z1/z0))/(Gm*Gm)
     Gm = 1.0_rp - ((9.4_rp*Ri)/Gm)
     Gh = 1.0_rp - ((9.4_rp*Ri)/Gh)
  end if
  !
  !***   ustar
  !
  ustar = log(z1/z0)
  ustar = k*umod*sqrt(Gm)/ustar
  !
  !***   thstar (Pr=1 assumed)
  !
  thstar = log(z1/z0)
  thstar = k*k*umod*(Thz1-Thz0)*Gh/(ustar*thstar*thstar)
  thstar = max(thstar,1e-4_rp)
  !
  !***   Monin-Obukhov
  !
  rmonin = ustar*ustar*0.5_rp*(Thz0+Thz1)/(k*g*thstar)
  !
  return
end subroutine get_par_ABL
