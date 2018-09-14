   subroutine merwind
  !********************************************************
  !*
  !*   Computes MER depending on the wind profile
  !*
  !********************************************************
  use KindType
  use InpOut
  use Master
  implicit none
  !
  integer(ip) :: iz,jz
  real   (rp) :: v_mean,N_mean
  real   (rp) :: Cao,Co,To,alfa,beta,z1,gprime
  real   (rp) :: H1,V1,Ws
  !
  SELECT CASE(MER_vs_h)
  case('ESTIMATE-DEGRUYTER')
    !
    !*** Estimates MER as in Degruyter and Bonadonna (2012)
    !
    Cao     =  998.0_rp    ! specific heat capacity at constant pressure of dry air (J kg^-1 K^-1)
    Co      = 1250.0_rp    ! specific heat capacity at constant pressure of solids (J kg^-1 K^-1)
    To      = 1200.0_rp    ! initial plume temperature (K)
    alfa    = 0.1_rp       ! radial entrainment coefficient
    beta    = 0.5_rp       ! wind entrainment coefficient
    z1      = 2.8_rp       ! maximum non-dimensional height (Morton et al. 1956)
    !
    gprime = g*(Co*To-Cao*Tair0)/(Cao*Tair0)
    !
    v_mean = 0.0_rp
    N_mean = 0.0_rp
    jz     = 2
    do iz = 2,nz
       if(zmodel(iz-1).le.HPlume) then
           v_mean = v_mean + 0.5_rp*(Vair(iz-1)+Vair(iz))*(zmodel(iz)-zmodel(iz-1))
           N_mean = N_mean + 0.5_rp*(Nair(iz-1)+Nair(iz))*(zmodel(iz)-zmodel(iz-1))
           jz = iz
       end if
    end do
    v_mean = v_mean/(zmodel(jz)-zmodel(1))
    N_mean = N_mean/(zmodel(jz)-zmodel(1))
    N_mean = sqrt(N_mean)
    !
    M0 = (2.0_rp**(5.0_rp/2.0_rp))*alfa*alfa*N_mean*N_mean*N_mean*HPlume*HPlume*HPlume*HPlume/(z1*z1*z1*z1)
    M0 = M0 + beta*beta*N_mean*N_mean*v_mean*HPlume*HPlume*HPlume/6.0_rp
    M0 = pi*Rair0*M0/gprime
    !
  case('ESTIMATE-WOODHOUSE')
    !
    !*** Estimates MER as in Woodhouse et al. (2012)
    !
    N_mean = 0.0_rp
    jz     = 2
    do iz = 2,nz
       if(zmodel(iz-1).le.HPlume) then
           N_mean = N_mean + 0.5_rp*(Nair(iz-1)+Nair(iz))*(zmodel(iz)-zmodel(iz-1))
           jz = iz
       end if
    end do
    N_mean = N_mean/(zmodel(jz)-zmodel(1))
    N_mean = sqrt(N_mean)
    !
    V1 = Vair(jz)                 ! reference velocity
    H1 = zmodel(jz) - zmodel(1)   ! reference height
    !
    Ws = 1.44_rp*V1/(N_mean*H1)
    Ws = 0.318_rp*(1.0_rp+1.373_rp*Ws)/(1.0_rp+4.266_rp*Ws+0.3527_rp*Ws*Ws)
    M0 = (HPlume/1d3/Ws)**(1.0_rp/0.253_rp)
    !
  END SELECT
  !
  return
  end subroutine merwind
