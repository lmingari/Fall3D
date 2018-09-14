subroutine addradialwind
  !*****************************************************************
  !
  !  This subroutine modifies the wind field adding a zero-divergence field.
  !  to account for gravity currents
  !
  ! Authors: Costa, Macedonio
  !
  !  TO CHECK:
  !   - The wind components vx,vy are corrected by H1m(i,j) inside
  !     the subroutine meteo after the call to this routine.
  !     For this reason, this correction is not implemented here (check).
  !
  !   - The optimal time step is recomputed in subroutine meteo, after
  !     the call to this subroutine (should be OK).
  !
  !   - This subroutine implements both LON-LAT and UTM coordinate models
  !     However, only LON-LAT was checked (Check for UTM to be done).
  !
  !   - The minimum radius of the gravity currente is set proportional
  !     to the cell size (see below)
  !*****************************************************************
  use KindType
  use InpOut,   only: nlst
  use parallel, only: mpime, root
  use Master
  use Numeric
  implicit none
  !
  integer(ip) :: kmin,kmax            ! Lower and upper levels of the cylinder (cell nums)
  integer(ip) :: i,j,k
  integer(ip) :: i_wind,j_wind,k_wind
  !
  real(rp) :: rvel                           ! Radial velocity (m/s)
  real(rp) :: deltax,deltay,radius           ! Distance form the vent
  real(rp) :: dist,store
  real(rp) :: wind_center_x, wind_center_y   ! Center of the radial wind
  real(rp) :: min_radius                     ! Minimum radius of the gravity current
  real(rp) :: max_radius                     ! Maximum radius of the gravity current
  real(rp) :: h_tot,h_umbr, hmin,hmax
  real(rp) :: v_wind,Tp,Tb
  real(rp) :: th_grav                        ! Thickness of the gravity current
  real(rp) :: c_gravr,c_gravu,radius_grav,c_grav_factor
  real(rp) :: time_eruption,vol_flow_rate,mass_flow_rate
  !
  !*** The minimum radius of the gravity current is important from the
  !*** numerical point of view. Too small radius results in strong asymmetry
  !*** of the deposit and/or large radial velocities (small dt)
  !
  min_radius = 1.5*sqrt(Hm1(1,1)*dX1*Hm1(1,1)*dX1+dX2*dX2)    ! Proportional to diagonal of the cell
  !
  SELECT CASE(coord_sys)
  case('LON-LAT')
     wind_center_x = radial_wind_lon
     wind_center_y = radial_wind_lat
  case('UTM')
     wind_center_x = radial_wind_x
     wind_center_y = radial_wind_y
  end SELECT
  !
  !*** Evaluate parameters for gravity current model
  !
  time_eruption  = time - beg_time                                           ! Time from start of eruption
  mass_flow_rate = sum(tmrat(1:nsou,1:nc))                                   ! MER (source already term scaled)
  !
  !*** Check if the source is on
  !
  if(mass_flow_rate.le.0.1) return
  !
  vol_flow_rate = c_flow_rate*sqrt(k_entrain)*mass_flow_rate**0.75_rp/ &
       (brunt_vaisala**0.625_rp)                                  ! Volumetric flow rate at NBL (q)
  !
  c_gravr     = 3.0_rp*lambda_grav*brunt_vaisala*vol_flow_rate/(2.0_rp*pi)   ! constant
  radius_grav = (c_gravr*time_eruption**2)**(0.33333_rp)                     ! Radius(t) of the front
  radius_grav = max(radius_grav,min_radius)                                  ! Avoid division by zero
  !
  c_gravu       = sqrt(2.0_rp*lambda_grav*brunt_vaisala*vol_flow_rate/(3.0_rp*pi))
  c_grav_factor = 0.75_rp*c_gravu*sqrt(radius_grav)
  th_grav       = c_gravu/(sqrt(radius_grav)*lambda_grav*brunt_vaisala)      ! Thickness
  !
  h_tot  = maxval(zs)                        ! Maximum column height (top of source)
  h_umbr = h_tot/1.32_rp                     ! Height of NBL (empirical). This is consistent with the plume model in SetSrc

  hmin   = max(h_umbr-0.5_rp*th_grav,0.5_rp*h_tot) ! Higher that H_tot/2 (empirical)
  hmax   = min(h_umbr+0.5_rp*th_grav,h_tot)     ! Lower than total column height
  !
  !*** Find nodes of the center of the umbrella
  !
  i_wind = (wind_center_x - xorigr)/dx
  j_wind = (wind_center_y - yorigr)/dy
  i_wind = max(1,i_wind)
  j_wind = max(1,j_wind)
  k_wind = 1
  store  = zlayer(nz)  ! A large value
  do k = 1,nz
     dist = abs(zlayer(k) - h_umbr)
     if(dist <= store) then
        k_wind = k
        store  = dist
     end if
  end do
  !
  !*** Wind at NBL above the crater
  !
  v_wind  = sqrt(vx(i_wind,j_wind,k_wind)**2 + vy(i_wind,j_wind,k_wind)**2)
  !
  !*** Radius_grav(Tp) where Tp=64/27*c_gravr/v^3. Transition at Ri=025 for passive transport transition
  !
  Tp = 64.0_rp/27.0_rp*c_gravr/max(v_wind**3,0.1_rp)
  Tb = Tp/8.0_rp
  max_radius = 1.7777_rp*c_gravr/max(v_wind**2,0.1_rp)
  !
  !*** log file
  !
  if(mpime == root ) then
     write(nlst,10) mass_flow_rate,vol_flow_rate,Tb/3600.,Tp/3600.,&
          max_radius/1e3_rp/4.0_rp,min(radius_grav/1e3_rp,max_radius/1e3_rp), &
          max_radius/1e3_rp,th_grav/1e3_rp,c_gravu/sqrt(radius_grav), &
          c_gravu/sqrt(min_radius)
10   format(/,&
          'GRAVITY CURRENT ',/, &
          '  Mass Flow Rate            (kg/s) : ',e13.6,/, &
          '  Volume flow rate          (m3/s) : ',e13.6,/, &
          '  Density driven time-scale (h   ) : ',e13.6,/, &
          '  Fully passive  time-scale (h   ) : ',e13.6,/, &
          '  Density driven GC radius  (km  ) : ',e13.6,/, &
          '  Current        GC radius  (km  ) : ',e13.6,/, &
          '  Fully passive  GC radius  (km  ) : ',e13.6,/, &
          '  Thickness                 (km  ) : ',e13.6,/, &
          '  Velocity front            (m/s ) : ',e13.6,/, &
          '  Maximum velocity          (m/s ) : ',e13.6)
  end if
  !
  !*** Find kmin and kmax (lower and upper z-index of the umbrella region)
  !
  kmin = 1
  store = zlayer(nz)  ! A large value
  do k = 1,nz
     dist = abs(zlayer(k) - hmin)
     if(dist <= store) then
        kmin = k
        store = dist
     end if
  end do
  kmax = nz
  store = zlayer(nz)  ! A large value
  do k = 1,nz
     dist = abs(zlayer(k) - hmax)
     if(dist <= store) then
        kmax = k
        store = dist
     end if
  end do
  !
  !*** Add radial wind
  !
  do j=1,ny
     do i=1,nx
        deltax = dX1*Hm1(i,j)*((xorigr + i*dx)-wind_center_x)/dx     ! distance (in m)
        deltay = dX2         *((yorigr + j*dy)-wind_center_y)/dy
        radius = sqrt(deltax*deltax + deltay*deltay)                 ! in meters
        !
        if(abs(radius) <= min_radius ) cycle ! Skip for R <= Rmin
        if(abs(radius) >  max_radius ) cycle ! Skip for R >  Rmin
        !
        rvel = c_grav_factor/radius*(1.0_rp+radius**2 /(3.0_rp*radius_grav**2))
        do k=kmin,kmax
           vx(i,j,k) = vx(i,j,k) + rvel * deltax/radius/Hm1(i,j)
           vy(i,j,k) = vy(i,j,k) + rvel * deltay/radius
        end do
     end do
  end do
  !
  return
end subroutine addradialwind
