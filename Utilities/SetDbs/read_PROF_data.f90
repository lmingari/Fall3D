subroutine read_PROF_data(fname,time_lag,time,timesec,nx,ny,nz,zlayer,lon,lat,topg,LDU, &
     u,v,w,T,Tv,Tp,qv,ro,p,pblh,ust,L)
  !********************************************************************
  !*
  !*   Reads PROF data from a profile. This routine MUST be called
  !*   after read_PROF_grid
  !*
  !********************************************************************
  use KindType, ONLY :  ip,rp
  use InpOut
  use MathFun
  use PROF_nc
  implicit none
  !
  character(len=*) :: fname
  integer(ip)      :: nx,ny,nz
  real(rp)         :: time_lag,time,timesec
  real(rp)         :: zlayer(nz),lat(nx,ny),lon(nx,ny),topg(nx,ny),LDU(nx,ny), &
       u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz),T(nx,ny,nz),Tv(nx,ny,nz), &
       Tp(nx,ny,nz),qv(nx,ny,nz),ro(nx,ny,nz),p(nx,ny,nz), &
       pblh(nx,ny),ust(nx,ny),L(nx,ny)
  !
  logical     :: go_on
  integer(ip) :: ix,iy,iz,np,jp
  real(rp)    :: itime1,itime2
  real(rp)    :: z,ux,uy,TT,PP,roo
  integer(ip), parameter :: npmax=100
  real(rp)    :: windz(npmax),windx(npmax),windy(npmax)
  !
  !***  Opens the file
  !
  open(90,FILE=TRIM(fname),STATUS='old',ERR=100)
  !
  !***  Reads until the required time
  !
  read(90,*,ERR=101) rvoid  ! xo,yo
  read(90,*,ERR=101) cvoid  ! YYYYMMDD
  !
  go_on = .true.
  do while(go_on)
     !
     read(90,*,ERR=101) itime1,itime2
     read(90,*,ERR=101) np
     if(np.gt.npmax) call runend('get_wind_profile_point: To many layers for wind profile. npmax=100')
     !
     !***    Read data for this time interval
     !
     do jp = 1,np
        read(90,*,ERR=101) windz(jp),windx(jp),windy(jp)
     end do
     !
     if((time_lag+timesec.ge.itime1).and.(time_lag+timesec.le.itime2)) go_on = .false.
     !
  end do
  close(90)
  !
  !***  Reads profile data (u,v) for the current dbs interval
  !***  Rest variables are interpolated in terrain following
  !
  do iz = 1,nz
     !
     z = zlayer(iz)   ! terrain following
     call get_wind_profile_point(np,windz,windx,windy,z,ux,uy)
     u(1:nx,1:ny,iz) = ux
     v(1:nx,1:ny,iz) = uy
     w(1:nx,1:ny,iz) = 0.0
     !
     do iy = 1,ny
        do ix = 1,nx
           z = zlayer(iz) + topg(ix,iy)   ! asl
           !
           !***            Standard pressure and density (standard atmosphere assumed)
           !

           call standard_atm(PP,TT,roo,z)
           !
           T (ix,iy,iz) = TT
           p (ix,iy,iz) = PP
           ro(ix,iy,iz) = roo
           !
           !***            Potential temperature
           !
           Tp(ix,iy,iz) =  T(ix,iy,iz)*( (1e5/p(ix,iy,iz))**(0.285) )
           !
        end do
     end do
  end do
  !
  ! Rotate wind if request
  if(rotate_wind) then
     allocate(unew(nx,ny))
     allocate(vnew(nx,ny))
     do iz=1,nz
        do iy=1,ny
           do ix = 1,nx
              call rotate2d(u(ix,iy,iz),v(ix,iy,iz),unew(ix,iy),vnew(ix,iy),rotation_angle)
           end do
        end do
        u(:,:,iz) = unew      ! Copy
        v(:,:,iz) = vnew      ! Copy
     end do
     deallocate(unew)
     deallocate(vnew)
  end if
  !
  !***  Other variables not given
  !
  LDU  = -999
  qv   = 0.     ! no water
  Tv   = T
  pblh = 1000.
  !
  !*** ust(nx,ny) and L(nx,ny). Estimated
  !
  do iy = 1,ny
     do ix = 1,nx
        call get_par_ABL(1.0_rp ,zlayer(2),T(ix,iy,1),T(ix,iy,2),P(ix,iy,1),P(ix,iy,2), &
             u(ix,iy,1),v(ix,iy,1),ust(ix,iy),L(ix,iy))
     end do
  end do
  !
  !***  Successful end
  !
  return
  !
  !***  List of errors
  !
100 call runend('get_wind_profile_point: Error opening the input file '//TRIM(fname))
101 call runend('get_wind_profile_point: Error reading winds file at time ')

  return
  !
end subroutine read_PROF_data
!
!
!
subroutine get_wind_profile_point(np,windz,windx,windy,z,ux,uy)
  !**************************************************************************
  !*
  !*    Gets the wind components and temperature at a certain height from
  !*    a wind profile file
  !*
  !*****************************************************************************
  use KindType, ONLY :  ip,rp
  implicit none
  !
  integer(ip) :: np,kp
  real(rp) :: z,ux,uy
  real(rp) :: windz(*),windx(*),windy(*)
  !
  logical  :: go_on
  real(rp) :: s
  !
  !***  Interpolate values at z
  !
  if(z.le.windz(1)) then
     kp = 1
     s  = 1.0
  else if(z.ge.windz(np)) then
     kp = np-1
     s  = 0.0
  else
     go_on = .true.
     kp = 0
     do while(go_on)
        kp = kp + 1
        if((z.ge.windz(kp)).and.(z.le.windz(kp+1))) then
           go_on = .false.
           s = (windz(kp+1)-z)/(windz(kp+1)-windz(kp))
        end if
     end do
  end if
  !
  !***  Interpolates
  !
  ux = s*windx(kp) + (1.0-s)*windx(kp+1)
  uy = s*windy(kp) + (1.0-s)*windy(kp+1)
  !
  return
end subroutine get_wind_profile_point
!
!
!
subroutine standard_atm(p,T,rho,z)
  !******************************************************************************
  !*
  !*    Returns p, T and rho at height z (asl) for the ICAO standard atmosphrere
  !*    Only three layers are used (valid up to 32km)
  !*
  !******************************************************************************
  use KindType
  implicit none
  !
  real(rp) :: p,T,rho,z
  !
  real(rp) :: goo = 9.80665_rp    ! standard gravity
  real(rp) :: Roo = 8.31432_rp    ! universal gas constant Nm / (molK)
  real(rp) :: Moo = 0.0289644_rp  ! molar mass of air   (kg/mol)
  !
  real(rp) :: zoo =  0.0_rp       ! First layer  (Troposphere)
  real(rp) :: Too = 288.15_rp
  real(rp) :: poo = 101325.0_rp
  real(rp) :: Loo = -0.0065_rp
  !
  real(rp) :: zo1 = 11000.0_rp     ! Second layer (Stratosphere)
  real(rp) :: To1 = 216.65_rp
  real(rp) :: po1 = 22632.0_rp
  !
  real(rp) :: zo2 = 20000.0_rp     ! Third layer (Stratosphere)
  real(rp) :: To2 = 216.65_rp
  real(rp) :: po2 = 5474.0_rp
  real(rp) :: Lo2 = 0.001_rp
  !
  !***  computes
  !
  if(z.lt.zo1) then
     T = Too + Loo*(z-zoo)
     P = poo*( (Too/T)**(goo*Moo/(Roo*Loo)) )
     rho = P*Moo/(Roo*T)
  else if(z.lt.zo2) then
     T = To1
     P = Po1*exp(-goo*Moo*(z-zo1)/(Roo*To1))
     rho = P*Moo/(Roo*T)
  else
     T = To2 + Lo2*(z-zo2)
     P = po2*( (To2/T)**(goo*Moo/(Roo*Lo2)) )
     rho = P*Moo/(Roo*T)
  end if
  !
  return
end  subroutine standard_atm
