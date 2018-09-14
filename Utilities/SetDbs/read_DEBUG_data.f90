subroutine read_DEBUG_data(fname,time_lag,time,timesec,nx,ny,nz,zlayer,lon,lat,topg,LDU, &
     u,v,w,T,Tv,Tp,qv,ro,p,pblh,ust,L)
  !********************************************************************
  !*
  !*   Set data for a DEBUG case (ad hoc)
  !*
  !********************************************************************
  use KindType, ONLY :  ip,rp
  use InpOut
  use MathFun
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
  logical     :: standard
  integer(ip) :: ix,iy,iz
  real(rp)    :: z,TT,pp,roo
  !
  !***  Topography
  !
  topg = 0.
  !
  !***  velocity
  !
  u = 10.
  v = 0.
  w = 0.
  !
  !***  Atmosphere (T,Tp,p,ro)
  !
  standard = .false.

  if(.not.standard) then
     T = 273.
     p = 1e5
     TP = T
     ro = 1.
  else
     do iz = 1,nz
        !
        z = zlayer(iz)   ! terrain following
        do iy = 1,ny
           do ix = 1,nx
              z = zlayer(iz) + topg(ix,iy)   ! asl
              !
              !***  Standard pressure and density (standard atmosphere assumed)
              !
              call standard_atm(PP,TT,roo,z)
              !
              T (ix,iy,iz) = TT
              p (ix,iy,iz) = PP
              ro(ix,iy,iz) = roo
              !
              !***  Potential temperature
              !
              Tp(ix,iy,iz) =  T(ix,iy,iz)*( (1e5/p(ix,iy,iz))**(0.285) )
              !
           end do
        end do
     end do
  end if
  !
  !***  Other variables not given
  !
  LDU  = -999
  qv   = 0.     ! no water
  Tv   = T
  pblh = 1000.
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
end subroutine read_DEBUG_data
