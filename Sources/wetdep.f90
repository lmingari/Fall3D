subroutine wetdep(c,cwet,ic,nx,ny,nz,nz_loc,iz_loc,dt,diam,zlayer,dz,zi,prate)
  !*************************************************************************
  !*
  !*    Computes wet deposition according to Jung and Shao (2006)
  !*      dC/dt = -L*C = a*(P**b)*C
  !*       a = 8.4d-5
  !*       b = 0.79
  !*       P precipitation rate in mmh-1
  !*    A critical cut-off size is assumed. As a first approach, wet deposition
  !*    is assumed below the PBL only. This is done because there is no
  !*    informtion about the height of the precipitation (only the total
  !*    rate is known).
  !*
  !*************************************************************************
  use KindType
  implicit none
  !
  integer(ip) :: ic,nx,ny,nz,nz_loc,iz_loc
  real   (rp) :: dt
  real   (rp) :: c(0:nx+1,0:ny+1,0:nz+1),cwet(nx,ny,0:*)
  real   (rp) :: diam(*),zlayer(nz),dz(0:nz),zi(nx,ny),prate(nx,ny)

  integer(ip) :: ix,iy,iz                           ! Work integers
  real   (rp) :: a,b,diam_max,lambda
  !
  !*** Check if wet deposition applies for this particle size
  !
  diam_max = 100.0_rp                                    !  cutt-off size in um
  diam_max = diam_max*1d-6                               !  cutt-off size in m
  if( diam(ic).gt.diam_max ) return
  !
  !*** Computes wet deposition
  !
  a = 8.4d-5
  b = 0.79_rp
  !
  do ix = 1,nx
     do iy = 1,ny
        !
        lambda = min(dt*a*(prate(ix,iy)**b),1.0_rp)
        !
        do iz = iz_loc,iz_loc + nz_loc -1
           !
           if(zlayer(iz).le.zi(ix,iy)) then
              cwet(ix,iy,ic) = cwet(ix,iy,ic) + dz(iz)*lambda*c(ix,iy,iz)
              c   (ix,iy,iz) = c   (ix,iy,iz)*(1.0_rp-lambda)
           end if
        end do
     end do
  end do
  !
  return
end subroutine wetdep
