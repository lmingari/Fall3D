  subroutine diffx(c,nx,ny,nz,nz_loc,iz_loc,dx,dt,rkhor,rho,work)
  !**********************************************************************
  !*
  !*    Computes diffusion (X)
  !*
  !**********************************************************************
  use KindType
  implicit none
  integer(ip) :: nx,ny,nz,nz_loc,iz_loc
  integer(ip) :: i,j,k                        ! work variables
  real   (rp) :: dx,dt
  real   (rp) :: dtdx2,rksr,rksl,gradr,gradl      ! work variables
  real   (rp) :: c(0:nx+1,0:ny+1,0:nz+1),rkhor(nx,ny,nz),rho(nx,ny,nz),work(nx)
  !
  real   (rp), allocatable :: l_rho(:)
  !
  allocate ( l_rho(nx) )
  dtdx2=dt/dx**2
  !
  do k=iz_loc, iz_loc + nz_loc - 1
     do j=1,ny
        !
        l_rho(1:nx) = 1.0_rp/rho(1:nx,j,k)       ! work array
        !
        !           Left boundary
        !
        rksr    = rkhor(2,j,k)*rho(2,j,k)
        rksl    = rkhor(1,j,k)*rho(1,j,k)
        gradr   = rksr*( c(2,j,k)*l_rho(2) - c(1,j,k)*l_rho(1) )
        gradl   = rksl*( c(1,j,k)*l_rho(1) - c(0,j,k)*l_rho(1) )   ! rho(0,i,j) undefined
        gradl   = max(0.0_rp,gradl)
        work(1) = dtdx2*(gradr-gradl)
        !
        do i=2,nx-1
           rksr    = 0.5_rp*(rho(i+1,j,k)*rkhor(i+1,j,k)+ rho(i  ,j,k)*rkhor(i  ,j,k))
           rksl    = 0.5_rp*(rho(i  ,j,k)*rkhor(i  ,j,k)+ rho(i-1,j,k)*rkhor(i-1,j,k))
           gradr   = rksr*( c(i+1,j,k)*l_rho(i+1) - c(i  ,j,k)*l_rho(i  ) )
           gradl   = rksl*( c(i  ,j,k)*l_rho(i  ) - c(i-1,j,k)*l_rho(i-1) )
           work(i) = dtdx2*(gradr-gradl)
        end do
        !
        !           Right boundary
        !
        rksr  = rkhor(nx  ,j,k)*rho(nx  ,j,k)
        rksl  = rkhor(nx-1,j,k)*rho(nx-1,j,k)
        gradr = rksr*( c(nx+1,j,k)*l_rho(nx) - c(nx  ,j,k)*l_rho(nx  ) )
        gradl = rksl*( c(nx  ,j,k)*l_rho(nx) - c(nx-1,j,k)*l_rho(nx-1) )
        gradr = min(0.0_rp,gradr)
        work(nx) = dtdx2*(gradr-gradl)
        !
        do i=1,nx
           c(i,j,k) = c(i,j,k) + work(i)
        enddo
     enddo
  enddo
  !
  deallocate ( l_rho )
  !
  return
  end subroutine diffx
