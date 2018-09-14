  subroutine diffy(c,nx,ny,nz,nz_loc,iz_loc,dy,dt,rkhor,rho,work)
  !***********************************************************************
  !*
  !*    Computes diffusion (Y)
  !*
  !**********************************************************************
  use KindType
  implicit none
  integer(ip) :: nx,ny,nz,nz_loc,iz_loc
  integer(ip) :: i,j,k                ! work variables
  real   (rp) :: dy,dt
  real   (rp) :: dtdy2,rksu,rksd,gradu,gradd     ! work variables
  real   (rp) :: c(0:nx+1,0:ny+1,0:nz+1),rkhor(nx,ny,nz),rho(nx,ny,nz),work(ny)
  !
  real   (rp), allocatable :: l_rho(:)
  !
  allocate ( l_rho(ny) )
  dtdy2=dt/dy**2
  !
  do k=iz_loc, iz_loc + nz_loc - 1
     do i=1,nx
        !
        l_rho(1:ny) = 1.0_rp/rho(i,1:ny,k)       ! work array
        !
        !           Down boundary
        !
        rksu    = rkhor(i,2,k)*rho(i,2,k)
        rksd    = rkhor(i,1,k)*rho(i,1,k)
        gradu   = rksu*( c(i,2,k)*l_rho(2) - c(i,1,k)*l_rho(1) )
        gradd   = rksd*( c(i,1,k)*l_rho(1) - c(i,0,k)*l_rho(1) ) ! rho(0,i,j) undefined
        gradd   = max(0.0_rp,gradd)
        work(1) = dtdy2*(gradu-gradd)
        !
        do j=2,ny-1
           rksu    = 0.5_rp*(rho(i,j+1,k)*rkhor(i,j+1,k)+ rho(i,j  ,k)*rkhor(i,  j,k))
           rksd    = 0.5_rp*(rho(i,j  ,k)*rkhor(i,j  ,k)+ rho(i,j-1,k)*rkhor(i,j-1,k))
           gradu   = rksu*( c(i,j+1,k)*l_rho(j+1) - c(i,j  ,k)*l_rho(j  ) )
           gradd   = rksd*( c(i,j  ,k)*l_rho(j  ) - c(i,j-1,k)*l_rho(j-1) )
           work(j) = dtdy2*(gradu-gradd)
        end do
        !
        !           Up boundary
        !
        rksu  = rkhor(i,ny  ,k)*rho(i,ny  ,k)
        rksd  = rkhor(i,ny-1,k)*rho(i,ny-1,k)
        gradu = rksu*( c(i,ny+1,k)*l_rho(ny) - c(i,ny  ,k)*l_rho(ny  ) )
        gradd = rksd*( c(i,ny  ,k)*l_rho(ny) - c(i,ny-1,k)*l_rho(ny-1) )
        gradu = min(0.0_rp,gradu)
        work(ny) = dtdy2*(gradu-gradd)
        !
        do j=1,ny
           c(i,j,k) = c(i,j,k) + work(j)
        enddo
     enddo
  enddo
  !
  deallocate ( l_rho )
  !
  return
  end subroutine diffy
