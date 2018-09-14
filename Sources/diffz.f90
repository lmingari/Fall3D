  subroutine diffz(c,nx,ny,nz,nz_loc,iz_loc,dz,dt,rkver,rho,work)
  !************************************************************************
  !*
  !*    Computes diffusion (Z)
  !*
  !************************************************************************
  use KindType
  implicit none
  !
  integer(ip) :: nx,ny,nz,nz_loc,iz_loc
  real   (rp) :: c(0:nx+1,0:ny+1,0:nz+1),dz(0:nz),rkver(nx,ny,nz),rho(nx,ny,nz),work(nz)
  !
  integer(ip) :: i,j,k
  real   (rp) :: dt,rksb,rkst,gradb,gradt
  real   (rp), allocatable :: l_dz(:),l_dzdz(:)              ! work arrays
  !
  !***  Allocates and computes temporal arrays
  !
  allocate ( l_dz  (0:nz))
  allocate ( l_dzdz(nz  ))
  !
  l_dz(0:nz) = 1.0_rp/dz(0:nz)
  do k=1,nz
     l_dzdz(k) = 2.0_rp/(dz(k)+dz(k-1))
  end do
  !
  do j=1,ny
     do i=1,nx
        !
        !           Bottom boundary
        !
        if(iz_loc == 1) then
           rkst    = rkver(i,j,2)*rho(i,j,2)
           rksb    = rkver(i,j,1)*rho(i,j,1)
           gradt   = rkst*( (c(i,j,2)/rho(i,j,2)) - (c(i,j,1)/rho(i,j,1)) )*l_dz(1)
           gradb   = rksb*( (c(i,j,1)/rho(i,j,1)) - (c(i,j,0)/rho(i,j,1)) )*l_dz(0)
           gradb   = max(0.0_rp,gradb)
           work(1) = dt*(gradt-gradb)*l_dzdz(1)
        end if
        !
        do k= MAX(2,iz_loc), MIN(iz_loc + nz_loc - 1, nz-1)
           rkst    = 0.5_rp*(rho(i,j,k+1)*rkver(i,j,k+1)+ rho(i,j,k  )*rkver(i,j,k  ))
           rksb    = 0.5_rp*(rho(i,j,k  )*rkver(i,j,k  )+ rho(i,j,k-1)*rkver(i,j,k-1))
           gradt   = rkst*( (c(i,j,k+1)/rho(i,j,k+1)) - (c(i,j,k  )/rho(i,j,k  )) ) *l_dz(k  )
           gradb   = rksb*( (c(i,j,k  )/rho(i,j,k  )) - (c(i,j,k-1)/rho(i,j,k-1)) ) *l_dz(k-1)
           work(k) = dt*(gradt-gradb)*l_dzdz(k)
        end do
        !
        !           Top boundary
        !
        if(iz_loc + nz_loc - 1  ==  nz) then
           rkst  = rkver(i,j,nz  )*rho(i,j,nz  )
           rksb  = rkver(i,j,nz-1)*rho(i,j,nz-1)
           gradt = rkst*( (c(i,j,nz+1)/rho(i,j,nz)) - (c(i,j,nz  )/rho(i,j,nz  )) ) *l_dz(nz  )
           gradb = rksb*( (c(i,j,nz  )/rho(i,j,nz)) - (c(i,j,nz-1)/rho(i,j,nz-1)) ) *l_dz(nz-1)
           gradt = min(0.0_rp,gradt)
           work(nz) = dt*(gradt-gradb)*l_dz(nz)
        end if
        !
        do k=iz_loc, iz_loc + nz_loc - 1
           c(i,j,k) = c(i,j,k) + work(k)
        enddo
     enddo
  enddo
  !
  deallocate( l_dz   )
  deallocate( l_dzdz )
  !
  return
  end subroutine diffz
