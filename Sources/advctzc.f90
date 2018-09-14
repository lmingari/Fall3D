  subroutine advctzc(c,nx,ny,nz,nz_loc,iz_loc,vz,vset,dz,dt,lmeth)
  !***************************************************************************
  !*
  !*    Advection in the Z direction with non-uniform grid (Lax-Wendroff)
  !*
  !*    c        - Concentration (Input/Output matrix)
  !*    nx,ny,nz - Number of cells in the X,Y,Z directions (Input)
  !*    vz       - Z-component of the wind direction (Input Matrix)
  !*    dz       - Array of sizes of the cells in Z direction (Input)
  !*    dt       - Time step (Input)
  !*    topg     - Topography (Input matrix)
  !*    lmeth    - Limiter method (Input)   *** NOT USED IN THIS VERSION ***
  !*
  !***************************************************************************
  use KindType
  implicit none
  integer(ip) :: nx,ny,nz,lmeth,nz_loc,iz_loc
  real   (rp) :: dt
  real   (rp) :: c(0:nx+1,0:ny+1,0:nz+1),vz(nx,ny,nz),vset(nx,ny,nz),dz(0:nz)
  !
  integer(ip) :: i,j,k,kmax,kmin
  real   (rp) :: dzz,vk,vk_1,vk_12,Crk,deltar,deltal,sigma,dzk,dzk_1,dznz
  !
  real(rp), allocatable :: work(:,:,:,:)
  !
  !***  allocates memory
  !
  allocate( work( nx, ny, nz+1, 2 ) )
  !
  !***  Main loop
  !
  do k = MAX( 1, iz_loc-1 ), MIN( nz, iz_loc + nz_loc )
     !
     dzz   = 2.0_rp/(dz(k)+dz(k-1))
     dzk   = 1.0_rp/dz(k)
     dzk_1 = 1.0_rp/dz(k-1)
     kmax  = max(1 ,k-1)
     do j=1,ny
        do i=1,nx
           !
           vk    = vz(i,j,k   ) + vset(i,j,k   )
           vk_1  = vz(i,j,kmax) + vset(i,j,kmax)
           vk_12 = 0.5_rp*(vk+vk_1)
           !
           work(i,j,k,1) = dt*(vk*c(i,j,k) - vk_1*c(i,j,k-1))*dzz   ! first order
           !
           deltar = (c(i,j,k+1)-c(i,j,k  ))*dzk
           deltal = (c(i,j,k  )-c(i,j,k-1))*dzk_1
           sigma  = 0.5_rp*(sign(1.0_rp,deltar)+sign(1.0_rp,deltal))*min(abs(deltar),abs(deltal))

           Crk    = abs(vk_12)*dt*dzk_1
           work(i,j,k,2) = dt*0.5_rp*dz(k-1)*(1.0_rp-Crk)*sigma*abs(vk_12)*dzz
           !
        end do
     end do
  end do
  !
  dznz = 1.0_rp/dz(nz)
  do j=1,ny
     do i=1,nx
        vk  = vz(i,j,nz) + vset(i,j,nz)
        work(i,j,nz+1,1) = dt*vk*(c(i,j,nz+1) - c(i,j,nz))*dznz
     end do
  end do
  !
  do k = iz_loc, iz_loc + nz_loc - 1
     kmax = max(1 ,k-1)
     kmin = min(nz,k+1)
     do j=1,ny
        do i=1,nx
           vk = vz(i,j,k) + vset(i,j,k)
           if(vk > 0.0_rp) then
              c(i,j,k) = c(i,j,k) - work(i,j,k  ,1) - (work(i,j,k,2   )-work(i,j,kmax,2))
           else
              c(i,j,k) = c(i,j,k) - work(i,j,k+1,1) - (work(i,j,kmin,2)-work(i,j,k,2))
           endif
        end do
     end do
  end do
  !
  deallocate( work )
  !
  return
  end subroutine advctzc
