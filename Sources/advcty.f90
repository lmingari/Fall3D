  subroutine advcty(c,nx,ny,nz,nz_loc,iz_loc,vy,dy,dt,work,lmeth)
  !***************************************************************************
  !*
  !*    Advection in the X direction with uniform grid (Lax-Wendroff)
  !*
  !*    c        - Concentration (Input/Output matrix)
  !*    nx,ny,nz - Number of cells in the X,Y,Z directions (Input)
  !*    vx       - X-component of the wind direction (Input Matrix)
  !*    dx       - Dimension of the cells in X direction (Input)
  !*    dt       - Time step (Input)
  !*    work     - Auxiliary vector (Input)
  !*    lmeth    - Limiter method (Input)   *** NOT USED IN THIS VERSION ***
  !*
  !***************************************************************************
  use KindType
  implicit none
  integer(ip) :: nx,ny,nz,lmeth,nz_loc,iz_loc
  real   (rp) :: dy,dt
  real   (rp) :: c(0:nx+1,0:ny+1,0:nz+1),vy(nx,ny,nz),work(ny+1,2)
  !
  integer(ip) :: i,j,k
  real   (rp) :: auxv,vj,vj_1,vj_12,Crj,deltar,deltal,sigma
  !
  auxv = dt/dy              ! Auxialiary variable
  !
  !***  Main loop
  !
  do k = iz_loc, iz_loc+nz_loc-1
     do i=1,nx
        do j=1,ny
           !
           vj    = vy(i,j          ,k)
           vj_1  = vy(i,max(1 ,j-1),k)
           vj_12 = 0.5_rp*(vj+vj_1)
           !
           work(j,1) = auxv*(vj*c(i,j,k) - vj_1*c(i,j-1,k))   ! first order
           !
           deltar = (c(i,j+1,k)-c(i,j  ,k))/dy
           deltal = (c(i,j  ,k)-c(i,j-1,k))/dy
           sigma  = 0.5_rp*(sign(1.0_rp,deltar)+sign(1.0_rp,deltal))*min(abs(deltar),abs(deltal))
           !
           Crj       = auxv*abs(vj_12)
           work(j,2) = dt*0.5_rp*(1.0_rp-Crj)*sigma*abs(vj_12)
        end do
        work(ny+1,1) = auxv*vj*(c(i,ny+1,k)-c(i,ny,k))
        !
        do j=1,ny
           if(vy(i,j,k) > 0.0_rp) then
              c(i,j,k) = c(i,j,k) - work(j  ,1) - (work(j,2)-work(max(1,j-1),2))
           else
              c(i,j,k) = c(i,j,k) - work(j+1,1) - (work(min(ny,j+1),2)-work(j,2))
           endif
        end do
     end do
  end do
  !
  return
  end subroutine advcty
