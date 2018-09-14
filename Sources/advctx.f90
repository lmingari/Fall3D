  subroutine advctx(c,nx,ny,nz,nz_loc,iz_loc,vx,dx,dt,work,lmeth)
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
  real   (rp) :: dx,dt
  real   (rp) :: c(0:nx+1,0:ny+1,0:nz+1),vx(nx,ny,nz),work(nx+1,2)
  !
  integer(ip) :: i,j,k
  real   (rp) :: auxv,vi,vi_1,vi_12,Cri,deltar,deltal,sigma
  !
  auxv = dt/dx              ! Auxialiary variable
  !
  !***  Main loop
  !
  do k = iz_loc, iz_loc + nz_loc - 1
     do j=1,ny
        do i=1,nx
           !
           vi    = vx(i          ,j,k)
           vi_1  = vx(max(1 ,i-1),j,k)
           vi_12 = 0.5_rp*(vi+vi_1)
           !
           work(i,1) = auxv*(vi*c(i,j,k) - vi_1*c(i-1,j,k))   ! first order
           !
           deltar = (c(i+1,j,k)-c(i  ,j,k))/dx
           deltal = (c(i  ,j,k)-c(i-1,j,k))/dx
           sigma  = 0.5_rp*(sign(1.0_rp,deltar)+sign(1.0_rp,deltal))*min(abs(deltar),abs(deltal))
           !
           Cri       = auxv*abs(vi_12)
           work(i,2) = dt*0.5_rp*(1.0_rp-Cri)*sigma*abs(vi_12)

        end do
        work(nx+1,1) = auxv*vi*(c(nx+1,j,k)-c(nx,j,k))
        !
        do i=1,nx
           if(vx(i,j,k) > 0.0_rp) then
              c(i,j,k) = c(i,j,k) - work(i  ,1) - (work(i,2)-work(max(1,i-1),2))
           else
              c(i,j,k) = c(i,j,k) - work(i+1,1) - (work(min(nx,i+1),2)-work(i,2))
           endif
        enddo
     enddo
  enddo
  !
  return
  end subroutine advctx
