subroutine setbcc(c,nx,ny,nz,nz_loc,iz_loc,vx,vy,vz,vset,vdep,rkver,dz)
  !*********************************************************************
  !*
  !*    Set the Boundary conditions
  !*
  !********************************************************************
  use KindType
  implicit none
  !
  integer(ip) :: nx,ny,nz,nz_loc,iz_loc
  integer(ip) :: i,j,k
  real   (rp) :: vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz)
  real   (rp) :: c(0:nx+1,0:ny+1,0:nz+1),vset(nx,ny,nz),vdep(nx,ny)
  real   (rp) :: rkver(nx,ny,nz)
  real   (rp) :: dz(0:nz)
  !
  !***  Bottom
  !
  if( iz_loc == 1 ) then
     do j=1,ny
        do i=1,nx
           if((vz(i,j,1)+vset(i,j,1)) >= 0.0_rp) then
              c(i,j,0) = 0.0_rp
           else
              !             c(i,j,0)    = c(i,j,1)
              !c(i,j,0)  = c(i,j,1)/ (1.0_rp + vdep(i,j)/vset(i,j,1))   ! dry deposition mechanism
              c(i,j,0) = c(i,j,1)*(1.0_rp+dz(0)*vdep(i,j)/rkver(i,j,1))
           endif
        enddo
     enddo
  else
     c(:,:,0) = 0.0_rp
  end if
  !
  !***  Top
  !
  if( ( iz_loc + nz_loc - 1 ) == nz ) then
     do j=1,ny
        do i=1,nx
           if((vz(i,j,nz)+vset(i,j,nz)) <= 0.0_rp) then
              c(i,j,nz+1) = 0.0_rp
           else
              c(i,j,nz+1) = c(i,j,nz)
           endif
        enddo
     enddo
  else
     c(:,:,nz+1) = 0.0_rp
  end if
  !
  !***  Left/Right (X)
  !
  do k= iz_loc, iz_loc + nz_loc - 1
     do j=1,ny
        if(vx(1,j,k) >= 0.0_rp) then
           c(0,j,k)    = 0.0_rp
        else
           c(0,j,k)    = c(1,j,k)
        endif
        if(vx(nx,j,k) <= 0.0_rp) then
           c(nx+1,j,k) = 0.0_rp
        else
           c(nx+1,j,k) = c(nx,j,k)
        endif
     enddo
  enddo
  !
  !***  Up/Down (Y)
  !
  do k= iz_loc, iz_loc + nz_loc - 1
     do i=1,nx
        if(vy(i,1,k) >= 0.0_rp) then
           c(i,0,k)    = 0.0_rp
        else
           c(i,0,k)    = c(i,1,k)
        endif
        if(vy(i,ny,k) <= 0.0_rp) then
           c(i,ny+1,k) = 0.0_rp
        else
           c(i,ny+1,k) = c(i,ny,k)
        endif
     enddo
  enddo
  !
  !***  Corners
  !
  if( iz_loc == 1 ) then
     c(0,0,0)          = c(1,1,1)
     c(nx+1,0,0)       = c(nx,1,1)
     c(0,ny+1,0)       = c(1,ny,1)
     c(nx+1,ny+1,0)    = c(nx,ny,1)
  else
     c(0,0,0)          = 0.0_rp
     c(nx+1,0,0)       = 0.0_rp
     c(0,ny+1,0)       = 0.0_rp
     c(nx+1,ny+1,0)    = 0.0_rp
  end if
  !
  if( ( iz_loc + nz_loc - 1 ) == nz ) then
     c(0,0,nz+1)       = c(1,1,nz)
     c(nx+1,0,nz+1)    = c(nx,1,nz)
     c(0,ny+1,nz+1)    = c(1,ny,nz)
     c(nx+1,ny+1,nz+1) = c(nx,ny,nz)
  else
     c(0,0,nz+1)       = 0.0_rp
     c(nx+1,0,nz+1)    = 0.0_rp
     c(0,ny+1,nz+1)    = 0.0_rp
     c(nx+1,ny+1,nz+1) = 0.0_rp
  end if
  !
  return
end subroutine setbcc
