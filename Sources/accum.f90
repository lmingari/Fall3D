subroutine accum(c,cdry,ic,nx,ny,nz,nz_loc,iz_loc,vz,dt,vset,vdep)
  !*************************************************************************
  !*
  !*    Computes ground accumulation by sedimentation and dry deposition
  !*
  !*************************************************************************
  use KindType
  implicit none
  !
  integer(ip) :: ic,nx,ny,nz,nz_loc,iz_loc
  integer(ip) :: i,j                        ! Work integers
  real   (rp) :: vz(nx,ny,nz)
  !
  real(rp) :: dt
  real(rp) :: c(0:nx+1,0:ny+1,0:nz+1),cdry(nx,ny,0:*),vset(nx,ny,nz),vdep(nx,ny)
  !
  if( iz_loc == 1 ) then
     do i=1,nx
        do j=1,ny
           !cdry(i,j,ic) = cdry(i,j,ic) - dt*0.5_rp*(c(i,j,1)*min(0.0_rp,vset(i,j,1))+ &
           !     c(i,j,0)*min(0.0_rp,vset(i,j,1)+vdep(i,j)))
           cdry(i,j,ic) = cdry(i,j,ic) - dt*c(i,j,1)*min(0.0_rp,vset(i,j,1)+vdep(i,j))
        enddo
     enddo
  end if
  !
  return
end subroutine accum
