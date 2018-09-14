subroutine corvver
  !************************************************************
  !*
  !*     Correct vertical velocity. This is a transformation
  !*     from the cartesian (x,y,z) to the terrain following
  !*     coordinate system (X;Y;Z). Note that J=1 is assumed because
  !*     we do not use the classical sigma maping but Z = z - zo
  !*     (i.e. the factor H/(H-zo) is not present)
  !*
  !*************************************************************
  use KindType
  use Numeric
  use Master
  implicit none
  !
  integer(ip) :: ix,iy,iz,ic
  real   (rp) :: dxinv,dyinv,vs,zx,zy
  !
  dyinv = 1.0_rp/dX2
  !
  do iy = 1,ny
     do ix = 1,nx
        dxinv = 1.0_rp/(Hm1(ix,iy)*dX1)
        zx    = dxinv*(topg(min(ix+1,nx),iy)-topg(ix,iy))
        zy    = dyinv*(topg(ix,min(iy+1,ny))-topg(ix,iy))
        do iz = 1,nz                                       ! not in the first layer
           vs = zx*vx(ix,iy,iz)+zy*vy(ix,iy,iz)
  !        do ic = 1,nc                                    ! correction affects to all classes
  !           vset(ix,iy,iz,ic) = vset(ix,iy,iz,ic) - vs
  !        end do
           vz(ix,iy,iz) = vz(ix,iy,iz) - vs
        end do
     end do
  end do
  !
  return
end subroutine corvver
