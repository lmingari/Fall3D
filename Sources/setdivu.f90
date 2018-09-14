  subroutine setdivu
  !************************************************************************
  !*
  !*  Computes the velocity field divergence (wind field plus possible
  !*  gravity current)
  !*
  !************************************************************************
  use KindType
  use Numeric
  use Master
  implicit none
  !
  integer(ip) :: i,j,k
  real   (rp) :: vi,vi_1,vj,vj_1,vk,vk_1
  !
  do i = 1,nx
     do j = 1,ny
        do k = 1,nz
           vi    = vx(i,j,k)
           vi_1  = vx(max(1,i-1),j,k)
           vj    = vy(i,j,k)
           vj_1  = vy(i,max(1,j-1),k)
           vk    = vz(i,j,k)           ! +vset(i,j,k,ic)
           vk_1  = vz(i,j,max(1,k-1))  ! +vset(i,j,max(1,k-1),ic)
           !
           divu(i,j,k) = (vi-vi_1)/(Hm1(i,j)*dX1) + (vj-vj_1)/dX2 + (vk-vk_1)/dz(k)
           !
        end do
     end do
  end do
  !
  return
  end subroutine setdivu
