  subroutine getnsrc(nsrc)
  !****************************************************
  !*
  !*    Returns the number of sources
  !*
  !****************************************************
  use KindType
  use Master
  implicit none
  !
  logical     :: found
  integer(ip) :: nsrc,ix,iy,iz,ic
  !
  nsrc=0
  do iz = 1,nz
     do iy = 1,ny
        do ix = 1,nx
           found = .false.
           do ic = 1,nc
              if(src(ic,ix,iy,iz).gt.0.0) found = .true.
           end do
           if(found) nsrc = nsrc + 1
        end do
     end do
  end do
  !
  return
  end subroutine getnsrc
