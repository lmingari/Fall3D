subroutine ter2asl(conc,topg,zlayer,nx,ny,nz)
  !***************************************************************************
  !*
  !*   Converts a 3d array from terrain following coordinates to absolute
  !*   coordinates
  !*
  !***************************************************************************
  use KindType
  implicit none
  integer(ip) :: nx,ny,nz
  real   (rp) :: conc(nx,ny,nz),topg(nx,ny),zlayer(nz)
  !
  logical     :: found
  integer(ip) :: ix,iy,iz,izz
  real   (rp) :: zo,zasl,s
  real   (rp), allocatable :: c(:,:,:)
  !
  !*** Allocates memory for temporal array and loads terrain following data
  !
  allocate(c(nx,ny,nz))
  c(1:nx,1:ny,1:nz) = conc(1:nx,1:ny,1:nz)
  !
  !*** Loop
  !
  do iy = 1,ny
     do ix = 1,nx
        zo = max(0.0_rp,topg(ix,iy))
        do iz = 1,nz
           zasl = zlayer(iz)    ! this is the point to find
           !
           if(zasl < (zo+zlayer(1))) then
              conc(ix,iy,iz) = 0.0_rp
           else
              found = .false.
              izz   = 0
              do while(.not.found)
                 izz = izz + 1
                 if( (zasl >= (zo+zlayer(izz))) .and. (zasl <= (zo+zlayer(izz+1))) ) then
                    found = .true.
                    s = (zo+zlayer(izz+1)-zasl)/(zlayer(izz+1)-zlayer(izz))
                    conc(ix,iy,iz) = s*c(ix,iy,izz) + (1.0_rp-s)*c(ix,iy,izz+1)
                 end if
                 if( (izz+1) == nz .and. (.not.found)) then
                    conc(ix,iy,iz) = 0.0_rp
                    found = .true.
                 end if
              end do
           end if
        end do
     end do
  end do
  !
  deallocate(c)
  return
end subroutine ter2asl
