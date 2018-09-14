  subroutine wrisrc_mesh(isec1,isec2)
  !***************************************************************
  !*
  !*    Writes the source file during a time step (isec1,isec2)
  !*
  !***************************************************************
  use KindType
  use InpOut
  use Master
  implicit none
  !
  integer(ip) :: isec1,isec2
  logical     :: found
  integer(ip) :: nsrc,ix,iy,iz,ic
  real(rp)    :: x,y,z,MFR
  !
  !***  Computes the source vector src(nc,nx,ny,nz)
  !
  call getsrc
  !
  !***  Gets the number of sources
  !
  call getnsrc(nsrc)
  !
  !***  Total mass flow rate
  !
  MFR = SUM(src)
  !
  !***  Writes header for the current interval
  !
  write(lusrc,10) isec1,isec2
  write(lusrc,11) nsrc,nc
  write(lusrc,12) MFR
10 format(i7,1x,i7)
11 format(i7,1x,i7)
12 format(e16.9)
  !
  !***  Writes the rest of file
  !
  do iz = 1,nz
     do iy = 1,ny
        do ix = 1,nx
           found = .false.
           do ic = 1,nc
              if(src(ic,ix,iy,iz).gt.0.0) found = .true.
           end do
           if(found) then
              x = xorigr + (ix-1)*dx
              y = yorigr + (iy-1)*dy
              z = zmodel(iz)
              SELECT CASE(coord_sys)
              case('LON-LAT')
                 write(lusrc,20) x,y,z,ix,iy,iz,(src(ic,ix,iy,iz),ic=1,nc)
20               format(2(1x,f11.6),2x,f9.0,2x,3(1x,i5),2x,100(e16.9,1x))
              case('UTM')
                 write(lusrc,21) x,y,z,ix,iy,iz,(src(ic,ix,iy,iz),ic=1,nc)
21               format(2(1x,f10.1),2x,f9.0,2x,3(1x,i5),2x,100(e16.9,1x))
              END SELECT
           end if
        end do
     end do
  end do
  !
  return
  end subroutine wrisrc_mesh
