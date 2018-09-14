  subroutine wrisrc_resu(isec1,isec2)
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
  !***  Gets the number of sources
  !
  nsrc = 0
  do iy = 1,ny
  do ix = 1,nx
     found = .false.
     do ic = 1,nc
        if(emis(ix,iy,ic).gt.0.0_rp) found = .true.
     end do
     if(found) nsrc = nsrc + 1
  end do
  end do
  nsrc = nz_emission*nsrc  ! vertical source distribution in nz_emission layers
  !
  !***  MFR
  !
  MFR = SUM(emis)
  !
  !***  Writes header for the current interval
  !
  if(nsrc.gt.0) then
     write(lusrc,10) isec1,isec2
     write(lusrc,11) nsrc,nc
     write(lusrc,12) MFR
10   format(i7,1x,i7)
11   format(i7,1x,i7)
12   format(e16.9)
    !
    !***  Writes the rest of file
    !
    do iy = 1,ny
    do ix = 1,nx
       found = .false.
       do ic = 1,nc
          if(emis(ix,iy,ic).gt.0d0) found = .true.
       end do
       if(found) then
          x = lon(ix)
          y = lat(iy)
          do iz = 1,nz_emission
             z = zmodel(iz)
             write(lusrc,20) x,y,z,ix,iy,iz,(emis(ix,iy,ic)/nz_emission,ic=1,nc)
20          format(2(1x,f11.6),2x,f9.0,2x,3(1x,i5),2x,100(e16.9,1x))
          end do
       end if
    end do
    end do
  !
  else
     write(lusrc,10) isec1,isec2
     write(lusrc,11) 1,nc
     write(lusrc,12) 0.0
     write(lusrc,20) xorigr,yorigr,0.0,1,1,1,(0.0,ic=1,nc)
  end if
  !
  return
  end subroutine wrisrc_resu
