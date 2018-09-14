  subroutine wriplumetem(isec1,isec2)
  !******************************************************************************
  !*
  !*     Writes the plume temperature file
  !*
  !******************************************************************************
  use KindType
  use InpOut
  use Master
  implicit none
  integer(ip) :: isec1,isec2
  !
  logical     :: found
  integer(ip) :: nsour,ix,iy,iz,ic,is,js
  real(rp)    :: z,zp,dist,distmin,Mpoint
  !
  !***  Gets the number of sources
  !
  call getnsrc(nsour)
  !
  !***  Writes header for the current interval
  !
  write(lutem,10) isec1,isec2
  write(lutem,11) nsour
10 format(i7,1x,i7)
11 format(i7)
  !
  !***  Loop over all points
  !
  do iz = 1,nz
     do iy = 1,ny
        do ix = 1,nx
           found = .false.
           do ic = 1,nc
              if(src(ic,ix,iy,iz).gt.0.0) found = .true.
           end do
           !
           !***           The point has associated mass; temperature is written
           !
           if(found) then
              z = zmodel(iz)               ! terrain
              !
              !***              Find the plume closest point (vertical coordinate only)
              !
              distmin = 1d20
              do is = 1,ns
                 zp = Zplum(is)-Z0-dZ0    ! terrain
                 dist = abs(zp-z)
                 if(dist.lt.distmin) then
                    js = is
                    distmin = dist
                 end if
              end do
              !
              !***   Mass of water proportional to the total mass
              !***   This is NOT the released mass to be transported (defined in the soruce
              !***   term file) but an estimation of the amount of water mass
              !***   to consider aggregation inside the plume
              !
              Mpoint = 0.
              do ic = 1,npart
                 Mpoint = Mpoint + src(ic,ix,iy,iz)
              end do
              !
              if(ngas.gt.0) then
                 write(lutem,20) ix,iy,iz,Tplum(js),Mpoint*fc(npart+1)
              else
                 write(lutem,20) ix,iy,iz,Tplum(js),Mpoint*w0
              end if
20            format(2(1x,i8),2x,i4,2x,f7.2,2x,e14.6)
              !
           end if
           !
        end do
     end do
  end do
  !
  return
  end subroutine wriplumetem
!
!
!
subroutine wriplumetem_void(isec1,isec2)
  !******************************************************************************
  !*
  !*     Writes the plume temperature file
  !*
  !******************************************************************************
  use KindType
  use InpOut
  use Master
  implicit none
  integer(ip) :: isec1,isec2
  !
  integer(ip) :: nsour
  !
  !***  Gets the number of sources
  !
  nsour=0
  !
  !***  Writes header for the current interval
  !
  write(lutem,10) isec1,isec2
  write(lutem,11) nsour
10 format(i7,1x,i7)
11 format(i7)
  !
  return
end subroutine wriplumetem_void
