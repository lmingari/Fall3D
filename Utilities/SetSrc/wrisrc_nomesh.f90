subroutine wrisrc_nomesh(isec1,isec2)
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
  integer(ip) :: is,ic,nsrc
  real(rp)    :: x,y,z,MFR
  !
  !***  Gets the number of sources
  !
  nsrc = ns
  !
  !***  Total mass flow rate
  !
  MFR = SUM(Mplum)
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
  do is = 1,ns
     x = Xplum(is)
     y = Yplum(is)
     if(terrain_following) then
        z = Zplum(is) - Z0 -dZ0  ! terrain following coordinates
     else
        z = Zplum(is)            ! a.s.l.
     end if
     SELECT CASE(coord_sys)
     case('LON-LAT')
        write(lusrc,20) x,y,z,(Mplum(ic,is),ic=1,nc)
20      format(2(1x,f11.6),2x,f9.0,2x,100(e16.9,1x))
     case('UTM')
        write(lusrc,21) x,y,z,(Mplum(ic,is),ic=1,nc)
21      format(2(1x,f10.1),2x,f9.0,2x,100(e16.9,1x))
     END SELECT
  end do
  !
  return
end subroutine wrisrc_nomesh
