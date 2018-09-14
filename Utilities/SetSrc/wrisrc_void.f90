  subroutine wrisrc_void(isec1,isec2)
  !***************************************************************
  !*
  !*    Writes the source file during the time between the eruption
  !*    end and the end of the run (i.e. with MFR=0)
  !*
  !***************************************************************
  use KindType
  use InpOut
  use Master
  implicit none
  !
  integer(ip) :: isec1,isec2
  integer(ip) :: nsrc,ic
  real   (rp) :: MFR
  !
  !***  Gets the number of sources
  !
  nsrc = 1
  !
  !***  Total mass flow rate
  !
  MFR = 0.0
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
  SELECT CASE(coord_sys)
  case('LON-LAT')
     write(lusrc,20) xorigr,yorigr,0.0,1,1,1,(MFR,ic=1,nc)
20   format(2(1x,f11.6),2x,f9.0,2x,3(1x,i5),2x,100(e16.9,1x))
  case('UTM')
     write(lusrc,21) xorigr,yorigr,0.0,1,1,1,(MFR,ic=1,nc)
21   format(2(1x,f10.1),2x,f9.0,2x,3(1x,i5),2x,100(e16.9,1x))
  END SELECT
  !
  return
  end subroutine wrisrc_void
