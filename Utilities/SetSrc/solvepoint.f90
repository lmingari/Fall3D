  subroutine solvepoint
  !********************************************************
  !*
  !*    Point source solution.
  !*     x : LON
  !*     y : LAT
  !*     z : height a.s.l
  !*
  !********************************************************
  use Master
  implicit none
  !
  !***  Source position
  !
  Xplum(1) = X0
  Yplum(1) = Y0
  Zplum(1) = Z0 + dZ0 + Hplume
  !
  !***  MFR
  !
  Mplum(1:nc,1) = M0*fc(1:nc)
  !
  return
  end subroutine solvepoint
