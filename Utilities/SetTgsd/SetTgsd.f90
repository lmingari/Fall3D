program SetTgsd
  !*****************************************************************************
  !*
  !*    AUTHOR : A.Folch
  !*    VERSION: FALL3D-7.0
  !*    DATE   : DEC-2012
  !*    PURPOSE: This program generates a granulometric distribution in the
  !*             Fall3d/Hazmap input file format
  !*
  !*****************************************************************************
  use InpOut
  implicit none
  !
  !***  Gets filenames from the call arguments
  !
  iarg = 1                          ! log   file
  call GETARG(iarg,lulogname)
  iarg = 2                          ! input file
  call GETARG(iarg,luinpname)
  iarg = 3                          ! granulometry file
  call GETARG(iarg,lugrnname)
  !
  !***  Opens the log file
  !
  call openinp
  !
  !***  Reads input data
  !
  call readinp
  !
  !***  Sets the granulometric distribution
  !
  call setfrac
  !
  !***  Writes the file
  !
  call wrigrn
  !
  !***  Ends
  !
  call runend('OK')
end program SetTgsd
